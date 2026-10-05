"""Hodgkin-Huxley, single compartment.

Classic four-variable formulation: ``V`` with sodium activation ``m``,
sodium inactivation ``h`` and potassium activation ``n``::

    C_m dV/dt = -g_Na m^3 h (V - E_Na) - g_K n^4 (V - E_K) - g_L (V - E_L) + I_ext
    dx/dt     = alpha_x(V) (1 - x) - beta_x(V) x        for x in {m, h, n}

Currents are densities in uA/cm^2 and conductances in mS/cm^2, which is the
convention the source parameters are published in. ``I_ext`` is a density too,
so anything injected -- including the opsin photocurrent, which arrives as an
absolute pA -- must be converted first. :func:`optosim.units.opsin_current_density`
is where that happens and it requires the membrane area to do it.

Two rate expressions are singular
---------------------------------
``alpha_n`` at V = -55 mV and ``alpha_m`` at V = -40 mV are both 0/0: the
numerator and the ``1 - exp(-x/10)`` denominator vanish together. The limits are
exact -- expanding the exponential to first order gives ``0.01x/(x/10) = 0.1``
and ``0.1x/(x/10) = 1.0`` -- and they are handled explicitly. These are not
obscure voltages: -55 mV is a few mV above rest and a depolarising stimulus
passes straight through it, so an unguarded implementation divides by zero on an
ordinary simulation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..parameters import HH_SQUID_AXON, ParameterSet

#: Voltages where alpha_m and alpha_n are 0/0, and the width of the window where
#: the series expansion is used instead. 1e-6 mV is far below any physical
#: resolution and well above float noise.
_SINGULARITY_WINDOW_MV = 1e-6


def alpha_m(voltage_mv: float) -> float:
    """Sodium activation forward rate, per ms. Singular at -40 mV, limit 1.0."""
    x = voltage_mv + 40.0
    if abs(x) < _SINGULARITY_WINDOW_MV:
        return 1.0
    return 0.1 * x / (1.0 - math.exp(-x / 10.0))


def beta_m(voltage_mv: float) -> float:
    """Sodium activation backward rate, per ms."""
    return 4.0 * math.exp(-(voltage_mv + 65.0) / 18.0)


def alpha_h(voltage_mv: float) -> float:
    """Sodium inactivation forward rate, per ms."""
    return 0.07 * math.exp(-(voltage_mv + 65.0) / 20.0)


def beta_h(voltage_mv: float) -> float:
    """Sodium inactivation backward rate, per ms."""
    return 1.0 / (1.0 + math.exp(-(voltage_mv + 35.0) / 10.0))


def alpha_n(voltage_mv: float) -> float:
    """Potassium activation forward rate, per ms. Singular at -55 mV, limit 0.1."""
    x = voltage_mv + 55.0
    if abs(x) < _SINGULARITY_WINDOW_MV:
        return 0.1
    return 0.01 * x / (1.0 - math.exp(-x / 10.0))


def beta_n(voltage_mv: float) -> float:
    """Potassium activation backward rate, per ms."""
    return 0.125 * math.exp(-(voltage_mv + 65.0) / 80.0)


def steady_state_gate(alpha: float, beta: float) -> float | None:
    """``alpha / (alpha + beta)``, or ``None`` if both rates vanish.

    Both vanishing means the gate has no defined equilibrium at that voltage.
    Returning 0.5, or 0.0, would be inventing one.
    """
    total = alpha + beta
    if total <= 0.0:
        return None
    return alpha / total


@dataclass(frozen=True, slots=True)
class HodgkinHuxley:
    """Single-compartment HH membrane.

    Args:
        parameters: An ``hh`` parameter set; defaults to the squid axon values.
        spike_threshold_mv: Voltage an upward crossing must pass to be counted
            as a spike. HH has no threshold parameter -- spiking is a dynamical
            property, not a comparison -- so this exists only for *counting*
            spikes after the fact and is deliberately not part of the model.
    """

    parameters: ParameterSet = HH_SQUID_AXON
    spike_threshold_mv: float = 0.0

    def __post_init__(self) -> None:
        if self.parameters.model != "hh":
            raise ValueError(
                f"expected an hh parameter set, got model {self.parameters.model!r}"
            )

    @property
    def model_name(self) -> str:
        return f"Hodgkin-Huxley ({self.parameters.name})"

    @property
    def state_names(self) -> tuple[str, ...]:
        return ("V", "m", "h", "n")

    def resting_state(self) -> tuple[float, float, float, float]:
        """Gates at their steady values for the resting voltage.

        Starting from anything else produces a large transient in the first
        milliseconds that is an artefact of initialisation, not a response to a
        stimulus -- and which a reader would otherwise mistake for one.
        """
        v = self.parameters.value_of("V_rest")
        gates = []
        for alpha, beta in (
            (alpha_m(v), beta_m(v)),
            (alpha_h(v), beta_h(v)),
            (alpha_n(v), beta_n(v)),
        ):
            value = steady_state_gate(alpha, beta)
            if value is None:
                raise ValueError(
                    f"no steady state for a gate at V = {v} mV; the parameter "
                    "set cannot be rested"
                )
            gates.append(value)
        return (v, gates[0], gates[1], gates[2])

    def ionic_current_density(
        self, state: tuple[float, float, float, float]
    ) -> tuple[float, float, float]:
        """``(I_Na, I_K, I_L)`` in uA/cm^2, outward positive."""
        v, m, h, n = state
        i_na = self.parameters.value_of("g_Na") * (m**3) * h * (v - self.parameters.value_of("E_Na"))
        i_k = self.parameters.value_of("g_K") * (n**4) * (v - self.parameters.value_of("E_K"))
        i_l = self.parameters.value_of("g_L") * (v - self.parameters.value_of("E_L"))
        return (i_na, i_k, i_l)

    def derivatives(
        self,
        state: tuple[float, float, float, float],
        external_current_density: float = 0.0,
    ) -> tuple[float, float, float, float]:
        """``(dV/dt, dm/dt, dh/dt, dn/dt)``.

        Args:
            external_current_density: uA/cm^2, inward negative. The photocurrent
                arrives here after conversion from absolute pA.
        """
        v, m, h, n = state
        i_na, i_k, i_l = self.ionic_current_density(state)
        dv = (-(i_na + i_k + i_l) - external_current_density) / self.parameters.value_of("C_m")
        dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
        dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h
        dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
        return (dv, dm, dh, dn)
