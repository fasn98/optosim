"""ChR2 as a Markov photocycle, three states with the door open for four.

The model
---------
Three states, closed/open/desensitised, with light driving activation and also
accelerating recovery::

    dC/dt = Gr(phi)*D - Ga(phi)*C
    dO/dt = Ga(phi)*C - Gd*O
    dD/dt = Gd*O       - Gr(phi)*D

    Ga(phi) = k_a * phi^p / (phi^p + phi_m^p)
    Gr(phi) = Gr0 + k_r * phi^q / (phi^q + phi_m^q)
    Gd      = constant

with ``C + O + D == 1`` by construction. The photocurrent carries ChR2's
inward rectification::

    f_v(V) = v1 * (1 - exp(-(V - E)/v0)) / (V - E)
    I      = g0 * O * f_v(V) * (V - E)

This is PyRhO's three-state formulation, implemented independently from its
published equations and parameter table rather than by calling it -- see
:mod:`optosim.parameters` for the citation and ``scripts/validate_opsin.py`` for
what that independence is worth.

Why three, and what four would add
----------------------------------
A correction worth recording, because the first version of this docstring had it
backwards. It claimed three states cannot produce ChR2's characteristic
peak-then-plateau photocurrent and that the transient needs the four-state
model. Measured on this implementation at the PyRhO parameters, that is false:

    1 mW/mm^2   peak O = 0.297 at 12.3 ms, plateau 0.134, ratio 2.22
    10 mW/mm^2  peak O = 0.627 at  4.7 ms, plateau 0.222, ratio 2.83

The peak is there, and the ratio is in the range reported for ChR2. The
mechanism is the slow dark recovery: ``Gr`` is 0.0002/ms in darkness and reaches
only ~0.1/ms under light, against ``Gd`` at 0.104/ms, so the closed state
empties into the open state faster than desensitised channels return to it. O
overshoots, then settles once C and D equilibrate. No second open state is
needed for that.

What four states actually buys is separate light- and dark-adapted branches,
which give bi-exponential off-kinetics and a wavelength dependence this model has
no representation for. Those are the reasons to want it -- not the peak.

:class:`OpsinModel` is the seam. A four-state model implements the same methods
and drops in without the coupling layer changing.

Numerics
--------
During a constant-irradiance step all three rates are constant, so the system is
linear with constant coefficients and has an exact solution by matrix
exponential. :func:`analytic_step_response` provides it. That is what the
integrator is tested against: an exact reference for the part of the problem
this module controls, independent of whether the model itself is right about
ChR2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol

from ..parameters import CHR2_3STATE, ParameterSet

#: Conversion for the one place two conventions meet. PyRhO publishes g0 in pS;
#: this package works in nS so that nS * mV == pA exactly.
PS_PER_NS = 1000.0


class OpsinModel(Protocol):
    """The interface the coupling layer depends on.

    Deliberately narrow, so a four-state or six-state photocycle can replace a
    three-state one without the neuron or the stimulus layer knowing.
    """

    @property
    def n_states(self) -> int: ...

    def initial_state(self) -> tuple[float, ...]: ...

    def derivatives(
        self, state: tuple[float, ...], photon_flux: float
    ) -> tuple[float, ...]: ...

    def open_fraction(self, state: tuple[float, ...]) -> float: ...

    def current_pa(self, state: tuple[float, ...], voltage_mv: float) -> float: ...


@dataclass(frozen=True, slots=True)
class ChR2ThreeState:
    """The three-state photocycle. Immutable; parameters come from a set.

    Args:
        parameters: A :class:`~optosim.parameters.ParameterSet` for the
            ``chr2_3state`` model. Defaults to the PyRhO fit.
    """

    parameters: ParameterSet = CHR2_3STATE

    def __post_init__(self) -> None:
        if self.parameters.model != "chr2_3state":
            raise ValueError(
                f"expected a chr2_3state parameter set, got model "
                f"{self.parameters.model!r} from {self.parameters.name!r}"
            )

    # -- identity --------------------------------------------------------

    @property
    def n_states(self) -> int:
        return 3

    @property
    def state_names(self) -> tuple[str, ...]:
        return ("C", "O", "D")

    @property
    def model_name(self) -> str:
        return f"ChR2 three-state ({self.parameters.name})"

    # -- rates -----------------------------------------------------------

    def activation_rate(self, photon_flux: float) -> float:
        """``Ga(phi)``, per ms. Zero in darkness, saturating at ``k_a``."""
        if photon_flux <= 0.0:
            return 0.0
        k_a = self.parameters.value_of("k_a")
        phi_m = self.parameters.value_of("phi_m")
        p = self.parameters.value_of("p")
        numerator = photon_flux**p
        return k_a * numerator / (numerator + phi_m**p)

    def recovery_rate(self, photon_flux: float) -> float:
        """``Gr(phi)``, per ms. ``Gr0`` in darkness, rising with light.

        Light-sensitive recovery is what makes the model's desensitisation
        depend on illumination history rather than only on time.
        """
        gr0 = self.parameters.value_of("Gr0")
        if photon_flux <= 0.0:
            return gr0
        k_r = self.parameters.value_of("k_r")
        phi_m = self.parameters.value_of("phi_m")
        q = self.parameters.value_of("q")
        numerator = photon_flux**q
        return gr0 + k_r * numerator / (numerator + phi_m**q)

    def decay_rate(self) -> float:
        """``Gd``, per ms. Light independent in this model."""
        return self.parameters.value_of("Gd")

    # -- dynamics --------------------------------------------------------

    def initial_state(self) -> tuple[float, ...]:
        """Fully dark-adapted: everything closed."""
        return (1.0, 0.0, 0.0)

    def derivatives(
        self, state: tuple[float, ...], photon_flux: float
    ) -> tuple[float, ...]:
        if len(state) != 3:
            raise ValueError(f"three-state model needs 3 occupancies, got {len(state)}")
        closed, opened, desensitised = state
        ga = self.activation_rate(photon_flux)
        gr = self.recovery_rate(photon_flux)
        gd = self.decay_rate()
        return (
            gr * desensitised - ga * closed,
            ga * closed - gd * opened,
            gd * opened - gr * desensitised,
        )

    def open_fraction(self, state: tuple[float, ...]) -> float:
        return state[1]

    # -- current ---------------------------------------------------------

    def rectification(self, voltage_mv: float) -> float:
        """``f_v(V)``, dimensionless.

        Singular at ``V == E``, where the limit is 1: expanding
        ``1 - exp(-x/v0)`` to first order gives ``x/v0``, so ``f_v -> v1/v0``.
        Handled explicitly rather than left to produce a division by zero at the
        one voltage a resting-at-reversal cell would sit at.
        """
        e = self.parameters.value_of("E")
        v0 = self.parameters.value_of("v0")
        v1 = self.parameters.value_of("v1")
        driving = voltage_mv - e
        if abs(driving) < 1e-9:
            return v1 / v0
        return v1 * (1.0 - math.exp(-driving / v0)) / driving

    def current_pa(self, state: tuple[float, ...], voltage_mv: float) -> float:
        """Photocurrent in pA. Negative is inward, the HH sign convention.

        ``g0`` is published in pS and converted here, once. ``nS * mV == pA``
        exactly, which is the identity the whole unit choice rests on.
        """
        g0_ns = self.parameters.value_of("g0") / PS_PER_NS
        e = self.parameters.value_of("E")
        return g0_ns * self.open_fraction(state) * self.rectification(voltage_mv) * (voltage_mv - e)

    # -- closed forms, for validation ------------------------------------

    def steady_state_open_fraction(self, photon_flux: float) -> float | None:
        """Open fraction as ``t -> inf`` at constant flux, in closed form.

        From PyRhO's documented steady state, ``I_SS`` proportional to
        ``Ga*Gr / (Gd*(Gr + Ga) + Ga*Gr)``. Derivable directly: setting all three
        derivatives to zero with ``C + O + D == 1`` gives this ratio.

        Returns:
            The open fraction, or ``None`` in darkness, where ``Ga == 0`` makes
            the open state unreachable and the steady state is simply the dark
            state -- an answer, but not one this ratio computes.
        """
        ga = self.activation_rate(photon_flux)
        if ga <= 0.0:
            return None
        gr = self.recovery_rate(photon_flux)
        gd = self.decay_rate()
        denominator = gd * (gr + ga) + ga * gr
        if denominator <= 0.0:
            return None
        return ga * gr / denominator


def analytic_step_response(
    model: ChR2ThreeState,
    photon_flux: float,
    times_ms: list[float],
    initial: tuple[float, float, float] | None = None,
) -> list[tuple[float, float, float]]:
    """Exact occupancies under a constant-flux step, by matrix exponential.

    Under constant illumination every rate is constant, so ``dx/dt = A x`` with
    constant ``A`` and the solution is ``x(t) = expm(A t) x0``. This is the
    reference the numerical integrator is checked against: it is exact to
    floating point, and it tests the integration without assuming the model is a
    correct description of ChR2.

    Requires numpy. Kept out of the model class so the class itself stays
    importable in a bare environment.
    """
    import numpy as np

    ga = model.activation_rate(photon_flux)
    gr = model.recovery_rate(photon_flux)
    gd = model.decay_rate()
    matrix = np.array(
        [
            [-ga, 0.0, gr],
            [ga, -gd, 0.0],
            [0.0, gd, -gr],
        ]
    )
    start = np.array(initial if initial is not None else model.initial_state(), dtype=float)

    # Eigendecomposition rather than scipy.linalg.expm, so this needs only numpy.
    values, vectors = np.linalg.eig(matrix)
    coefficients = np.linalg.solve(vectors, start)
    out: list[tuple[float, float, float]] = []
    for t in times_ms:
        evolved = vectors @ (coefficients * np.exp(values * t))
        real = np.real_if_close(evolved, tol=1e6)
        out.append((float(real[0]), float(real[1]), float(real[2])))
    return out
