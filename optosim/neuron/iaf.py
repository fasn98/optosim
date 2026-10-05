"""Leaky integrate-and-fire: the pessimistic baseline.

Paired with Hodgkin-Huxley for the reason chemdisco pairs a scaffold split with
a random split. Not because it is better -- it is strictly less -- but because
the gap between the two says how much of a result is coming from the biophysical
detail and how much from the setup. A stimulus protocol that drives both models
to the same spike count has not demonstrated anything about channel dynamics.

The model::

    C_m dV/dt = -g_L (V - E_L) + I_ext        while V < V_threshold
    on crossing:  emit a spike, V <- V_reset, hold for t_refractory

What it cannot do, stated because the comparison is only honest if the baseline's
limits are explicit:

* No spike shape. The "spike" is a recorded time and a reset, so any measure
  that depends on the waveform is meaningless here.
* No sodium inactivation, so no refractoriness beyond the imposed dead time and
  no accommodation to slow depolarisation.
* The threshold is a free parameter with no source (see
  :data:`optosim.parameters.LIF_BASELINE`). Spike *counts* from this model are
  therefore comparisons against an arbitrary choice, which is why the gap
  between the models is informative but its sign is not, on its own, a result.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..parameters import LIF_BASELINE, ParameterSet


@dataclass(frozen=True, slots=True)
class LeakyIntegrateAndFire:
    """Single-compartment LIF, in the same uA/cm^2 convention as the HH model."""

    parameters: ParameterSet = LIF_BASELINE

    def __post_init__(self) -> None:
        if self.parameters.model != "lif":
            raise ValueError(
                f"expected a lif parameter set, got model {self.parameters.model!r}"
            )

    @property
    def model_name(self) -> str:
        return f"Leaky integrate-and-fire ({self.parameters.name})"

    @property
    def state_names(self) -> tuple[str, ...]:
        return ("V",)

    def resting_state(self) -> float:
        return self.parameters.value_of("E_L")

    def derivative(self, voltage_mv: float, external_current_density: float = 0.0) -> float:
        """``dV/dt`` in mV/ms. ``external_current_density`` is uA/cm^2, inward negative."""
        leak = self.parameters.value_of("g_L") * (voltage_mv - self.parameters.value_of("E_L"))
        return (-leak - external_current_density) / self.parameters.value_of("C_m")

    def has_crossed_threshold(self, voltage_mv: float) -> bool:
        return voltage_mv >= self.parameters.value_of("V_threshold")

    def reset_voltage(self) -> float:
        return self.parameters.value_of("V_reset")

    def refractory_ms(self) -> float:
        return self.parameters.value_of("t_refractory")

    @property
    def maximum_firing_rate_hz(self) -> float:
        """The rate ceiling the refractory period imposes, in Hz.

        An artefact of an uncited parameter, exposed so it can be reported
        rather than discovered in a saturating f-I curve.
        """
        refractory = self.refractory_ms()
        return 1000.0 / refractory if refractory > 0 else float("inf")
