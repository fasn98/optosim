"""The coupled dynamics, and an exact anchor for the part that has one.

State vector
------------
Four dimensionless variables, each read as a fraction except the last::

    P  pluripotency-programme activation   in [0, 1]   (0 = somatic, 1 = fully on)
    I  somatic identity                    in [0, 1]   (1 = fully differentiated)
    A  epigenetic-age marker               in [0, 1]   (1 = old at start, 0 = young)
    R  cumulative risk index               >= 0        (ordinal, not a probability)

Equations
---------
With factor expression ``f(t)`` in ``[0, 1]`` from the protocol::

    dP/dt = k_on * potency * f * (1 - P)  -  k_off * P
    dI/dt = -k_dediff * P * I   +  k_rediff * (1 - P) * (1 - I)
    dA/dt = -k_rejuv * P * A
    dR/dt =  k_risk * ( relu(P - theta_commit) + relu(theta_identity - I) )

Why this shape, and no finer. ``P`` is a driven leaky switch: it rises toward a
flux-dependent steady state under induction and relaxes when induction stops,
with ``k_off < k_on`` giving the switch its hysteresis. ``potency`` multiplies
the drive, so an enhanced factor variant is represented as the same protocol
acting harder -- see the parameter's own note for why a published fold-change in
marker expression cannot be entered here as a value for it. ``I`` is eroded by active
``P`` and restores when ``P`` is low -- the recovery that cyclic induction
exploits. ``A`` falls only while ``P`` is active and cannot go below zero (the
``A`` factor makes 0 a fixed point), so rejuvenation is bought with the very
variable that threatens identity. ``R`` simply accumulates time spent past either
threshold; it is a comparator between protocols, not a hazard rate.

The exact anchor
----------------
Under constant ``f`` the ``P`` equation is linear with constant coefficients, so
it has a closed-form solution. :func:`pluripotency_step_response` returns it, and
the test suite checks the integrator against it to machine precision -- the same
move the parent package makes with the opsin photocycle. It validates the
*integration of P*, and nothing about whether the model is right about biology.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..parameters import ParameterSet
from .parameters import REPROGRAM_CONCEPTUAL


def _relu(x: float) -> float:
    return x if x > 0.0 else 0.0


@dataclass(frozen=True, slots=True)
class ConceptualReprogramming:
    """The four-variable conceptual model. Immutable; parameters from a set."""

    parameters: ParameterSet = REPROGRAM_CONCEPTUAL

    def __post_init__(self) -> None:
        if self.parameters.model != "reprogram_conceptual":
            raise ValueError(
                f"expected a reprogram_conceptual parameter set, got model "
                f"{self.parameters.model!r} from {self.parameters.name!r}"
            )

    @property
    def n_states(self) -> int:
        return 4

    @property
    def state_names(self) -> tuple[str, ...]:
        return ("P", "I", "A", "R")

    @property
    def model_name(self) -> str:
        return f"conceptual reprogramming ({self.parameters.name})"

    def initial_state(self, *, identity: float = 1.0, epi_age: float = 1.0) -> tuple[float, float, float, float]:
        """A somatic cell before induction: no pluripotency, no risk yet.

        Args:
            identity: starting somatic identity, default fully differentiated.
            epi_age: starting epigenetic-age marker, default ``1.0`` (old), so a
                run reports rejuvenation as the fraction of this that is removed.
        """
        if not 0.0 <= identity <= 1.0:
            raise ValueError(f"identity is a fraction in [0, 1], got {identity}")
        if not 0.0 <= epi_age <= 1.0:
            raise ValueError(f"epi_age is a fraction in [0, 1], got {epi_age}")
        return (0.0, identity, epi_age, 0.0)

    def derivatives(
        self, state: tuple[float, float, float, float], expression: float
    ) -> tuple[float, float, float, float]:
        if len(state) != 4:
            raise ValueError(f"model needs 4 state variables, got {len(state)}")
        p, i, a, _r = state
        k_on = self.parameters.value_of("k_on")
        k_off = self.parameters.value_of("k_off")
        k_dediff = self.parameters.value_of("k_dediff")
        k_rediff = self.parameters.value_of("k_rediff")
        k_rejuv = self.parameters.value_of("k_rejuv")
        theta_commit = self.parameters.value_of("theta_commit")
        theta_identity = self.parameters.value_of("theta_identity")
        k_risk = self.parameters.value_of("k_risk")

        potency = self.parameters.value_of("potency")

        dp = k_on * potency * expression * (1.0 - p) - k_off * p
        di = -k_dediff * p * i + k_rediff * (1.0 - p) * (1.0 - i)
        da = -k_rejuv * p * a
        dr = k_risk * (_relu(p - theta_commit) + _relu(theta_identity - i))
        return (dp, di, da, dr)

    def pluripotency_steady_state(self, expression: float) -> float:
        """``P_inf`` under constant ``f``: ``k_on c f / (k_on c f + k_off)``."""
        k_on = self.parameters.value_of("k_on")
        k_off = self.parameters.value_of("k_off")
        drive = k_on * self.parameters.value_of("potency") * expression
        denom = drive + k_off
        return drive / denom if denom > 0.0 else 0.0

    def pluripotency_step_response(
        self, expression: float, p0: float, time_ms: float
    ) -> float:
        """Exact ``P(t)`` under constant ``f`` from initial ``p0``.

        The ``P`` equation is linear with constant coefficients when ``f`` is
        held, so ``P(t) = P_inf + (p0 - P_inf) * exp(-(k_on f + k_off) t)``.
        This is the reference the integrator is validated against.
        """
        k_on = self.parameters.value_of("k_on")
        k_off = self.parameters.value_of("k_off")
        rate = k_on * self.parameters.value_of("potency") * expression + k_off
        p_inf = self.pluripotency_steady_state(expression)
        return p_inf + (p0 - p_inf) * math.exp(-rate * time_ms)
