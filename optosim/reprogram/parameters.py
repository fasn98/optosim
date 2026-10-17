"""The parameters of the conceptual reprogramming model -- every one uncited.

This file is the honest heart of the branch. In the parent package the
parameter report's interesting column is "which of these did anyone measure?",
and the answer is "most of them". Here the answer is **none of them**, and
that is not a defect to be fixed later -- it is the correct description of a
model at this level of abstraction.

There is no published parameter set for a four-variable phenomenological model
of partial reprogramming. The quantities below are not reduced-order summaries
of a fitted network; they are the coefficients of a deliberately simple set of
equations chosen to reproduce one qualitative fact -- the trade-off between
rejuvenation and loss of identity -- and nothing finer. Giving any of them an
``Origin.MEASURED`` with a citation would be laundering a modelling choice into
a result, which is the one thing optosim is built to make impossible.

So every parameter is ``uncited(...)`` with a ``why`` that says what it stands
for and why no source backs the number. :func:`optosim.report.parameter_report`
will list all of them, every run. A reader who wants a *calibrated* model knows
immediately that this is not one.
"""

from __future__ import annotations

from ..parameters import ParameterSet, uncited

#: Shared prefix for the "why uncited" strings, so the reason the whole set is
#: heuristic is stated once and not diluted into eight near-copies.
_CONCEPTUAL = (
    "no source: this is a conceptual phenomenological model with no fitted "
    "parameter set. The number sets the qualitative behaviour, not a measured "
    "rate"
)

REPROGRAM_CONCEPTUAL = ParameterSet(
    name="Conceptual partial-reprogramming dynamics (uncalibrated)",
    model="reprogram_conceptual",
    parameters=(
        uncited(
            "k_on",
            0.08,
            "1/ms",
            "rate the pluripotency programme activates under factor expression",
            f"{_CONCEPTUAL}; it sets how fast induction takes hold relative to "
            "the pulse length, which is itself in uncalibrated model-time",
        ),
        uncited(
            "k_off",
            0.03,
            "1/ms",
            "rate the pluripotency programme relaxes once factors are withdrawn",
            f"{_CONCEPTUAL}; k_off < k_on encodes the well-documented hysteresis "
            "of the reprogramming switch -- it turns on faster than it turns off "
            "-- but the specific ratio is a choice, not a measurement",
        ),
        uncited(
            "k_dediff",
            0.05,
            "1/ms",
            "rate active pluripotency erodes somatic identity",
            f"{_CONCEPTUAL}; this is the coupling that makes strong or long "
            "induction dangerous. Its size relative to k_rejuv is what sets "
            "whether a protocol can rejuvenate before it dedifferentiates",
        ),
        uncited(
            "k_rediff",
            0.02,
            "1/ms",
            "rate somatic identity recovers while pluripotency is low",
            f"{_CONCEPTUAL}; it is what lets identity restore in the gaps of a "
            "cyclic protocol. Set smaller than k_dediff so recovery is slower "
            "than loss, a pessimistic and deliberate choice",
        ),
        uncited(
            "k_rejuv",
            0.04,
            "1/ms",
            "rate the epigenetic-age marker falls while pluripotency is active",
            f"{_CONCEPTUAL}; it represents the rejuvenation effect partial "
            "reprogramming is sought for. That it is driven by the SAME "
            "pluripotency variable that erodes identity is the whole point: the "
            "benefit and the hazard share a cause",
        ),
        uncited(
            "theta_commit",
            0.6,
            "",
            "pluripotency level past which commitment (teratoma) risk accrues",
            f"{_CONCEPTUAL}; a conceptual commitment point. Above it the cell is "
            "treated as crossing toward pluripotency and the risk index climbs. "
            "The real commitment is a network bistability, not a scalar threshold",
        ),
        uncited(
            "theta_identity",
            0.4,
            "",
            "identity level below which loss-of-identity risk accrues",
            f"{_CONCEPTUAL}; the floor under which a cell is treated as having "
            "lost its somatic identity. A single number standing in for a cell "
            "fate decision",
        ),
        uncited(
            "potency",
            1.0,
            "",
            "multiplier on factor drive: 1 = reference factors, >1 = enhanced variants",
            f"{_CONCEPTUAL}. It exists to ask one question -- what does a MORE "
            "POTENT factor (e.g. an engineered SOX2/KLF4 variant) do to the safe "
            "window? Reported gains for such variants (e.g. '>50x' for "
            "RetroSOX/RetroKLF, OpenAI/Retro Biosciences 2025) are fold-changes in "
            "pluripotency-MARKER EXPRESSION in vitro, not a rate constant, and "
            "must NOT be entered here as though they were one",
        ),
        uncited(
            "k_risk",
            0.5,
            "1/ms",
            "rate the dimensionless risk index accumulates past the thresholds",
            f"{_CONCEPTUAL}; it converts time-spent-over-threshold into a risk "
            "index. The index is an ordinal signal for comparing protocols, NOT "
            "a probability of tumour formation and must never be read as one",
        ),
    ),
)
