"""What a reprogramming run is standing on, and what it must not be read as.

The parent package's report answers "which of these constants did anyone
measure?". Here that question has a one-word answer -- none -- so this module's
work is different: it names the **regime** a protocol landed in, and it states
the limits that make the regime the only honest output.

The regime, not the number, is the result. The time axis is uncalibrated and
every rate is a modelling choice, so "rejuvenation 0.41" means nothing on its
own. What survives those caveats is the *ordering*: that a cyclic protocol
reaches a given rejuvenation at lower risk than a sustained one, and that past
some drive the cell stops coming back. Those are properties of the equations'
shape, and they are what the model is for.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..parameters import ParameterSet
from .model import ConceptualReprogramming
from .simulate import ReprogramResult

#: The independent things that can go wrong, each with the threshold it tests.
#: They are reported together rather than collapsed into one ordered label, and
#: that is a correction rather than a style choice: an earlier version returned
#: a single worst-wins label, which printed "identity lost" for a protocol that
#: was *better* on both identity and risk than one it printed "committed" for.
#: A severity ranking over conditions that are not mutually exclusive hides the
#: comparison the model exists to make.
HAZARDS = (
    "committed toward pluripotency",
    "identity lost",
    "identity not recovered",
)


def regime_flags(
    result: ReprogramResult, *, model: ConceptualReprogramming | None = None
) -> tuple[str, ...]:
    """Every hazard condition this run met, in :data:`HAZARDS` order.

    Each is a trajectory read against one of the model's own conceptual
    thresholds, and inherits that threshold's arbitrariness exactly. An empty
    tuple means no condition was met, which is not the same as "safe".
    """
    model = model or ConceptualReprogramming()
    peak_p = result.peak_pluripotency().value
    floor_i = result.identity_floor().value
    final_i = result.identity_recovered().value
    if peak_p is None or floor_i is None or final_i is None:
        return ()

    flags: list[str] = []
    if peak_p > model.parameters.value_of("theta_commit"):
        flags.append(HAZARDS[0])
    if floor_i < model.parameters.value_of("theta_identity"):
        flags.append(HAZARDS[1])
    if final_i < 0.9 * result.identity.values[0]:
        flags.append(HAZARDS[2])
    return tuple(flags)


def classify_regime(result: ReprogramResult, *, model: ConceptualReprogramming | None = None) -> str:
    """A one-line regime label: every hazard met, or the clean outcome.

    A label for comparing protocols within this model -- not a clinical
    category and not a prediction about a cell.
    """
    peak_p = result.peak_pluripotency().value
    rejuv = result.rejuvenation().value
    if peak_p is None or rejuv is None:
        return "not classifiable: the run produced no trajectory"
    if peak_p <= 0.0:
        return "uninduced"
    flags = regime_flags(result, model=model)
    if flags:
        return "; ".join(flags)
    return "rejuvenating, no threshold crossed"


@dataclass(frozen=True, slots=True)
class ReprogramReport:
    """A printable account of one run: outcome, regime, and every limit on it."""

    result: ReprogramResult

    @property
    def parameter_sets(self) -> tuple[ParameterSet, ...]:
        return self.result.parameter_sets

    def outcome(self) -> str:
        """The four numbers a run produces, each with its provenance label."""
        rows = (
            ("rejuvenation (fraction of age marker removed)", self.result.rejuvenation()),
            ("peak pluripotency activation", self.result.peak_pluripotency()),
            ("minimum somatic identity", self.result.identity_floor()),
            ("final somatic identity", self.result.identity_recovered()),
            ("accumulated risk index (ORDINAL)", self.result.final_risk()),
        )
        lines = ["=" * 78, "Outcome", "=" * 78]
        lines.extend(f"  {label:<46} {quantity.label(digits=4)}" for label, quantity in rows)
        lines.append(f"\n  regime: {classify_regime(self.result)}")
        lines.append(f"  integrator: {self.result.integrator}")
        return "\n".join(lines)

    def caveats(self) -> str:
        """Everything that must be read before any number above is quoted.

        Deliberately not a footnote. In the parent package the caveats narrow
        what a trace can be compared with; here they are most of what the run
        establishes, so they are printed at the same weight as the result.
        """
        return "\n".join(
            [
                "=" * 78,
                "What this run is NOT",
                "=" * 78,
                "  - NOT a validated model. It is not fitted to any dataset and has never",
                "    been compared with a reprogramming experiment. There is no published",
                "    parameter set for a model at this level of abstraction, so every",
                "    constant behind these numbers is UNCITED and listed as such.",
                "",
                "  - NOT a gene-regulatory network. Four lumped variables stand in for the",
                "    OSK/OSKM programme, somatic identity, an epigenetic-age marker and a",
                "    risk index. No gene, no transcription factor, no methylation site and",
                "    no chromatin state is represented.",
                "",
                "  - NOT calibrated in time. The clock is arbitrary model-time carried in",
                "    'ms' only for unit compatibility. The model makes no claim about hours,",
                "    days or cell cycles, so a 'duration' here maps to no protocol.",
                "",
                "  - NOT a risk estimate. The risk index is accumulated time past two",
                "    conceptual thresholds, times an uncited rate. It orders protocols",
                "    within this model. It is not a probability of tumour formation, not a",
                "    hazard rate, and not an incidence.",
                "",
                "  - NOT an epigenetic clock. The age marker is a dimensionless fraction. It",
                "    is not a methylation clock readout and not convertible into years.",
                "",
                "  - NOT a dose. Expression is specified at the cell and uniform across it.",
                "    Delivery, transduction heterogeneity, mosaic expression and trafficking",
                "    are absent, and one representative cell stands in for a population.",
                "",
                "  - NOT a safety claim about any therapy, factor variant or trial. Published",
                "    fold-changes in marker expression are in-vitro measurements of a",
                "    different quantity and cannot be entered as a rate in these equations.",
                "",
                "  What it IS: a visualisation of the trade-off the equations encode -- that",
                "  the variable which lowers the age marker is the same one that erodes",
                "  identity, so rejuvenation and commitment risk share a cause. The regime",
                "  label and the ORDERING of protocols are the output. The digits are not.",
            ]
        )

    def unmodelled(self) -> str:
        """Mechanisms absent from the model, listed so they are not inferred as handled."""
        return "\n".join(
            [
                "=" * 78,
                "Not modelled at all",
                "=" * 78,
                "  - Cell cycle, proliferation and clonal expansion. A reprogrammed cell",
                "    that divides is how a damaged clone grows out; there is no division here.",
                "  - Pre-existing oncogenic mutations. The documented interaction in which",
                "    reprogramming accelerates cancer in an already-mutated cell has no",
                "    representation, so this model CANNOT see that pathway.",
                "  - Cell-type dependence. Liver and intestine are reported as the tissues",
                "    where in vivo reprogramming goes wrong first; one generic cell here.",
                "  - Immune response, senescence, paracrine signalling and tissue context.",
                "  - Network bistability. Commitment is a scalar threshold crossing in this",
                "    model, where in reality it is a self-sustaining network state.",
                "  - Stochasticity. Deliberately none; two identical runs agree bit for bit.",
            ]
        )

    def full(self) -> str:
        """Outcome, caveats and omissions together. The only form fit to quote."""
        return "\n\n".join([self.outcome(), self.caveats(), self.unmodelled()])
