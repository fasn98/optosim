"""Provenance-carrying quantities.

The single most important module in this package.

PORTED VERBATIM from the chemdisco project (``chemdisco/provenance.py``), whose
discipline this simulator inherits. It is pure Python with no domain content, so
the two copies will drift; when a third project needs it, this should become a
shared dependency rather than a third copy. Changes here should be mirrored
there, and vice versa, until that happens.

Why a simulator needs this more than a data pipeline does. chemdisco mostly
handled measurements and had to avoid laundering a heuristic into one. Here the
situation is inverted: **essentially every output of this package is
PREDICTED**, because that is what a simulator produces. A membrane voltage trace
is not a recording. A spike time is not an observation. The risk is not that a
guess gets mistaken for a measurement by accident -- it is that a whole
simulated trace, plotted and tabulated and compared, quietly acquires the
authority of data.

So the rule here: model *constants* may be MEASURED when they carry a citation,
or HEURISTIC when they do not, and every *output* computed from them is
PREDICTED, naming the model and the parameter set that produced it. A parameter
with no citation is reported as such, never silently promoted.

Design rules enforced here:

1. A quantity whose value is unknown is ``value=None``. There is no fallback,
   no imputed default, no "safe" placeholder. Downstream code renders that as
   "not computed", never as a number.
2. ``Origin`` is mandatory and carries no default. Forgetting it is a
   ``TypeError``, not a silently-measured value.
3. ``Origin.MEASURED`` requires a non-empty ``source``. You cannot claim
   something was measured without saying by whom.
4. Quantities are immutable. A transformation returns a new quantity whose
   origin degrades (never upgrades) along MEASURED -> DERIVED -> PREDICTED ->
   HEURISTIC, so provenance cannot be laundered by passing a value through a
   pipeline.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field, replace
from enum import Enum
from typing import Any


class Origin(Enum):
    """Where a number actually came from.

    The ordering of the members is the trust ordering, strongest first. It is
    used by :func:`weakest` to decide the origin of a derived quantity.
    """

    MEASURED = "measured"
    """A value read from an experimental record or a published table: a fitted
    Hodgkin-Huxley rate constant with its paper, a ChR2 transition rate with the
    table it came from. Requires a source. Note what this does NOT mean here --
    a model constant fitted by its authors to their recordings is MEASURED in the
    sense that it is traceable to published work, not in the sense that anyone
    measured it in the preparation being simulated."""

    DERIVED = "derived"
    """Computed deterministically from measured inputs by a published,
    unambiguous transformation -- unit conversion, IC50 -> pIC50, a counted
    graph property such as rotatable-bond count. Reproducible to the digit."""

    PREDICTED = "predicted"
    """Output of a model: every voltage trace, photocurrent, open-state
    occupancy and spike time this package computes. Must name the model and the
    parameter set. This is the origin of almost everything a simulator emits,
    and labelling it honestly is the point of the package."""

    HEURISTIC = "heuristic"
    """A value with no citation behind it: a plausible conductance chosen to
    make a demonstration run, a rate of the right order of magnitude taken from
    nobody in particular. Legitimate as a starting point and never to be
    presented as sourced. :func:`optosim.report.parameter_report` lists every
    one of these explicitly."""

    @property
    def rank(self) -> int:
        return _ORIGIN_RANK[self]

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


_ORIGIN_RANK: dict[Origin, int] = {
    Origin.MEASURED: 0,
    Origin.DERIVED: 1,
    Origin.PREDICTED: 2,
    Origin.HEURISTIC: 3,
}


def weakest(origins: Iterable[Origin]) -> Origin:
    """Return the least trustworthy origin in ``origins``.

    A quantity derived from a measurement and a prediction is a prediction. A
    chain is only as strong as its weakest link, which is why provenance can
    degrade through a pipeline but never improve.

    Raises:
        ValueError: if ``origins`` is empty -- there is no sensible answer, and
            returning a default here would be exactly the kind of invented
            value this module exists to prevent.
    """
    ranked = sorted(origins, key=lambda o: o.rank)
    if not ranked:
        raise ValueError("weakest() requires at least one origin")
    return ranked[-1]


class ProvenanceError(ValueError):
    """Raised when a quantity would misrepresent where its value came from."""


@dataclass(frozen=True, slots=True)
class Quantity:
    """A number that knows where it came from, or an explicit absence.

    Attributes:
        value: The magnitude, or ``None`` when genuinely unknown. ``None`` is a
            first-class outcome, not an error state.
        unit: Physical unit as a short string (``"nM"``, ``"kcal/mol"``,
            ``"Da"``). ``None`` for dimensionless quantities such as pIC50,
            counts and probabilities.
        origin: Mandatory :class:`Origin`.
        source: Free-text citation precise enough to find the number again: a
            ChEMBL assay id, a PDB code, a model name and version. Required
            for ``MEASURED``.
        uncertainty: One standard deviation in the same unit as ``value``, when
            the producer can state one honestly. ``None`` means "not
            quantified", which is different from zero.
        in_domain: For ``PREDICTED`` quantities, whether the input fell inside
            the model's applicability domain. ``None`` when not assessed.
        notes: Caveats a reader needs: censored relation, aggregated replicates,
            assay type mismatch.
    """

    value: float | None
    unit: str | None
    origin: Origin
    source: str = ""
    uncertainty: float | None = None
    in_domain: bool | None = None
    notes: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not isinstance(self.origin, Origin):
            raise ProvenanceError(
                f"origin must be an Origin, got {type(self.origin).__name__}"
            )
        if self.origin is Origin.MEASURED and not self.source.strip():
            raise ProvenanceError(
                "Origin.MEASURED requires a non-empty source: a measurement "
                "without a citation is not a measurement"
            )
        if self.value is not None:
            if not isinstance(self.value, (int, float)) or isinstance(self.value, bool):
                raise ProvenanceError(
                    f"value must be a real number or None, got {self.value!r}"
                )
            if math.isnan(self.value) or math.isinf(self.value):
                raise ProvenanceError(
                    f"value must be finite; got {self.value!r}. Use None for "
                    "unknown rather than NaN, which propagates silently."
                )
        if self.uncertainty is not None:
            if self.uncertainty < 0:
                raise ProvenanceError("uncertainty cannot be negative")
            if self.value is None:
                raise ProvenanceError(
                    "uncertainty without a value is meaningless"
                )
        if self.in_domain is not None and self.origin is not Origin.PREDICTED:
            raise ProvenanceError(
                "in_domain only applies to Origin.PREDICTED quantities"
            )
        object.__setattr__(self, "notes", tuple(self.notes))

    # -- constructors ----------------------------------------------------

    @classmethod
    def measured(
        cls,
        value: float | None,
        unit: str | None,
        source: str,
        *,
        uncertainty: float | None = None,
        notes: Sequence[str] = (),
    ) -> Quantity:
        return cls(
            value=value,
            unit=unit,
            origin=Origin.MEASURED,
            source=source,
            uncertainty=uncertainty,
            notes=tuple(notes),
        )

    @classmethod
    def derived(
        cls,
        value: float | None,
        unit: str | None,
        source: str,
        *,
        uncertainty: float | None = None,
        notes: Sequence[str] = (),
    ) -> Quantity:
        return cls(
            value=value,
            unit=unit,
            origin=Origin.DERIVED,
            source=source,
            uncertainty=uncertainty,
            notes=tuple(notes),
        )

    @classmethod
    def predicted(
        cls,
        value: float | None,
        unit: str | None,
        model: str,
        *,
        uncertainty: float | None = None,
        in_domain: bool | None = None,
        notes: Sequence[str] = (),
    ) -> Quantity:
        return cls(
            value=value,
            unit=unit,
            origin=Origin.PREDICTED,
            source=model,
            uncertainty=uncertainty,
            in_domain=in_domain,
            notes=tuple(notes),
        )

    @classmethod
    def heuristic(
        cls,
        value: float | None,
        unit: str | None,
        rule: str,
        *,
        notes: Sequence[str] = (),
    ) -> Quantity:
        return cls(
            value=value,
            unit=unit,
            origin=Origin.HEURISTIC,
            source=rule,
            notes=tuple(notes),
        )

    @classmethod
    def unknown(
        cls, unit: str | None, reason: str, *, origin: Origin = Origin.DERIVED
    ) -> Quantity:
        """An explicit, documented absence.

        This is the correct return value when a calculation fails. The
        predecessor returned ``0.0`` or a mid-range default in these cases,
        which is indistinguishable from a real result.
        """
        return cls(
            value=None,
            unit=unit,
            origin=origin,
            source=reason,
            notes=("value unavailable",),
        )

    # -- interrogation ---------------------------------------------------

    @property
    def is_known(self) -> bool:
        return self.value is not None

    @property
    def is_trustworthy_for_ranking(self) -> bool:
        """Whether this quantity should be allowed to order a candidate list.

        A known prediction that fell outside its applicability domain is not;
        an unassessed prediction is not either, because an unchecked domain is
        an unknown domain. Heuristics are allowed to rank only as tie-breakers,
        so they are excluded here deliberately.
        """
        if not self.is_known:
            return False
        if self.origin in (Origin.MEASURED, Origin.DERIVED):
            return True
        if self.origin is Origin.PREDICTED:
            return self.in_domain is True
        return False

    def require(self) -> float:
        """Return the value, or raise rather than substitute a default."""
        if self.value is None:
            raise ProvenanceError(
                f"quantity is unknown ({self.source or 'no reason recorded'})"
            )
        return float(self.value)

    def or_else(self, fallback: float) -> float:
        """Return the value, or ``fallback``.

        Use only where a caller genuinely tolerates a substitute -- sorting a
        display table, say -- never to feed a model. The verbose name is
        deliberate: it should be obvious in review that a default leaked in.
        """
        return float(self.value) if self.value is not None else float(fallback)

    # -- transformation --------------------------------------------------

    def map_value(
        self,
        func: Callable[[float], float],
        *,
        unit: str | None,
        source: str,
        origin: Origin | None = None,
        notes: Sequence[str] = (),
    ) -> Quantity:
        """Apply ``func`` to the value, carrying provenance forward.

        An unknown quantity maps to an unknown quantity: ``func`` is never
        called on ``None`` and no value is conjured. The resulting origin is
        the weaker of this quantity's origin and ``origin``, so a conversion
        cannot upgrade a prediction into a measurement.
        """
        new_origin = self.origin if origin is None else weakest([self.origin, origin])
        combined_notes = tuple(self.notes) + tuple(notes)
        if self.value is None:
            return Quantity(
                value=None,
                unit=unit,
                origin=new_origin if new_origin is not Origin.MEASURED else Origin.DERIVED,
                source=f"{source} <- {self.source}" if self.source else source,
                notes=combined_notes + ("propagated unknown",),
            )
        return Quantity(
            value=func(float(self.value)),
            unit=unit,
            origin=new_origin,
            source=f"{source} <- {self.source}" if self.source else source,
            in_domain=self.in_domain if new_origin is Origin.PREDICTED else None,
            notes=combined_notes,
        )

    def with_notes(self, *notes: str) -> Quantity:
        return replace(self, notes=tuple(self.notes) + notes)

    # -- rendering -------------------------------------------------------

    def label(self, *, digits: int = 3) -> str:
        """Human-readable rendering that always discloses provenance.

        Every UI surface must use this instead of formatting ``value``
        directly, so a reader can never mistake a heuristic for a measurement.
        """
        if self.value is None:
            return "not computed"
        text = f"{self.value:.{digits}f}"
        if self.uncertainty is not None:
            text += f" ± {self.uncertainty:.{digits}f}"
        if self.unit:
            text += f" {self.unit}"
        text += f" [{self.origin.value}"
        if self.origin is Origin.PREDICTED and self.in_domain is False:
            text += ", OUT OF DOMAIN"
        text += "]"
        return text

    def to_dict(self) -> dict[str, Any]:
        """Serialise losslessly, provenance included.

        Exports use this so a downstream consumer -- a journal reviewer, a
        collaborator, a later version of this pipeline -- receives the
        provenance along with the number instead of a bare float.
        """
        return {
            "value": self.value,
            "unit": self.unit,
            "origin": self.origin.value,
            "source": self.source,
            "uncertainty": self.uncertainty,
            "in_domain": self.in_domain,
            "notes": list(self.notes),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> Quantity:
        return cls(
            value=payload.get("value"),
            unit=payload.get("unit"),
            origin=Origin(payload["origin"]),
            source=payload.get("source", ""),
            uncertainty=payload.get("uncertainty"),
            in_domain=payload.get("in_domain"),
            notes=tuple(payload.get("notes", ())),
        )


UNKNOWN_DIMENSIONLESS = Quantity.unknown(None, "not calculated")
