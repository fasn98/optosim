"""The stimulus: reprogramming-factor expression as a function of time.

This is the direct analogue of :mod:`optosim.light`. There the stimulus is
irradiance at the opsin; here it is the **expression level of the reprogramming
factors** (OSK / OSKM) at the cell, as a dimensionless fraction in ``[0, 1]``
where ``1`` is the maximum the protocol drives to and ``0`` is no induction.

Cyclic induction -- the standard strategy for keeping partial reprogramming on
the safe side of the commitment point -- is a pulse train, exactly as a train of
light pulses is the standard optogenetic drive. :meth:`FactorProtocol.train`
builds one.

Scope, stated plainly because the omissions matter
--------------------------------------------------
Two gaps are deliberate and are not to be inferred from a missing function:

1. **Expression is the value AT THE CELL, uniform across it.** Delivery
   efficiency, transduction heterogeneity between cells, mosaic expression and
   trafficking are not modelled. A real doxycycline dose or AAV titre produces a
   *distribution* of expression across a population; this is one representative
   cell at one expression level. Relating a protocol here to an administered
   dose is a layer that does not exist.

2. **The time axis is UNCALIBRATED.** Times and rates are carried in ``ms`` and
   ``1/ms`` only because those are the units optosim's registry recognises. In
   this conceptual model one time unit is **arbitrary model-time**: it is not a
   millisecond, and the model makes no claim about how many hours, days or cell
   cycles a run represents. Only the *shape* of the trajectories and the
   *ordering* of outcomes are meaningful; the clock is not. Every report says so.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True)
class FactorPulse:
    """A rectangular window of constant factor expression.

    Args:
        onset_ms: When induction starts (uncalibrated model-time; see module).
        duration_ms: How long expression is held. Positive; a zero-length pulse
            is not an induction and is refused rather than silently ignored.
        expression: Factor expression level while on, a fraction in ``[0, 1]``.
    """

    onset_ms: float
    duration_ms: float
    expression: float

    def __post_init__(self) -> None:
        if self.duration_ms <= 0:
            raise ValueError(
                f"pulse duration must be positive, got {self.duration_ms}. A "
                "zero-length induction is refused rather than silently doing "
                "nothing."
            )
        if not 0.0 <= self.expression <= 1.0:
            raise ValueError(
                f"expression is a fraction in [0, 1], got {self.expression}. "
                "This is a normalised level, not a dose or a titre."
            )
        if self.onset_ms < 0:
            raise ValueError("pulse onset cannot be before the start of the run")

    @property
    def offset_ms(self) -> float:
        return self.onset_ms + self.duration_ms

    def contains(self, time_ms: float) -> bool:
        """Half-open ``[onset, offset)`` so abutting pulses do not double up."""
        return self.onset_ms <= time_ms < self.offset_ms

    def describe(self) -> str:
        return (
            f"expression {self.expression:g} from {self.onset_ms:g} to "
            f"{self.offset_ms:g} ({self.duration_ms:g} model-time)"
        )


@dataclass(frozen=True, slots=True)
class FactorProtocol:
    """An induction protocol: a train of factor-expression pulses."""

    pulses: tuple[FactorPulse, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        ordered = sorted(self.pulses, key=lambda p: p.onset_ms)
        for earlier, later in zip(ordered, ordered[1:], strict=False):
            if later.onset_ms < earlier.offset_ms:
                raise ValueError(
                    f"pulses overlap: {earlier.describe()} and {later.describe()}. "
                    "Overlapping inductions would sum to an expression the "
                    "protocol never delivered, so they are refused rather than "
                    "added."
                )

    @classmethod
    def single(
        cls, onset_ms: float, duration_ms: float, expression: float
    ) -> FactorProtocol:
        """One sustained induction. The way to over-reprogramme on purpose."""
        return cls(pulses=(FactorPulse(onset_ms, duration_ms, expression),))

    @classmethod
    def train(
        cls,
        n_pulses: int,
        onset_ms: float,
        duration_ms: float,
        period_ms: float,
        expression: float,
    ) -> FactorProtocol:
        """Cyclic induction: a regular pulse train.

        The strategy for keeping partial reprogramming on the safe side of the
        commitment point -- short inductions with recovery gaps between them, so
        identity restores between cycles while the age marker ratchets down.
        """
        if n_pulses < 1:
            raise ValueError("a cyclic protocol needs at least one pulse")
        if period_ms < duration_ms:
            raise ValueError(
                f"period {period_ms} is shorter than the pulse duration "
                f"{duration_ms}, which would overlap"
            )
        return cls(
            pulses=tuple(
                FactorPulse(onset_ms + index * period_ms, duration_ms, expression)
                for index in range(n_pulses)
            ),
        )

    def expression_at(self, time_ms: float) -> float:
        """Factor expression at a time; ``0.0`` when no pulse is active."""
        for pulse in self.pulses:
            if pulse.contains(time_ms):
                return pulse.expression
        return 0.0

    @property
    def switch_times_ms(self) -> tuple[float, ...]:
        """Every instant expression changes.

        The integrator is told these and lands on them exactly, for the same
        reason the parent package lands on light-pulse edges: a fixed step that
        straddles an onset smears it.
        """
        edges: set[float] = set()
        for pulse in self.pulses:
            edges.add(pulse.onset_ms)
            edges.add(pulse.offset_ms)
        return tuple(sorted(edges))

    def delivery_note(self) -> str:
        """The caveat every report carries while delivery is unmodelled."""
        return (
            "Expression is specified AT THE CELL and uniform across it. "
            "Delivery, transduction heterogeneity, mosaic expression and "
            "trafficking are NOT modelled, so a protocol here cannot be read as "
            "a dose or an AAV titre. The time axis is UNCALIBRATED model-time "
            "(carried in ms only for unit compatibility): trajectory shape and "
            "outcome ordering are meaningful, the clock is not."
        )

    def describe(self) -> str:
        if not self.pulses:
            return "  no induction: the whole run is uninduced"
        lines = [f"  {len(self.pulses)} induction pulse(s):"]
        lines.extend(f"    {pulse.describe()}" for pulse in self.pulses)
        lines.append(f"  {self.delivery_note()}")
        return "\n".join(lines)
