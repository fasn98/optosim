"""Integrate the conceptual reprogramming model and return results that know what they are.

Structurally this is the parent package's :mod:`optosim.simulate` with the
membrane swapped for a regulatory state and light swapped for factor
expression. The deliberate reuse is the point: every series that leaves here is
an :class:`optosim.simulate.Trace`, which can only be ``PREDICTED``, so a
reprogramming trajectory cannot be handed onward looking like a measurement any
more than a voltage trace can.

Two things are inherited exactly because the reasons carry over unchanged:

- **The grid lands on every switch time.** Induction onsets and withdrawals are
  discontinuities in ``f(t)``, so a step that straddles one smears it.
- **Fixed-step RK4 by default.** Reproducible to the sample, and there is no
  stochasticity anywhere (``tests/test_no_fabrication.py`` fails the build if a
  module reaches an RNG), so two identical runs agree bit for bit.

What is *not* inherited is any claim about the clock. The parent package's
milliseconds are milliseconds. Here they are uncalibrated model-time, and
:meth:`ReprogramResult.summary` says so on every run.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..parameters import ParameterSet
from ..provenance import Origin, Quantity
from ..simulate import Trace
from .factors import FactorProtocol
from .model import ConceptualReprogramming


@dataclass(frozen=True, slots=True)
class ReprogramResult:
    """One run: four trajectories and the protocol that produced them."""

    pluripotency: Trace
    identity: Trace
    epigenetic_age: Trace
    risk_index: Trace
    model: str
    parameter_sets: tuple[ParameterSet, ...]
    protocol: FactorProtocol
    integrator: str

    @property
    def traces(self) -> tuple[Trace, ...]:
        return (self.pluripotency, self.identity, self.epigenetic_age, self.risk_index)

    def _predicted(self, value: float | None, unit: str, what: str) -> Quantity:
        if value is None:
            return Quantity.unknown(unit, f"{what} of an empty run", origin=Origin.PREDICTED)
        return Quantity.predicted(
            value,
            unit,
            self.model,
            notes=(
                f"{what}, from the simulated trajectory",
                "UNCALIBRATED model-time: shape and ordering are meaningful, the clock is not",
            ),
        )

    def rejuvenation(self) -> Quantity:
        """Fraction of the starting epigenetic-age marker removed by the run.

        ``0`` is no rejuvenation, ``1`` would be the marker driven to zero. A
        fraction of a conceptual marker -- not years, not a methylation clock
        readout, and not comparable to a published epigenetic age.
        """
        values = self.epigenetic_age.values
        if not values:
            return self._predicted(None, "", "rejuvenation")
        start = values[0]
        if start <= 0.0:
            return Quantity.unknown(
                "", "started at zero epigenetic age; no fraction to remove",
                origin=Origin.PREDICTED,
            )
        return self._predicted((start - values[-1]) / start, "", "fraction of the age marker removed")

    def identity_floor(self) -> Quantity:
        """Lowest somatic identity reached. The cost side of the trade-off."""
        values = self.identity.values
        return self._predicted(min(values) if values else None, "", "minimum somatic identity")

    def identity_recovered(self) -> Quantity:
        """Identity at the end of the run: whether the cell came back."""
        values = self.identity.values
        return self._predicted(values[-1] if values else None, "", "final somatic identity")

    def peak_pluripotency(self) -> Quantity:
        """Highest pluripotency activation reached."""
        values = self.pluripotency.values
        return self._predicted(max(values) if values else None, "", "peak pluripotency activation")

    def final_risk(self) -> Quantity:
        """Accumulated risk index.

        **An ordinal comparator between protocols, not a probability.** It is
        time spent past a conceptual threshold, times an uncited rate. Nothing
        converts it into a tumour incidence, and no run here licenses a
        statement about how likely a teratoma is.
        """
        values = self.risk_index.values
        quantity = self._predicted(values[-1] if values else None, "", "accumulated risk index")
        return quantity.with_notes(
            "ORDINAL ONLY: this index is not a probability, a hazard rate or an incidence"
        )


def _sample_grid(duration_ms: float, dt_ms: float, protocol: FactorProtocol) -> list[float]:
    """A uniform grid with every induction switch time inserted exactly."""
    if dt_ms <= 0:
        raise ValueError(f"dt must be positive, got {dt_ms}")
    if duration_ms <= 0:
        raise ValueError(f"duration must be positive, got {duration_ms}")
    n = int(round(duration_ms / dt_ms))
    grid = {round(index * dt_ms, 9) for index in range(n + 1)}
    grid.update(t for t in protocol.switch_times_ms if 0.0 <= t <= duration_ms)
    return sorted(grid)


def _rk4_step(derivatives, state: tuple[float, ...], dt: float, expression: float) -> tuple[float, ...]:
    """One classical RK4 step at constant factor expression.

    Expression is held across the step rather than evaluated mid-step. Since the
    grid lands on every switch time, within a step it really is constant, so
    holding it is exact rather than an approximation.
    """
    k1 = derivatives(state, expression)
    half = tuple(s + 0.5 * dt * d for s, d in zip(state, k1, strict=True))
    k2 = derivatives(half, expression)
    half2 = tuple(s + 0.5 * dt * d for s, d in zip(state, k2, strict=True))
    k3 = derivatives(half2, expression)
    full = tuple(s + dt * d for s, d in zip(state, k3, strict=True))
    k4 = derivatives(full, expression)
    return tuple(
        s + (dt / 6.0) * (a + 2.0 * b + 2.0 * c + d)
        for s, a, b, c, d in zip(state, k1, k2, k3, k4, strict=True)
    )


def simulate_reprogramming(
    protocol: FactorProtocol,
    *,
    duration_ms: float = 400.0,
    dt_ms: float = 0.05,
    model: ConceptualReprogramming | None = None,
    identity: float = 1.0,
    epi_age: float = 1.0,
) -> ReprogramResult:
    """Integrate the four-variable conceptual model under an induction protocol.

    Args:
        protocol: The induction protocol. An empty one is a legitimate run --
            the uninduced control, which should show no rejuvenation and no risk.
        duration_ms: Run length in uncalibrated model-time.
        dt_ms: Fixed step, in the same uncalibrated units.
        identity: Starting somatic identity.
        epi_age: Starting epigenetic-age marker.

    Returns:
        A :class:`ReprogramResult` whose every trace is ``PREDICTED``.
    """
    model = model or ConceptualReprogramming()
    state = model.initial_state(identity=identity, epi_age=epi_age)
    grid = _sample_grid(duration_ms, dt_ms, protocol)

    times: list[float] = [grid[0]]
    series: tuple[list[float], ...] = ([state[0]], [state[1]], [state[2]], [state[3]])

    for start, end in zip(grid, grid[1:], strict=False):
        step = end - start
        if step <= 0:
            continue
        state = _rk4_step(model.derivatives, state, step, protocol.expression_at(start))
        times.append(end)
        for index, value in enumerate(state):
            series[index].append(value)

    names = (model.parameters.name,)
    note = (
        protocol.delivery_note(),
        "conceptual model, no fitted parameters: every constant behind this is UNCITED",
    )
    frozen_times = tuple(times)

    def trace(name: str, values: list[float]) -> Trace:
        return Trace(
            name=name,
            unit="",
            times_ms=frozen_times,
            values=tuple(values),
            model=model.model_name,
            parameter_sets=names,
            notes=note,
        )

    return ReprogramResult(
        pluripotency=trace("P", series[0]),
        identity=trace("I", series[1]),
        epigenetic_age=trace("A", series[2]),
        risk_index=trace("R", series[3]),
        model=model.model_name,
        parameter_sets=(model.parameters,),
        protocol=protocol,
        integrator=f"fixed-step RK4, dt {dt_ms} model-time units",
    )
