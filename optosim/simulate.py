"""Couple an opsin to a membrane, integrate, and return results that know what they are.

Every array this module returns is a :class:`Trace`, which is PREDICTED and
names the models and parameter sets that produced it. There is no path by which
a bare list of floats leaves here, because a bare list of floats is what gets
plotted next to a recording and read as one.

The joint system
----------------
The opsin photocycle does not depend on voltage in the three-state model, but the
photocurrent does, so the two are integrated as one state vector rather than in
sequence::

    HH:  (V, m, h, n, C, O, D)
    LIF: (V, C, O, D)

Light steps are discontinuities. Both integrators are told the switch times from
the protocol and land exactly on them: a fixed step that straddles an onset
smears it over a step, and an adaptive controller that straddles one wastes
rejected steps discovering it.

The LIF reset is not integrable, so it is applied between steps -- which is the
standard treatment and is also why LIF spike times are quantised to the step in a
way HH spike times are not.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .light import IrradianceProtocol
from .neuron.hh import HodgkinHuxley
from .neuron.iaf import LeakyIntegrateAndFire
from .opsin.chr2 import ChR2ThreeState
from .parameters import DEFAULT_COMPARTMENT, ParameterSet
from .provenance import Origin, Quantity
from .units import opsin_current_density


@dataclass(frozen=True, slots=True)
class Trace:
    """A time series that carries its provenance.

    Attributes:
        name: What it is, e.g. ``"V"``.
        unit: Physical unit, checked against the registry.
        times_ms: Sample times.
        values: Values at those times, same length.
        model: The model that produced it, named.
        parameter_sets: Names of every parameter set involved.
        origin: PREDICTED, always. A simulated series is not a measurement and
            there is no constructor here that claims otherwise.
        notes: Caveats that travel with the data.
    """

    name: str
    unit: str
    times_ms: tuple[float, ...]
    values: tuple[float, ...]
    model: str
    parameter_sets: tuple[str, ...]
    origin: Origin = Origin.PREDICTED
    notes: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if len(self.times_ms) != len(self.values):
            raise ValueError(
                f"{self.name}: {len(self.times_ms)} times against "
                f"{len(self.values)} values"
            )
        if self.origin is not Origin.PREDICTED:
            raise ValueError(
                f"{self.name}: a simulated trace is PREDICTED. Nothing in this "
                "package may label a computed series as measured."
            )

    def __len__(self) -> int:
        return len(self.values)

    def _summary(self, value: float | None, what: str) -> Quantity:
        if value is None:
            return Quantity.unknown(self.unit, f"{what} of an empty {self.name} trace",
                                    origin=Origin.PREDICTED)
        return Quantity.predicted(
            value,
            self.unit,
            f"{self.model} [{', '.join(self.parameter_sets)}]",
            notes=self.notes + (f"{what} of the simulated {self.name} trace",),
        )

    def peak(self) -> Quantity:
        """Most positive value, as a PREDICTED quantity."""
        return self._summary(max(self.values) if self.values else None, "maximum")

    def trough(self) -> Quantity:
        """Most negative value, as a PREDICTED quantity."""
        return self._summary(min(self.values) if self.values else None, "minimum")

    def at(self, time_ms: float) -> Quantity:
        """Value at the sample nearest ``time_ms``."""
        if not self.times_ms:
            return self._summary(None, "value")
        index = min(
            range(len(self.times_ms)), key=lambda i: abs(self.times_ms[i] - time_ms)
        )
        return self._summary(self.values[index], f"value at t={self.times_ms[index]:g} ms")


@dataclass(frozen=True, slots=True)
class SimulationResult:
    """Everything one run produced, each piece carrying its provenance."""

    voltage: Trace
    open_fraction: Trace
    photocurrent_pa: Trace
    spike_times_ms: tuple[float, ...]
    model: str
    parameter_sets: tuple[ParameterSet, ...]
    protocol: IrradianceProtocol
    integrator: str
    stopped_early: str = ""

    @property
    def n_spikes(self) -> Quantity:
        """Spike count as a PREDICTED quantity.

        A count, not an observation. The notes say which threshold produced it,
        because for HH the threshold is a counting convention applied after the
        fact and for LIF it is the model's own free parameter -- two different
        things with the same name.
        """
        return Quantity.predicted(
            float(len(self.spike_times_ms)),
            None,
            self.model,
            notes=(
                "simulated spike count, not a recording",
                f"integrator: {self.integrator}",
            ),
        )

    def mean_rate_hz(self) -> Quantity:
        """Mean firing rate over the simulated window."""
        if not self.voltage.times_ms:
            return Quantity.unknown("", "no trace", origin=Origin.PREDICTED)
        span_ms = self.voltage.times_ms[-1] - self.voltage.times_ms[0]
        if span_ms <= 0:
            return Quantity.unknown("", "zero-length run", origin=Origin.PREDICTED)
        return Quantity.predicted(
            len(self.spike_times_ms) * 1000.0 / span_ms,
            "",
            self.model,
            notes=("mean over the whole window, including any pre-stimulus darkness",),
        )


def _sample_grid(duration_ms: float, dt_ms: float, protocol: IrradianceProtocol) -> list[float]:
    """A uniform grid with every light switch time inserted exactly."""
    if dt_ms <= 0:
        raise ValueError(f"dt must be positive, got {dt_ms}")
    if duration_ms <= 0:
        raise ValueError(f"duration must be positive, got {duration_ms}")
    n = int(round(duration_ms / dt_ms))
    grid = {round(index * dt_ms, 9) for index in range(n + 1)}
    grid.update(t for t in protocol.switch_times_ms if 0.0 <= t <= duration_ms)
    return sorted(grid)


def _rk4_step(
    derivatives, state: tuple[float, ...], dt: float, flux_at_start: float
) -> tuple[float, ...]:
    """One classical RK4 step at constant flux.

    Flux is held at its value for the step rather than evaluated mid-step. The
    grid lands on every switch time, so within a step the irradiance really is
    constant and holding it is exact rather than an approximation.
    """
    k1 = derivatives(state, flux_at_start)
    half = tuple(s + 0.5 * dt * d for s, d in zip(state, k1, strict=True))
    k2 = derivatives(half, flux_at_start)
    half2 = tuple(s + 0.5 * dt * d for s, d in zip(state, k2, strict=True))
    k3 = derivatives(half2, flux_at_start)
    full = tuple(s + dt * d for s, d in zip(state, k3, strict=True))
    k4 = derivatives(full, flux_at_start)
    return tuple(
        s + (dt / 6.0) * (a + 2.0 * b + 2.0 * c + d)
        for s, a, b, c, d in zip(state, k1, k2, k3, k4, strict=True)
    )


def simulate_hh(
    protocol: IrradianceProtocol,
    *,
    duration_ms: float = 200.0,
    dt_ms: float = 0.01,
    neuron: HodgkinHuxley | None = None,
    opsin: ChR2ThreeState | None = None,
    geometry: ParameterSet = DEFAULT_COMPARTMENT,
    adaptive: bool = False,
    tolerance_mv: float = 1e-3,
    max_dt_ms: float = 0.1,
) -> SimulationResult:
    """Integrate the coupled ChR2 + Hodgkin-Huxley system.

    Args:
        adaptive: Use step-doubling error control instead of a fixed step. The
            fixed step is the default because it is reproducible to the sample;
            adaptive is there to confirm the fixed step is small enough, which is
            a question about the integration and not about the biology.
        tolerance_mv: Per-step voltage error target for the adaptive integrator.
    """
    neuron = neuron or HodgkinHuxley()
    opsin = opsin or ChR2ThreeState()
    area_cm2 = geometry.value_of("area")

    def derivatives(state: tuple[float, ...], flux: float) -> tuple[float, ...]:
        v, m, h, n, closed, opened, desens = state
        current_pa = opsin.current_pa((closed, opened, desens), v)
        density = opsin_current_density(current_pa, area_cm2)
        if density is None:
            raise ValueError(
                f"membrane area {area_cm2} cm^2 cannot convert a photocurrent; "
                "refusing rather than dropping the opsin current silently"
            )
        dv, dm, dh, dn = neuron.derivatives((v, m, h, n), density)
        dc, do, dd = opsin.derivatives((closed, opened, desens), flux)
        return (dv, dm, dh, dn, dc, do, dd)

    v0, m0, h0, n0 = neuron.resting_state()
    c0, o0, d0 = opsin.initial_state()
    state: tuple[float, ...] = (v0, m0, h0, n0, c0, o0, d0)

    times, voltages, opens, currents, spikes = _run(
        derivatives,
        state,
        protocol,
        duration_ms,
        dt_ms,
        adaptive=adaptive,
        tolerance=tolerance_mv,
        max_dt_ms=max_dt_ms,
        voltage_index=0,
        open_index=5,
        current=lambda s: opsin.current_pa((s[4], s[5], s[6]), s[0]),
        spike_threshold=neuron.spike_threshold_mv,
    )

    sets = (neuron.parameters, opsin.parameters, geometry)
    names = tuple(s.name for s in sets)
    model = f"{neuron.model_name} + {opsin.model_name}"
    integrator = (
        f"adaptive RK4 step doubling, tol {tolerance_mv} mV, max dt {max_dt_ms} ms"
        if adaptive
        else f"fixed-step RK4, dt {dt_ms} ms"
    )
    note = (protocol.depth_note(),)
    return SimulationResult(
        voltage=Trace("V", "mV", times, voltages, model, names, notes=note),
        open_fraction=Trace("O", "", times, opens, model, names, notes=note),
        photocurrent_pa=Trace("I_opsin", "pA", times, currents, model, names, notes=note),
        spike_times_ms=spikes,
        model=model,
        parameter_sets=sets,
        protocol=protocol,
        integrator=integrator,
    )


def simulate_lif(
    protocol: IrradianceProtocol,
    *,
    duration_ms: float = 200.0,
    dt_ms: float = 0.01,
    neuron: LeakyIntegrateAndFire | None = None,
    opsin: ChR2ThreeState | None = None,
    geometry: ParameterSet = DEFAULT_COMPARTMENT,
) -> SimulationResult:
    """Integrate the coupled ChR2 + LIF system, the pessimistic baseline.

    Fixed step only. An adaptive controller would be false precision here: the
    reset is applied between steps, so spike times are quantised to ``dt_ms``
    whatever the controller does.
    """
    neuron = neuron or LeakyIntegrateAndFire()
    opsin = opsin or ChR2ThreeState()
    area_cm2 = geometry.value_of("area")

    grid = _sample_grid(duration_ms, dt_ms, protocol)
    voltage = neuron.resting_state()
    opsin_state = opsin.initial_state()
    refractory_until = float("-inf")

    times: list[float] = [grid[0]]
    voltages: list[float] = [voltage]
    opens: list[float] = [opsin_state[1]]
    currents: list[float] = [opsin.current_pa(opsin_state, voltage)]
    spikes: list[float] = []

    for start, end in zip(grid, grid[1:], strict=False):
        dt = end - start
        if dt <= 0:
            continue
        flux = protocol.photon_flux_at(start)

        def joint(state: tuple[float, ...], f: float) -> tuple[float, ...]:
            v, c, o, d = state
            current_pa = opsin.current_pa((c, o, d), v)
            density = opsin_current_density(current_pa, area_cm2)
            if density is None:
                raise ValueError(f"membrane area {area_cm2} cm^2 cannot convert a photocurrent")
            dv = neuron.derivative(v, density)
            dc, do, dd = opsin.derivatives((c, o, d), f)
            return (dv, dc, do, dd)

        advanced = _rk4_step(joint, (voltage, *opsin_state), dt, flux)
        voltage, opsin_state = advanced[0], (advanced[1], advanced[2], advanced[3])

        # The reset is a discrete event, applied between steps. This is the
        # standard treatment and is why LIF spike times are quantised to dt.
        if end >= refractory_until and neuron.has_crossed_threshold(voltage):
            spikes.append(end)
            voltage = neuron.reset_voltage()
            refractory_until = end + neuron.refractory_ms()
        elif end < refractory_until:
            voltage = neuron.reset_voltage()

        times.append(end)
        voltages.append(voltage)
        opens.append(opsin_state[1])
        currents.append(opsin.current_pa(opsin_state, voltage))

    sets = (neuron.parameters, opsin.parameters, geometry)
    names = tuple(s.name for s in sets)
    model = f"{neuron.model_name} + {opsin.model_name}"
    note = (
        protocol.depth_note(),
        "LIF baseline: no spike waveform, no sodium inactivation, and an uncited "
        "threshold. Spike counts are comparisons against an arbitrary choice.",
    )
    return SimulationResult(
        voltage=Trace("V", "mV", tuple(times), tuple(voltages), model, names, notes=note),
        open_fraction=Trace("O", "", tuple(times), tuple(opens), model, names, notes=note),
        photocurrent_pa=Trace("I_opsin", "pA", tuple(times), tuple(currents), model, names, notes=note),
        spike_times_ms=tuple(spikes),
        model=model,
        parameter_sets=sets,
        protocol=protocol,
        integrator=f"fixed-step RK4, dt {dt_ms} ms",
    )


def _run(
    derivatives,
    state: tuple[float, ...],
    protocol: IrradianceProtocol,
    duration_ms: float,
    dt_ms: float,
    *,
    adaptive: bool,
    tolerance: float,
    max_dt_ms: float,
    voltage_index: int,
    open_index: int,
    current,
    spike_threshold: float,
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    """Integrate, recording traces and upward threshold crossings."""
    grid = _sample_grid(duration_ms, dt_ms, protocol)
    times = [grid[0]]
    voltages = [state[voltage_index]]
    opens = [state[open_index]]
    currents = [current(state)]
    spikes: list[float] = []
    previous_voltage = state[voltage_index]

    for start, end in zip(grid, grid[1:], strict=False):
        span = end - start
        if span <= 0:
            continue
        flux = protocol.photon_flux_at(start)
        if not adaptive:
            state = _rk4_step(derivatives, state, span, flux)
        else:
            state = _adaptive_segment(
                derivatives, state, span, flux, tolerance, max_dt_ms, voltage_index
            )
        voltage = state[voltage_index]
        # Upward crossing only. Counting both directions would double every
        # spike, and counting "above threshold" would count a plateau as a train.
        if previous_voltage < spike_threshold <= voltage:
            spikes.append(end)
        previous_voltage = voltage
        times.append(end)
        voltages.append(voltage)
        opens.append(state[open_index])
        currents.append(current(state))

    return tuple(times), tuple(voltages), tuple(opens), tuple(currents), tuple(spikes)


def _adaptive_segment(
    derivatives,
    state: tuple[float, ...],
    span: float,
    flux: float,
    tolerance: float,
    max_dt_ms: float,
    voltage_index: int,
) -> tuple[float, ...]:
    """Cross one grid interval with step-doubling error control.

    One RK4 step of size h against two of h/2 differ by about 15 times the error
    of the finer pair, which is the Richardson estimate for a fourth-order
    method. Only the voltage is used to judge the error: the gating variables and
    occupancies are bounded in [0, 1] and never set the accuracy, while the
    voltage moves 100 mV in a millisecond during a spike.
    """
    remaining = span
    dt = min(span, max_dt_ms)
    while remaining > 1e-12:
        dt = min(dt, remaining)
        coarse = _rk4_step(derivatives, state, dt, flux)
        middle = _rk4_step(derivatives, state, 0.5 * dt, flux)
        fine = _rk4_step(derivatives, middle, 0.5 * dt, flux)
        error = abs(fine[voltage_index] - coarse[voltage_index]) / 15.0
        if error > tolerance and dt > 1e-9:
            dt *= 0.5
            continue
        state = fine
        remaining -= dt
        if error < 0.1 * tolerance:
            dt = min(2.0 * dt, max_dt_ms)
    return state
