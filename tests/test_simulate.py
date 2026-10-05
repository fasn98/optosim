"""The coupled run, its provenance, and the light protocol driving it."""

from __future__ import annotations

import unittest

from optosim.light import IrradianceProtocol, LightPulse
from optosim.provenance import Origin
from optosim.simulate import Trace, simulate_hh, simulate_lif


class TestLightProtocol(unittest.TestCase):
    def test_a_zero_length_pulse_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            LightPulse(10.0, 0.0, 1.0)

    def test_negative_irradiance_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            LightPulse(10.0, 5.0, -1.0)

    def test_overlapping_pulses_are_refused_not_summed(self) -> None:
        # Summing them would deliver an irradiance no source in the protocol
        # actually produced.
        with self.assertRaises(ValueError) as context:
            IrradianceProtocol(pulses=(LightPulse(0.0, 10.0, 1.0), LightPulse(5.0, 10.0, 1.0)))
        self.assertIn("overlap", str(context.exception))

    def test_abutting_pulses_are_allowed(self) -> None:
        # Half-open intervals, so [0,10) and [10,20) do not collide.
        protocol = IrradianceProtocol(
            pulses=(LightPulse(0.0, 10.0, 1.0), LightPulse(10.0, 10.0, 2.0))
        )
        self.assertEqual(protocol.irradiance_at(9.999), 1.0)
        self.assertEqual(protocol.irradiance_at(10.0), 2.0)

    def test_darkness_outside_every_pulse(self) -> None:
        protocol = IrradianceProtocol.single(10.0, 5.0, 3.0)
        self.assertEqual(protocol.irradiance_at(0.0), 0.0)
        self.assertEqual(protocol.irradiance_at(15.0), 0.0)
        self.assertEqual(protocol.irradiance_at(12.0), 3.0)

    def test_switch_times_are_exposed_for_the_integrator(self) -> None:
        # A light step is a discontinuity. An integrator that straddles one
        # smears the onset, so it has to be told where they are.
        protocol = IrradianceProtocol.train(3, 10.0, 5.0, 20.0, 1.0)
        self.assertEqual(protocol.switch_times_ms, (10.0, 15.0, 30.0, 35.0, 50.0, 55.0))

    def test_a_train_with_too_short_a_period_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            IrradianceProtocol.train(3, 0.0, 10.0, 5.0, 1.0)

    def test_the_depth_caveat_is_always_available(self) -> None:
        # The omission has to be stated, not inferred from a missing function.
        note = IrradianceProtocol.single(0.0, 1.0, 1.0).depth_note()
        self.assertIn("not modelled", note.lower())
        self.assertIn("AT THE OPSIN", note)


class TestTraceProvenance(unittest.TestCase):
    def test_a_trace_must_have_matching_times_and_values(self) -> None:
        with self.assertRaises(ValueError):
            Trace("V", "mV", (0.0, 1.0), (-65.0,), "m", ("p",))

    def test_summaries_are_predicted_and_name_the_model(self) -> None:
        trace = Trace("V", "mV", (0.0, 1.0), (-65.0, -40.0), "HH", ("squid",))
        peak = trace.peak()
        self.assertIs(peak.origin, Origin.PREDICTED)
        self.assertEqual(peak.value, -40.0)
        self.assertIn("HH", peak.source)
        self.assertIn("squid", peak.source)

    def test_an_empty_trace_summarises_as_unknown_not_zero(self) -> None:
        trace = Trace("V", "mV", (), (), "HH", ("squid",))
        self.assertFalse(trace.peak().is_known)
        self.assertIsNone(trace.peak().value)

    def test_value_at_a_time_picks_the_nearest_sample(self) -> None:
        trace = Trace("V", "mV", (0.0, 1.0, 2.0), (-65.0, -50.0, -20.0), "HH", ("squid",))
        self.assertEqual(trace.at(1.1).value, -50.0)
        self.assertEqual(trace.at(1.9).value, -20.0)


class TestCoupledRun(unittest.TestCase):
    def test_darkness_produces_no_photocurrent_and_no_spikes(self) -> None:
        result = simulate_hh(IrradianceProtocol(), duration_ms=50.0, dt_ms=0.02)
        self.assertEqual(result.spike_times_ms, ())
        for current in result.photocurrent_pa.values:
            self.assertEqual(current, 0.0)
        for opened in result.open_fraction.values:
            self.assertEqual(opened, 0.0)

    def test_light_opens_the_channel_and_depolarises(self) -> None:
        result = simulate_hh(
            IrradianceProtocol.single(10.0, 30.0, 1.0), duration_ms=60.0, dt_ms=0.02
        )
        self.assertGreater(max(result.open_fraction.values), 0.0)
        self.assertLess(min(result.photocurrent_pa.values), 0.0)
        self.assertGreater(max(result.voltage.values), -65.0)

    def test_the_photocurrent_is_inward_while_the_cell_is_below_reversal(self) -> None:
        result = simulate_hh(
            IrradianceProtocol.single(10.0, 30.0, 0.1), duration_ms=60.0, dt_ms=0.02
        )
        for voltage, current in zip(
            result.voltage.values, result.photocurrent_pa.values, strict=True
        ):
            if voltage < -5.0 and current != 0.0:
                self.assertLess(current, 0.0)

    def test_occupancies_stay_physical_through_a_coupled_run(self) -> None:
        result = simulate_hh(
            IrradianceProtocol.train(4, 10.0, 5.0, 20.0, 2.0), duration_ms=100.0, dt_ms=0.02
        )
        for opened in result.open_fraction.values:
            self.assertGreaterEqual(opened, -1e-9)
            self.assertLessEqual(opened, 1.0 + 1e-9)

    def test_the_grid_lands_exactly_on_every_light_switch(self) -> None:
        # Otherwise the onset is smeared across a step and the measured latency
        # is an artefact of dt.
        protocol = IrradianceProtocol.train(2, 10.0, 5.0, 20.0, 1.0)
        result = simulate_hh(protocol, duration_ms=60.0, dt_ms=0.03)
        for switch in protocol.switch_times_ms:
            self.assertIn(switch, result.voltage.times_ms, f"grid misses {switch} ms")

    def test_the_adaptive_integrator_agrees_with_the_fixed_step(self) -> None:
        # Not a biology test: it confirms the default step is small enough, which
        # is the only thing the adaptive path is for.
        protocol = IrradianceProtocol.single(10.0, 20.0, 1.0)
        fixed = simulate_hh(protocol, duration_ms=50.0, dt_ms=0.005)
        adaptive = simulate_hh(protocol, duration_ms=50.0, dt_ms=0.005, adaptive=True)
        self.assertEqual(len(fixed.voltage), len(adaptive.voltage))
        worst = max(
            abs(a - b)
            for a, b in zip(fixed.voltage.values, adaptive.voltage.values, strict=True)
        )
        self.assertLess(worst, 0.5, f"fixed and adaptive differ by {worst} mV")

    def test_the_integrator_is_named_in_the_result(self) -> None:
        protocol = IrradianceProtocol.single(5.0, 10.0, 1.0)
        self.assertIn("fixed-step", simulate_hh(protocol, duration_ms=20.0).integrator)
        self.assertIn(
            "adaptive", simulate_hh(protocol, duration_ms=20.0, adaptive=True).integrator
        )

    def test_a_nonsense_timestep_is_refused(self) -> None:
        for bad in (0.0, -0.01):
            with self.assertRaises(ValueError):
                simulate_hh(IrradianceProtocol(), duration_ms=10.0, dt_ms=bad)

    def test_a_nonsense_duration_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            simulate_hh(IrradianceProtocol(), duration_ms=0.0, dt_ms=0.01)


class TestBaselineComparison(unittest.TestCase):
    """The pairing exists so the gap can be read."""

    def test_the_baseline_overestimates_spiking_under_sustained_light(self) -> None:
        # HH settles at a depolarised equilibrium and stops firing; LIF has no
        # delayed rectifier and no sodium inactivation, so it integrates and
        # resets indefinitely. Measured: 1 against 39 on this protocol.
        protocol = IrradianceProtocol.single(50.0, 100.0, 5.0)
        hh = simulate_hh(protocol, duration_ms=200.0, dt_ms=0.01)
        lif = simulate_lif(protocol, duration_ms=200.0, dt_ms=0.01)
        self.assertGreater(len(lif.spike_times_ms), len(hh.spike_times_ms))
        self.assertLessEqual(len(hh.spike_times_ms), 3)

    def test_the_hh_spike_count_is_not_monotonic_in_irradiance(self) -> None:
        # More light is not more spikes: the stronger the sustained current, the
        # more depolarised the equilibrium and the more sodium is inactivated.
        counts = [
            len(
                simulate_hh(
                    IrradianceProtocol.single(50.0, 100.0, irradiance),
                    duration_ms=200.0,
                    dt_ms=0.01,
                ).spike_times_ms
            )
            for irradiance in (0.1, 0.3, 1.0, 5.0)
        ]
        self.assertNotEqual(counts, sorted(counts), f"counts were {counts}")

    def test_the_comparison_report_refuses_to_declare_a_winner(self) -> None:
        from optosim.report import model_comparison_report

        protocol = IrradianceProtocol.single(50.0, 100.0, 5.0)
        text = model_comparison_report(
            simulate_hh(protocol, duration_ms=200.0, dt_ms=0.02),
            simulate_lif(protocol, duration_ms=200.0, dt_ms=0.02),
        )
        self.assertIn("does NOT say which model is closer", text)
        self.assertIn("uncited", text)

    def test_the_baseline_refractory_ceiling_is_reported(self) -> None:
        from optosim.neuron.iaf import LeakyIntegrateAndFire

        protocol = IrradianceProtocol.single(10.0, 100.0, 10.0)
        lif = simulate_lif(protocol, duration_ms=120.0, dt_ms=0.01)
        ceiling = LeakyIntegrateAndFire().maximum_firing_rate_hz
        self.assertLessEqual(lif.mean_rate_hz().require(), ceiling + 1e-9)


class TestReport(unittest.TestCase):
    def test_the_parameter_report_lists_every_uncited_parameter(self) -> None:
        from optosim.parameters import DEFAULT_COMPARTMENT, LIF_BASELINE
        from optosim.report import parameter_report

        text = parameter_report(LIF_BASELINE, DEFAULT_COMPARTMENT)
        self.assertIn("UNCITED", text)
        for name in ("V_threshold", "V_reset", "t_refractory", "area"):
            self.assertIn(name, text)

    def test_a_fully_cited_set_says_what_cited_does_not_mean(self) -> None:
        from optosim.parameters import CHR2_3STATE
        from optosim.report import parameter_report

        text = parameter_report(CHR2_3STATE)
        self.assertIn("All 11 parameters carry a citation", text)
        self.assertIn("not measured in the cell being simulated", text)

    def test_the_unmodelled_report_names_tissue_optics_and_temperature(self) -> None:
        from optosim.report import unmodelled_report

        text = unmodelled_report()
        for expected in ("Tissue optics", "Temperature", "Spatial structure", "expression"):
            self.assertIn(expected, text)

    def test_the_run_report_states_that_nothing_is_a_measurement(self) -> None:
        from optosim.report import run_report

        text = run_report(
            simulate_hh(IrradianceProtocol.single(5.0, 10.0, 1.0), duration_ms=30.0, dt_ms=0.05)
        )
        self.assertIn("not a recording", text)
        self.assertIn("no measurement", text)


class TestAreaProvenanceAndSensitivity(unittest.TestCase):
    """The membrane area: derived from capacitance, and its influence measured.

    It was the package's most consequential uncited parameter, because
    photocurrent density is exactly linear in 1/area. It is now derived from a
    whole-cell capacitance and a specific capacitance, both cited. The value did
    not change -- 1e-4 cm^2 -- which is the point: the justification changed, not
    the number.
    """

    def test_the_area_is_derived_from_a_capacitance_and_cited(self) -> None:
        from optosim.parameters import DEFAULT_COMPARTMENT

        area = DEFAULT_COMPARTMENT["area"]
        self.assertTrue(area.is_cited, "the area must no longer be HEURISTIC")
        self.assertIn("C_whole / C_m", area.description)
        self.assertIn("pF", area.quantity.source)
        self.assertEqual(DEFAULT_COMPARTMENT.uncited, ())

    def test_the_derivation_arithmetic_is_right(self) -> None:
        # 100 pF / (1 uF/cm^2) = 1e-10 F / 1e-6 F/cm^2 = 1e-4 cm^2.
        from optosim.parameters import area_from_capacitance

        self.assertAlmostEqual(area_from_capacitance(100.0), 1e-4, places=12)
        self.assertAlmostEqual(area_from_capacitance(250.0), 2.5e-4, places=12)
        self.assertAlmostEqual(area_from_capacitance(50.0, 2.0), 2.5e-5, places=12)

    def test_a_nonsense_capacitance_is_refused(self) -> None:
        from optosim.parameters import area_from_capacitance

        for bad in (0.0, -100.0):
            with self.assertRaises(ValueError):
                area_from_capacitance(bad)
        with self.assertRaises(ValueError):
            area_from_capacitance(100.0, 0.0)

    def test_a_user_supplied_capacitance_carries_its_own_source(self) -> None:
        from optosim.parameters import compartment_from_capacitance

        compartment = compartment_from_capacitance(250.0, source="cell 7, 2026-02-01")
        self.assertAlmostEqual(compartment.value_of("area"), 2.5e-4, places=12)
        self.assertIn("cell 7", compartment["area"].quantity.source)
        self.assertTrue(compartment["area"].is_cited)

    def test_the_density_is_exactly_linear_in_one_over_area(self) -> None:
        # The part that genuinely is linear, and the reason the area matters.
        from optosim.units import opsin_current_density

        one = opsin_current_density(1000.0, 1e-4)
        half = opsin_current_density(1000.0, 5e-5)
        assert one is not None and half is not None
        self.assertAlmostEqual(half / one, 2.0, places=12)

    def test_near_threshold_the_spike_count_depends_on_the_assumed_area(self) -> None:
        # The honest warning, as a test. At low irradiance a 16x area range moves
        # the count from several spikes to none, so the result is a statement
        # about the assumed cell size.
        from optosim.report import area_sensitivity

        counts = [
            row[2]
            for row in area_sensitivity(
                IrradianceProtocol.single(50.0, 100.0, 0.1),
                duration_ms=200.0,
                dt_ms=0.02,
            )
        ]
        self.assertGreater(max(counts), 1)
        self.assertEqual(min(counts), 0)

    def test_under_saturating_light_the_count_is_robust_to_the_area(self) -> None:
        # The other half: the sensitivity is regime dependent, so neither
        # "it matters" nor "it doesn't" is true on its own.
        from optosim.report import area_sensitivity

        counts = {
            row[2]
            for row in area_sensitivity(
                IrradianceProtocol.single(50.0, 100.0, 5.0),
                duration_ms=200.0,
                dt_ms=0.02,
            )
        }
        self.assertEqual(len(counts), 1, f"expected one count across areas, got {counts}")

    def test_every_run_report_states_the_geometry_sensitivity(self) -> None:
        from optosim.report import run_report
        from optosim.simulate import simulate_hh

        text = run_report(
            simulate_hh(IrradianceProtocol.single(5.0, 10.0, 1.0), duration_ms=30.0, dt_ms=0.05)
        )
        self.assertIn("Geometry sensitivity", text)
        self.assertIn("not even monotonic", text)
        self.assertIn("C_whole / C_m", text)
