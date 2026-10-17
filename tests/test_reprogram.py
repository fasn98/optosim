"""The conceptual reprogramming branch.

Three kinds of test, and the distinction matters more here than anywhere else in
the package:

1. **Integration correctness.** ``P`` under constant drive has a closed form, so
   the integrator is checked against it to machine precision. This is a real
   verification and it verifies only the arithmetic.
2. **Structural invariants.** Fractions stay in range, the risk index never
   decreases, an uninduced run does nothing. These pin the equations' shape.
3. **The qualitative claims the module makes in prose.** If the README says
   cyclic induction reaches equal rejuvenation at lower peak pluripotency, a
   test asserts it, so the claim cannot quietly stop being true.

What is deliberately NOT tested, because no test could establish it: that any of
this describes a cell. Every parameter is uncited and the time axis is
uncalibrated; the tests check the implementation of a conceptual model.
"""

from __future__ import annotations

import unittest

from optosim.parameters import ParameterSet
from optosim.provenance import Origin
from optosim.reprogram import (
    ConceptualReprogramming,
    FactorProtocol,
    FactorPulse,
    ReprogramReport,
    classify_regime,
    regime_flags,
    simulate_reprogramming,
)

RECOVERY = 150.0
ONSET = 10.0


def run_matched(protocol: FactorProtocol, span: float, **kwargs):
    """A run with the same recovery window every protocol gets."""
    return simulate_reprogramming(
        protocol, duration_ms=ONSET + span + RECOVERY, dt_ms=0.05, **kwargs
    )


class TestFactorProtocol(unittest.TestCase):
    def test_a_zero_length_induction_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            FactorPulse(0.0, 0.0, 1.0)

    def test_expression_outside_the_unit_interval_is_refused(self) -> None:
        for bad in (-0.1, 1.1, 50.0):
            with self.assertRaises(ValueError, msg=f"{bad} was accepted"):
                FactorPulse(0.0, 10.0, bad)

    def test_overlapping_pulses_are_refused_rather_than_summed(self) -> None:
        with self.assertRaises(ValueError):
            FactorProtocol(pulses=(FactorPulse(0.0, 10.0, 0.5), FactorPulse(5.0, 10.0, 0.5)))

    def test_abutting_pulses_do_not_double_up(self) -> None:
        protocol = FactorProtocol(
            pulses=(FactorPulse(0.0, 10.0, 0.4), FactorPulse(10.0, 10.0, 0.9))
        )
        self.assertEqual(protocol.expression_at(9.999), 0.4)
        self.assertEqual(protocol.expression_at(10.0), 0.9)

    def test_a_train_with_a_period_shorter_than_the_pulse_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            FactorProtocol.train(3, 0.0, 10.0, 5.0, 1.0)

    def test_switch_times_include_every_edge(self) -> None:
        protocol = FactorProtocol.train(3, 10.0, 5.0, 20.0, 1.0)
        self.assertEqual(protocol.switch_times_ms, (10.0, 15.0, 30.0, 35.0, 50.0, 55.0))

    def test_expression_is_zero_outside_every_pulse(self) -> None:
        protocol = FactorProtocol.single(10.0, 5.0, 1.0)
        self.assertEqual(protocol.expression_at(0.0), 0.0)
        self.assertEqual(protocol.expression_at(15.0), 0.0)

    def test_the_delivery_note_names_what_is_unmodelled(self) -> None:
        note = FactorProtocol().delivery_note()
        self.assertIn("UNCALIBRATED", note)
        for word in ("Delivery", "titre", "uniform"):
            self.assertIn(word, note)


class TestExactStepResponse(unittest.TestCase):
    """The one part of this model with a closed form, and the integrator against it."""

    def test_the_integrator_matches_the_exact_solution_to_machine_precision(self) -> None:
        model = ConceptualReprogramming()
        for expression in (0.25, 0.5, 1.0):
            duration = 200.0
            result = simulate_reprogramming(
                FactorProtocol.single(0.0, duration, expression),
                duration_ms=duration,
                dt_ms=0.05,
            )
            worst = max(
                abs(value - model.pluripotency_step_response(expression, 0.0, time_ms))
                for time_ms, value in zip(
                    result.pluripotency.times_ms, result.pluripotency.values, strict=True
                )
            )
            self.assertLess(
                worst,
                1e-9,
                f"integration of P disagrees with its closed form at f={expression}: {worst:.3e}",
            )

    def test_the_steady_state_formula_matches_a_long_run(self) -> None:
        model = ConceptualReprogramming()
        expression = 1.0
        duration = 2000.0
        result = simulate_reprogramming(
            FactorProtocol.single(0.0, duration, expression), duration_ms=duration, dt_ms=0.1
        )
        self.assertAlmostEqual(
            result.pluripotency.values[-1],
            model.pluripotency_steady_state(expression),
            places=9,
        )

    def test_potency_scales_the_drive_in_the_closed_form_too(self) -> None:
        """A potency change must move the analytic anchor, not just the integration."""
        base = ConceptualReprogramming()
        potent = ConceptualReprogramming(parameters=base.parameters.with_overrides(potency=4.0))
        self.assertGreater(
            potent.pluripotency_steady_state(0.2), base.pluripotency_steady_state(0.2)
        )
        duration = 300.0
        result = simulate_reprogramming(
            FactorProtocol.single(0.0, duration, 0.2),
            duration_ms=duration,
            dt_ms=0.05,
            model=potent,
        )
        worst = max(
            abs(value - potent.pluripotency_step_response(0.2, 0.0, time_ms))
            for time_ms, value in zip(
                result.pluripotency.times_ms, result.pluripotency.values, strict=True
            )
        )
        self.assertLess(worst, 1e-9, f"potency broke the analytic anchor: {worst:.3e}")


class TestStructuralInvariants(unittest.TestCase):
    def test_an_uninduced_run_does_nothing_at_all(self) -> None:
        result = run_matched(FactorProtocol(), 40.0)
        self.assertEqual(result.peak_pluripotency().value, 0.0)
        self.assertEqual(result.rejuvenation().value, 0.0)
        self.assertEqual(result.final_risk().value, 0.0)
        self.assertEqual(result.identity_floor().value, 1.0)
        self.assertEqual(classify_regime(result), "uninduced")

    def test_every_fraction_stays_inside_its_range(self) -> None:
        result = run_matched(FactorProtocol.single(ONSET, 80.0, 1.0), 80.0)
        for trace in (result.pluripotency, result.identity, result.epigenetic_age):
            self.assertGreaterEqual(min(trace.values), -1e-12, f"{trace.name} went negative")
            self.assertLessEqual(max(trace.values), 1.0 + 1e-12, f"{trace.name} exceeded 1")

    def test_the_risk_index_never_decreases(self) -> None:
        """It is an accumulation. A fall would mean risk being paid back."""
        result = run_matched(FactorProtocol.train(6, ONSET, 10.0, 30.0, 1.0), 160.0)
        values = result.risk_index.values
        for earlier, later in zip(values, values[1:], strict=False):
            self.assertGreaterEqual(later, earlier - 1e-12)

    def test_the_age_marker_is_monotonically_non_increasing(self) -> None:
        """Nothing in these equations ages a cell, so A must never rise."""
        result = run_matched(FactorProtocol.train(4, ONSET, 10.0, 30.0, 0.8), 100.0)
        values = result.epigenetic_age.values
        for earlier, later in zip(values, values[1:], strict=False):
            self.assertLessEqual(later, earlier + 1e-12)

    def test_zero_is_a_fixed_point_of_the_age_marker(self) -> None:
        result = simulate_reprogramming(
            FactorProtocol.single(0.0, 100.0, 1.0), duration_ms=100.0, dt_ms=0.05, epi_age=0.0
        )
        self.assertEqual(max(abs(v) for v in result.epigenetic_age.values), 0.0)

    def test_pluripotency_relaxes_after_withdrawal(self) -> None:
        result = run_matched(FactorProtocol.single(ONSET, 40.0, 1.0), 40.0)
        peak = result.peak_pluripotency().value
        self.assertLess(result.pluripotency.values[-1], 0.05 * peak)

    def test_the_switch_turns_on_faster_than_it_turns_off(self) -> None:
        """k_off < k_on is the hysteresis the model claims; assert it holds."""
        parameters = ConceptualReprogramming().parameters
        self.assertLess(parameters.value_of("k_off"), parameters.value_of("k_on"))

    def test_a_longer_induction_removes_more_of_the_age_marker(self) -> None:
        previous = -1.0
        for exposure in (10.0, 20.0, 40.0, 80.0):
            value = run_matched(
                FactorProtocol.single(ONSET, exposure, 1.0), exposure
            ).rejuvenation().value
            self.assertGreater(value, previous)
            previous = value

    def test_the_integration_is_insensitive_to_the_step(self) -> None:
        coarse = simulate_reprogramming(
            FactorProtocol.train(4, ONSET, 10.0, 30.0, 1.0), duration_ms=250.0, dt_ms=0.1
        )
        fine = simulate_reprogramming(
            FactorProtocol.train(4, ONSET, 10.0, 30.0, 1.0), duration_ms=250.0, dt_ms=0.01
        )
        self.assertAlmostEqual(
            coarse.rejuvenation().value, fine.rejuvenation().value, places=6
        )
        self.assertAlmostEqual(
            coarse.final_risk().value, fine.final_risk().value, places=4
        )


class TestTheQualitativeClaims(unittest.TestCase):
    """The trade-off the module exists to show. If these fail, the prose is wrong."""

    def test_stronger_drive_buys_rejuvenation_with_identity(self) -> None:
        rejuvenation: list[float] = []
        floors: list[float] = []
        for expression in (0.1, 0.3, 0.5, 0.75, 1.0):
            result = run_matched(FactorProtocol.single(ONSET, 40.0, expression), 40.0)
            rejuvenation.append(result.rejuvenation().value)
            floors.append(result.identity_floor().value)
        self.assertEqual(rejuvenation, sorted(rejuvenation), "rejuvenation should rise with drive")
        self.assertEqual(
            floors, sorted(floors, reverse=True), "identity floor should fall as drive rises"
        )

    def test_cyclic_induction_beats_sustained_at_matched_exposure(self) -> None:
        """More rejuvenation, lower peak pluripotency, lower risk. The core claim.

        Both arms deliver 40 units of induction and both get the same recovery
        window afterwards. Without matching the recovery the comparison measures
        the schedule's bookkeeping -- an earlier version of the demo script did
        exactly that and the cyclic arm looked worse for no biological reason.
        """
        sustained = run_matched(FactorProtocol.single(ONSET, 40.0, 1.0), 40.0)
        pulse = 40.0 / 4
        period = pulse * 3.0
        span = period * 3 + pulse
        cyclic = run_matched(FactorProtocol.train(4, ONSET, pulse, period, 1.0), span)

        self.assertGreater(cyclic.rejuvenation().value, sustained.rejuvenation().value)
        self.assertLess(
            cyclic.peak_pluripotency().value, sustained.peak_pluripotency().value
        )
        self.assertLess(cyclic.final_risk().value, sustained.final_risk().value)
        self.assertGreater(cyclic.identity_floor().value, sustained.identity_floor().value)

    def test_a_more_potent_factor_can_push_a_benign_protocol_over_the_line(self) -> None:
        """The question the coverage raises: potency is not safety.

        At the reference potency this protocol crosses no threshold. Multiply the
        drive and the same protocol commits. This is a property of the equations,
        not evidence about any engineered factor variant.
        """
        base = ConceptualReprogramming()
        benign = run_matched(FactorProtocol.single(ONSET, 40.0, 0.3), 40.0)
        self.assertEqual(regime_flags(benign), ())
        self.assertEqual(benign.final_risk().value, 0.0)

        potent_model = ConceptualReprogramming(
            parameters=base.parameters.with_overrides(potency=5.0)
        )
        potent = run_matched(FactorProtocol.single(ONSET, 40.0, 0.3), 40.0, model=potent_model)
        self.assertIn("committed toward pluripotency", regime_flags(potent))
        self.assertGreater(potent.final_risk().value, 0.0)

    def test_the_regime_label_reports_every_hazard_not_just_the_worst(self) -> None:
        """A worst-wins label hid identity loss behind commitment; it must not.

        This is a regression test for a real defect: ranking non-exclusive
        conditions by severity printed a *better* protocol as the worse-sounding
        one, because the better arm's commitment flag cleared while its identity
        flag stayed set.
        """
        result = run_matched(FactorProtocol.single(ONSET, 40.0, 1.0), 40.0)
        flags = regime_flags(result)
        self.assertIn("committed toward pluripotency", flags)
        self.assertIn("identity lost", flags)
        label = classify_regime(result)
        for flag in flags:
            self.assertIn(flag, label)


class TestProvenance(unittest.TestCase):
    def test_every_parameter_in_the_set_is_uncited(self) -> None:
        """The honest state of a conceptual model, asserted so it stays honest.

        If someone adds a cited parameter to this set, that is a real change and
        this test should be the thing that makes them justify it.
        """
        parameters = ConceptualReprogramming().parameters
        self.assertEqual(
            len(parameters.uncited),
            len(parameters),
            "a parameter in the conceptual set claims a citation; a model with no "
            "fitted parameter set cannot have one",
        )
        for parameter in parameters.uncited:
            self.assertTrue(parameter.quantity.source.strip())

    def test_every_trace_is_predicted(self) -> None:
        result = run_matched(FactorProtocol.single(ONSET, 20.0, 1.0), 20.0)
        for trace in result.traces:
            self.assertIs(trace.origin, Origin.PREDICTED)

    def test_every_summary_quantity_is_predicted(self) -> None:
        result = run_matched(FactorProtocol.single(ONSET, 20.0, 1.0), 20.0)
        for quantity in (
            result.rejuvenation(),
            result.peak_pluripotency(),
            result.identity_floor(),
            result.identity_recovered(),
            result.final_risk(),
        ):
            self.assertIs(quantity.origin, Origin.PREDICTED)

    def test_the_risk_index_carries_its_ordinal_warning(self) -> None:
        result = run_matched(FactorProtocol.single(ONSET, 20.0, 1.0), 20.0)
        notes = " ".join(result.final_risk().notes)
        self.assertIn("not a probability", notes)

    def test_every_trace_carries_the_uncalibrated_time_warning(self) -> None:
        result = run_matched(FactorProtocol.single(ONSET, 20.0, 1.0), 20.0)
        for trace in result.traces:
            self.assertIn("UNCALIBRATED", " ".join(trace.notes))

    def test_overriding_potency_keeps_the_set_heuristic(self) -> None:
        parameters = ConceptualReprogramming().parameters.with_overrides(potency=3.0)
        self.assertIs(parameters["potency"].quantity.origin, Origin.HEURISTIC)
        self.assertEqual(len(parameters.uncited), len(parameters))

    def test_a_mismatched_parameter_set_is_refused(self) -> None:
        wrong = ParameterSet(name="not this model", model="hodgkin_huxley", parameters=())
        with self.assertRaises(ValueError):
            ConceptualReprogramming(parameters=wrong)

    def test_an_unknown_parameter_name_raises(self) -> None:
        with self.assertRaises(KeyError):
            ConceptualReprogramming().parameters.value_of("k_nonexistent")

    def test_the_caveats_name_what_the_model_cannot_see(self) -> None:
        report = ReprogramReport(run_matched(FactorProtocol.single(ONSET, 20.0, 1.0), 20.0))
        caveats = report.full()
        for phrase in (
            "NOT a validated model",
            "NOT calibrated in time",
            "NOT a risk estimate",
            "NOT an epigenetic clock",
            "Pre-existing oncogenic mutations",
            "Cell cycle",
        ):
            self.assertIn(phrase, caveats, f"the report stopped saying: {phrase}")

    def test_two_identical_runs_are_bit_identical(self) -> None:
        protocol = FactorProtocol.train(3, ONSET, 10.0, 25.0, 0.8)
        first = run_matched(protocol, 60.0)
        second = run_matched(protocol, 60.0)
        for a, b in zip(first.traces, second.traces, strict=True):
            self.assertEqual(a.values, b.values)


class TestRefusals(unittest.TestCase):
    def test_a_non_positive_step_is_refused(self) -> None:
        for bad in (0.0, -0.1):
            with self.assertRaises(ValueError):
                simulate_reprogramming(FactorProtocol(), duration_ms=10.0, dt_ms=bad)

    def test_a_non_positive_duration_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            simulate_reprogramming(FactorProtocol(), duration_ms=0.0, dt_ms=0.05)

    def test_initial_fractions_outside_range_are_refused(self) -> None:
        for kwargs in ({"identity": 1.5}, {"identity": -0.1}, {"epi_age": 2.0}):
            with self.assertRaises(ValueError, msg=f"{kwargs} was accepted"):
                ConceptualReprogramming().initial_state(**kwargs)

    def test_a_state_vector_of_the_wrong_length_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            ConceptualReprogramming().derivatives((0.0, 1.0, 1.0), 1.0)


if __name__ == "__main__":
    unittest.main()
