"""The ChR2 photocycle, and the exact reference its integration is checked against."""

from __future__ import annotations

import math
import unittest

from optosim.opsin.chr2 import ChR2ThreeState, analytic_step_response
from optosim.parameters import CHR2_3STATE, HH_SQUID_AXON
from optosim.units import irradiance_to_photon_flux

SATURATING_FLUX = irradiance_to_photon_flux(10.0)
MODERATE_FLUX = irradiance_to_photon_flux(1.0)


class TestPhotocycleStructure(unittest.TestCase):
    def test_occupancies_are_conserved_by_the_derivatives(self) -> None:
        # C + O + D == 1 is a structural invariant, so the derivatives must sum
        # to zero exactly. If they do not, probability leaks and every current
        # downstream is wrong by a drifting amount.
        model = ChR2ThreeState()
        for state in [(1.0, 0.0, 0.0), (0.3, 0.4, 0.3), (0.0, 0.0, 1.0)]:
            for flux in (0.0, MODERATE_FLUX, SATURATING_FLUX):
                total = sum(model.derivatives(state, flux))
                self.assertAlmostEqual(total, 0.0, places=12, msg=f"{state} at {flux}")

    def test_it_starts_dark_adapted(self) -> None:
        self.assertEqual(ChR2ThreeState().initial_state(), (1.0, 0.0, 0.0))

    def test_a_wrong_parameter_set_is_refused(self) -> None:
        # Handing the opsin the neuron's parameters would otherwise fail later,
        # inside an equation, as a confusing KeyError.
        with self.assertRaises(ValueError) as context:
            ChR2ThreeState(parameters=HH_SQUID_AXON)
        self.assertIn("chr2_3state", str(context.exception))

    def test_a_three_state_model_rejects_a_state_of_the_wrong_size(self) -> None:
        with self.assertRaises(ValueError):
            ChR2ThreeState().derivatives((0.5, 0.5), MODERATE_FLUX)


class TestRates(unittest.TestCase):
    def test_darkness_cannot_open_the_channel(self) -> None:
        self.assertEqual(ChR2ThreeState().activation_rate(0.0), 0.0)

    def test_activation_saturates_at_k_a(self) -> None:
        model = ChR2ThreeState()
        k_a = CHR2_3STATE.value_of("k_a")
        enormous = model.activation_rate(1e30)
        self.assertLess(enormous, k_a)
        self.assertGreater(enormous, 0.99 * k_a)

    def test_activation_is_monotonic_in_light(self) -> None:
        model = ChR2ThreeState()
        rates = [model.activation_rate(f) for f in (0.0, 1e15, 1e16, 1e17, 1e18)]
        self.assertEqual(rates, sorted(rates))

    def test_recovery_in_darkness_is_the_dark_rate(self) -> None:
        model = ChR2ThreeState()
        self.assertAlmostEqual(model.recovery_rate(0.0), CHR2_3STATE.value_of("Gr0"))

    def test_light_accelerates_recovery(self) -> None:
        # The feature that makes desensitisation depend on illumination history
        # rather than only on elapsed time.
        model = ChR2ThreeState()
        self.assertGreater(model.recovery_rate(MODERATE_FLUX), model.recovery_rate(0.0))


class TestCurrent(unittest.TestCase):
    def test_a_closed_channel_passes_no_current(self) -> None:
        self.assertEqual(ChR2ThreeState().current_pa((1.0, 0.0, 0.0), -70.0), 0.0)

    def test_the_current_is_inward_below_the_reversal_potential(self) -> None:
        # ChR2 reverses near 0 mV, so at rest the photocurrent depolarises:
        # negative by the HH convention.
        self.assertLess(ChR2ThreeState().current_pa((0.8, 0.2, 0.0), -70.0), 0.0)

    def test_the_current_reverses_above_the_reversal_potential(self) -> None:
        self.assertGreater(ChR2ThreeState().current_pa((0.8, 0.2, 0.0), 40.0), 0.0)

    def test_rectification_is_finite_at_the_reversal_potential(self) -> None:
        # f_v(V) is 0/0 at V == E. The limit is v1/v0; without the explicit case
        # this is a ZeroDivisionError at exactly the voltage a cell clamped to
        # the reversal potential would sit at.
        model = ChR2ThreeState()
        expected = CHR2_3STATE.value_of("v1") / CHR2_3STATE.value_of("v0")
        self.assertAlmostEqual(model.rectification(0.0), expected, places=10)

    def test_rectification_is_continuous_across_the_singularity(self) -> None:
        model = ChR2ThreeState()
        at = model.rectification(0.0)
        for epsilon in (1e-6, 1e-4, 1e-2):
            self.assertAlmostEqual(model.rectification(epsilon), at, places=3)
            self.assertAlmostEqual(model.rectification(-epsilon), at, places=3)

    def test_the_current_scales_linearly_with_open_fraction(self) -> None:
        model = ChR2ThreeState()
        one = model.current_pa((0.9, 0.1, 0.0), -70.0)
        two = model.current_pa((0.8, 0.2, 0.0), -70.0)
        self.assertAlmostEqual(two / one, 2.0, places=9)


class TestClosedFormAgreesWithTheModel(unittest.TestCase):
    """The validation anchor for the numerics, independent of the biology.

    Under constant illumination all three rates are constant, so the system is
    linear and has an exact matrix-exponential solution. Checking the closed-form
    steady state against that exact evolution tests the algebra; checking the
    integrator against it tests the integration. Neither tests whether the model
    describes real ChR2 -- that is what the PyRhO comparison is for.
    """

    def test_the_steady_state_formula_matches_the_exact_evolution(self) -> None:
        model = ChR2ThreeState()
        closed_form = model.steady_state_open_fraction(MODERATE_FLUX)
        assert closed_form is not None
        # Far enough out that the slowest eigenvalue has decayed: Gr0 is 2e-4/ms,
        # so the slow mode has a ~5 s time constant.
        late = analytic_step_response(model, MODERATE_FLUX, [500_000.0])[0]
        self.assertAlmostEqual(late[1], closed_form, places=6)

    def test_the_exact_solution_conserves_probability(self) -> None:
        model = ChR2ThreeState()
        for state in analytic_step_response(
            model, MODERATE_FLUX, [0.0, 1.0, 10.0, 100.0, 1000.0]
        ):
            self.assertAlmostEqual(sum(state), 1.0, places=10)

    def test_the_exact_solution_starts_where_it_was_told_to(self) -> None:
        model = ChR2ThreeState()
        first = analytic_step_response(model, MODERATE_FLUX, [0.0])[0]
        self.assertAlmostEqual(first[0], 1.0, places=12)
        self.assertAlmostEqual(first[1], 0.0, places=12)

    def test_occupancies_stay_physical_under_the_exact_solution(self) -> None:
        model = ChR2ThreeState()
        for state in analytic_step_response(
            model, SATURATING_FLUX, [0.1 * i for i in range(200)]
        ):
            for occupancy in state:
                self.assertGreaterEqual(occupancy, -1e-12)
                self.assertLessEqual(occupancy, 1.0 + 1e-12)

    def test_darkness_has_no_steady_state_open_fraction(self) -> None:
        # Not zero -- unreachable. Returning 0.0 would be a number where the
        # honest answer is that the ratio does not apply.
        self.assertIsNone(ChR2ThreeState().steady_state_open_fraction(0.0))


class TestPeakThenPlateau(unittest.TestCase):
    """The transient peak, which this model does produce.

    These tests replaced an assertion that the response was monotonic. That
    assertion failed, and it was the test that was wrong: at the PyRhO
    parameters the three-state model peaks and then decays to a plateau, with a
    ratio in the range reported for ChR2. The mechanism is slow dark recovery --
    Gd (0.104/ms) against Gr (0.0002/ms dark, ~0.1/ms lit) -- so C empties into
    O faster than D refills it and O overshoots. A second open state is not
    required for a peak.
    """

    def _step(self, flux: float) -> list[float]:
        trace = analytic_step_response(
            ChR2ThreeState(), flux, [0.05 * i for i in range(20_000)]
        )
        return [state[1] for state in trace]

    def test_the_response_peaks_and_then_declines(self) -> None:
        opens = self._step(SATURATING_FLUX)
        peak = max(opens)
        index = opens.index(peak)
        self.assertGreater(index, 0, "a peak at t=0 would mean no rise at all")
        self.assertLess(index, len(opens) - 1, "a peak at the end is a plateau, not a peak")
        self.assertLess(opens[-1], peak, "the trace must settle below its peak")

    def test_the_peak_to_plateau_ratio_is_in_the_reported_range(self) -> None:
        # ChR2 peak:plateau is usually quoted around 2. This pins the model's
        # behaviour; it is NOT evidence the model is right about ChR2, which is
        # what the PyRhO comparison would establish.
        model = ChR2ThreeState()
        for flux, expected in ((MODERATE_FLUX, 2.22), (SATURATING_FLUX, 2.83)):
            plateau = model.steady_state_open_fraction(flux)
            assert plateau is not None
            ratio = max(self._step(flux)) / plateau
            self.assertAlmostEqual(ratio, expected, places=1, msg=f"at flux {flux:.3e}")
            self.assertGreater(ratio, 1.0)

    def test_brighter_light_peaks_sooner(self) -> None:
        # Higher Ga empties the closed state faster, so the overshoot arrives
        # earlier. A consequence of the mechanism, so worth locking.
        moderate = self._step(MODERATE_FLUX)
        saturating = self._step(SATURATING_FLUX)
        self.assertLess(saturating.index(max(saturating)), moderate.index(max(moderate)))

    def test_the_plateau_is_the_closed_form_steady_state(self) -> None:
        model = ChR2ThreeState()
        plateau = model.steady_state_open_fraction(SATURATING_FLUX)
        assert plateau is not None
        late = analytic_step_response(model, SATURATING_FLUX, [500_000.0])[0][1]
        self.assertAlmostEqual(late, plateau, places=6)

    def test_the_model_names_itself_including_its_parameter_set(self) -> None:
        name = ChR2ThreeState().model_name
        self.assertIn("three-state", name)
        self.assertIn("PyRhO", name)


class TestUnitIdentity(unittest.TestCase):
    def test_nS_times_mV_is_pA(self) -> None:
        # The whole unit choice rests on this being exact.
        from optosim.units import nS_times_mV_is_pA

        self.assertTrue(nS_times_mV_is_pA())
        self.assertAlmostEqual(1e-9 * 1e-3, 1e-12, places=24)

    def test_the_conductance_conversion_is_applied_once(self) -> None:
        # g0 is published in pS; a current computed as if it were nS would be a
        # thousand times too large. Pin the magnitude at a known operating point.
        model = ChR2ThreeState()
        current = model.current_pa((0.0, 1.0, 0.0), -70.0)
        g0_ns = CHR2_3STATE.value_of("g0") / 1000.0
        expected = g0_ns * model.rectification(-70.0) * (-70.0 - CHR2_3STATE.value_of("E"))
        self.assertAlmostEqual(current, expected, places=9)
        self.assertLess(abs(current), 1e5, "a current this large means pS was read as nS")

    def test_irradiance_converts_to_the_expected_photon_flux(self) -> None:
        # 1 mW/mm^2 at 470 nm is ~2.37e15 photons/mm^2/s. Hand-checkable from
        # E = hc/lambda = 4.228e-19 J.
        flux = irradiance_to_photon_flux(1.0)
        self.assertAlmostEqual(flux / 1e15, 2.366, places=2)
        energy = 6.62607015e-34 * 299792458.0 / 470e-9
        self.assertAlmostEqual(flux, 1e-3 / energy, places=0)

    def test_negative_light_is_refused(self) -> None:
        from optosim.units import UnitError

        with self.assertRaises(UnitError):
            irradiance_to_photon_flux(-1.0)

    def test_zero_irradiance_is_zero_flux_not_an_error(self) -> None:
        self.assertEqual(irradiance_to_photon_flux(0.0), 0.0)

    def test_flux_and_irradiance_round_trip(self) -> None:
        from optosim.units import photon_flux_to_irradiance

        self.assertAlmostEqual(
            photon_flux_to_irradiance(irradiance_to_photon_flux(2.5)), 2.5, places=9
        )

    def test_current_density_needs_a_usable_area(self) -> None:
        from optosim.units import opsin_current_density

        self.assertIsNone(opsin_current_density(100.0, 0.0))
        self.assertIsNone(opsin_current_density(100.0, -1.0))
        self.assertIsNone(opsin_current_density(100.0, math.inf))

    def test_current_density_converts_pA_to_uA_per_cm2(self) -> None:
        from optosim.units import opsin_current_density

        # 1000 pA over 1e-4 cm^2 is 1e-3 uA / 1e-4 cm^2 = 10 uA/cm^2.
        self.assertAlmostEqual(opsin_current_density(1000.0, 1e-4), 10.0, places=9)
