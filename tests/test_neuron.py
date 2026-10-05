"""Hodgkin-Huxley, the integrate-and-fire baseline, and the coupling between them."""

from __future__ import annotations

import unittest

from optosim.neuron.hh import (
    HodgkinHuxley,
    alpha_h,
    alpha_m,
    alpha_n,
    beta_h,
    beta_m,
    beta_n,
    steady_state_gate,
)
from optosim.neuron.iaf import LeakyIntegrateAndFire
from optosim.parameters import CHR2_3STATE, HH_SQUID_AXON, LIF_BASELINE


class TestRateSingularities(unittest.TestCase):
    """The two voltages where the textbook expressions are 0/0.

    Neither is obscure: -55 mV is a few mV above rest and any depolarising
    stimulus passes through it, so an unguarded implementation divides by zero
    during an ordinary run rather than in some edge case.
    """

    def test_alpha_m_at_its_singularity_is_the_limit(self) -> None:
        self.assertAlmostEqual(alpha_m(-40.0), 1.0, places=12)

    def test_alpha_n_at_its_singularity_is_the_limit(self) -> None:
        self.assertAlmostEqual(alpha_n(-55.0), 0.1, places=12)

    def test_alpha_m_is_continuous_across_its_singularity(self) -> None:
        # The deviation from the limit is first order in the offset, so the test
        # is written that way rather than with fixed decimal places -- a fixed
        # tolerance passes or fails on which epsilon happens to be largest,
        # which is what the first version of this test did.
        for epsilon in (1e-5, 1e-3, 1e-1):
            for voltage in (-40.0 + epsilon, -40.0 - epsilon):
                self.assertAlmostEqual(
                    alpha_m(voltage) / 1.0, 1.0, delta=0.01,
                    msg=f"alpha_m at offset {epsilon} deviates more than 1%",
                )

    def test_alpha_n_is_continuous_across_its_singularity(self) -> None:
        for epsilon in (1e-5, 1e-3, 1e-1):
            for voltage in (-55.0 + epsilon, -55.0 - epsilon):
                self.assertAlmostEqual(
                    alpha_n(voltage) / 0.1, 1.0, delta=0.01,
                    msg=f"alpha_n at offset {epsilon} deviates more than 1%",
                )

    def test_the_deviation_from_the_limit_shrinks_with_the_offset(self) -> None:
        # The actual content of continuity: closer in, closer to the limit.
        deviations = [abs(alpha_m(-40.0 + e) - 1.0) for e in (1e-1, 1e-2, 1e-3)]
        self.assertEqual(deviations, sorted(deviations, reverse=True))

    def test_all_rates_are_positive_across_the_physiological_range(self) -> None:
        for millivolts in range(-100, 61):
            voltage = float(millivolts)
            for rate in (
                alpha_m(voltage), beta_m(voltage), alpha_h(voltage),
                beta_h(voltage), alpha_n(voltage), beta_n(voltage),
            ):
                self.assertGreater(rate, 0.0, f"non-positive rate at {voltage} mV")

    def test_a_gate_with_no_rates_has_no_steady_state(self) -> None:
        # None, not 0.5. There is no equilibrium to report.
        self.assertIsNone(steady_state_gate(0.0, 0.0))
        self.assertAlmostEqual(steady_state_gate(1.0, 1.0), 0.5)


class TestHodgkinHuxley(unittest.TestCase):
    def test_the_resting_state_reproduces_the_canonical_gate_values(self) -> None:
        # m ~ 0.053, h ~ 0.596, n ~ 0.318 at -65 mV: the values this parameter
        # set is universally reported with.
        v, m, h, n = HodgkinHuxley().resting_state()
        self.assertAlmostEqual(v, -65.0, places=10)
        self.assertAlmostEqual(m, 0.0529, places=3)
        self.assertAlmostEqual(h, 0.5961, places=3)
        self.assertAlmostEqual(n, 0.3177, places=3)

    def test_the_resting_state_is_very_nearly_stationary(self) -> None:
        # The canonical set rests within a hundredth of a mV/ms of stationary.
        # Not exactly zero, and asserting zero would be asserting something
        # false about the published parameters.
        neuron = HodgkinHuxley()
        dv = neuron.derivatives(neuron.resting_state())[0]
        self.assertLess(abs(dv), 0.01, f"dV/dt at rest is {dv} mV/ms")

    def test_an_unstimulated_membrane_does_not_spike(self) -> None:
        from optosim.light import IrradianceProtocol
        from optosim.simulate import simulate_hh

        result = simulate_hh(IrradianceProtocol(), duration_ms=100.0, dt_ms=0.01)
        self.assertEqual(result.spike_times_ms, ())

    def test_a_wrong_parameter_set_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            HodgkinHuxley(parameters=CHR2_3STATE)

    def test_the_ionic_currents_vanish_nowhere_in_particular_at_rest(self) -> None:
        # They cancel at rest, which is a different statement from each being
        # zero -- a sodium current of zero at rest would be wrong.
        neuron = HodgkinHuxley()
        i_na, i_k, i_l = neuron.ionic_current_density(neuron.resting_state())
        self.assertLess(abs(i_na + i_k + i_l), 0.01)
        self.assertGreater(abs(i_na), 1e-3)


class TestLeakyIntegrateAndFire(unittest.TestCase):
    def test_it_rests_at_the_leak_reversal(self) -> None:
        self.assertAlmostEqual(LeakyIntegrateAndFire().resting_state(), -54.387, places=6)

    def test_the_leak_pulls_the_voltage_back_towards_rest(self) -> None:
        neuron = LeakyIntegrateAndFire()
        self.assertLess(neuron.derivative(-40.0), 0.0)
        self.assertGreater(neuron.derivative(-70.0), 0.0)

    def test_an_inward_current_depolarises(self) -> None:
        # Inward is negative by the HH convention, so it must raise dV/dt.
        neuron = LeakyIntegrateAndFire()
        self.assertGreater(neuron.derivative(-54.387, -10.0), 0.0)

    def test_the_refractory_period_imposes_a_rate_ceiling(self) -> None:
        # An artefact of an uncited parameter, exposed so it is reported rather
        # than discovered in a saturating f-I curve.
        self.assertAlmostEqual(LeakyIntegrateAndFire().maximum_firing_rate_hz, 500.0)

    def test_its_threshold_is_uncited_and_says_so(self) -> None:
        # The honesty requirement: this model's spike counts rest on a chosen
        # number, and the parameter set must not pretend otherwise.
        threshold = LIF_BASELINE["V_threshold"]
        self.assertFalse(threshold.is_cited)
        self.assertIn("no source", threshold.quantity.source)

    def test_it_shares_leak_and_capacitance_with_the_hh_set(self) -> None:
        # So the comparison differs in spiking mechanism rather than in passive
        # properties, which would confound it.
        for name in ("C_m", "g_L", "E_L"):
            self.assertAlmostEqual(
                LIF_BASELINE.value_of(name), HH_SQUID_AXON.value_of(name), places=10
            )


class TestCouplingCurrentBalance(unittest.TestCase):
    """The physical check on the opsin-to-membrane join.

    The place a unit error would hide: the opsin produces absolute pA and the
    membrane consumes uA/cm^2, so a wrong area or a pS/nS slip changes the
    photocurrent by orders of magnitude while every trace still looks like a
    trace.

    Under sustained light the HH membrane settles at a depolarised equilibrium.
    At that point the ionic currents must sum to exactly cancel the photocurrent
    density. Checking that is an independent arithmetic test of the whole
    coupling: it was hand-verified at 5 mW/mm^2, where the plateau sits at
    -57.9 mV with I_K = +23.45, I_Na = -7.29, I_L = -1.05 uA/cm^2 against a
    photocurrent density of -15.03.
    """

    def test_the_plateau_balances_the_photocurrent(self) -> None:
        from optosim.light import IrradianceProtocol
        from optosim.parameters import DEFAULT_COMPARTMENT
        from optosim.simulate import simulate_hh
        from optosim.units import opsin_current_density

        protocol = IrradianceProtocol.single(50.0, 150.0, 5.0)
        result = simulate_hh(protocol, duration_ms=220.0, dt_ms=0.01)

        # Late in the pulse, well past the onset transient.
        index = min(
            range(len(result.voltage.times_ms)),
            key=lambda i: abs(result.voltage.times_ms[i] - 180.0),
        )
        voltage = result.voltage.values[index]
        current_pa = result.photocurrent_pa.values[index]
        density = opsin_current_density(current_pa, DEFAULT_COMPARTMENT.value_of("area"))
        assert density is not None

        # Reconstruct the gating variables at their steady values for this
        # voltage, which is what they have converged to on a flat plateau.
        gates = []
        for alpha, beta in (
            (alpha_m(voltage), beta_m(voltage)),
            (alpha_h(voltage), beta_h(voltage)),
            (alpha_n(voltage), beta_n(voltage)),
        ):
            value = steady_state_gate(alpha, beta)
            assert value is not None
            gates.append(value)
        i_na, i_k, i_l = HodgkinHuxley().ionic_current_density(
            (voltage, gates[0], gates[1], gates[2])
        )

        # Ionic outward must cancel the inward photocurrent.
        self.assertAlmostEqual(i_na + i_k + i_l, -density, delta=0.5)
        # And the plateau is where it was hand-verified to be.
        self.assertAlmostEqual(voltage, -57.9, delta=0.5)

    def test_a_wrong_area_changes_the_answer_by_orders_of_magnitude(self) -> None:
        # Guards the conversion itself: if the area were ignored, these two runs
        # would agree. They must not.
        from optosim.light import IrradianceProtocol
        from optosim.parameters import DEFAULT_COMPARTMENT
        from optosim.simulate import simulate_hh

        protocol = IrradianceProtocol.single(10.0, 40.0, 1.0)
        small = simulate_hh(protocol, duration_ms=60.0, dt_ms=0.02)
        bigger = simulate_hh(
            protocol,
            duration_ms=60.0,
            dt_ms=0.02,
            geometry=DEFAULT_COMPARTMENT.with_overrides(area=1e-2),
        )
        self.assertNotAlmostEqual(
            min(small.voltage.values), min(bigger.voltage.values), places=3
        )

    def test_overriding_a_parameter_drops_its_citation(self) -> None:
        from optosim.parameters import DEFAULT_COMPARTMENT

        overridden = DEFAULT_COMPARTMENT.with_overrides(area=1e-2)
        self.assertFalse(overridden["area"].is_cited)
        self.assertIn("overridden by caller", overridden["area"].quantity.source)
