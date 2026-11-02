"""Agreement with PyRhO, against the stored reference traces.

The reference traces in ``runs/pyrho_reference_*.json`` were produced by PyRhO
0.9.4 itself, in a separate Python 3.9 environment, integrating its own equations
with its own parameters and its own ``odeint`` call. They are committed so this
test can run anywhere without PyRhO installed -- PyRhO does not import on 3.12.

What this covers that the steady-state check cannot. The steady state agrees to
1e-12, but it is one equation in three rates: an error that shifted ``Ga`` and
``Gr`` in compensating directions would satisfy it and show up only in the
transient. ``Ga`` and ``Gr`` have different Hill exponents (p = 0.8, q = 0.25),
so a cancellation at one irradiance would not survive a tenfold change in flux --
which is why both 1 and 10 mW/mm^2 are checked.
"""

from __future__ import annotations

import json
import pathlib
import unittest

from optosim.opsin.chr2 import ChR2ThreeState, analytic_step_response

REFERENCE_DIR = pathlib.Path(__file__).resolve().parent.parent / "runs"

#: The brief's bar for "correct transcription". The measured deviations are five
#: orders below it and sit at the scale of odeint's own default tolerance.
AGREEMENT_TOLERANCE = 1e-3


def load(irradiance: int) -> dict:
    path = REFERENCE_DIR / f"pyrho_reference_{irradiance}mw.json"
    return json.loads(path.read_text())


class TestReferenceTracesArePresentAndGenuine(unittest.TestCase):
    def test_both_reference_traces_are_committed(self) -> None:
        # If these go missing the agreement claim in the README is unsupported,
        # so their absence must fail rather than silently skip.
        for irradiance in (1, 10):
            path = REFERENCE_DIR / f"pyrho_reference_{irradiance}mw.json"
            self.assertTrue(path.exists(), f"missing reference trace: {path}")

    def test_the_references_say_they_came_from_pyrho(self) -> None:
        for irradiance in (1, 10):
            reference = load(irradiance)
            self.assertEqual(reference["source"], "PyRhO")
            self.assertEqual(reference["pyrho_version"], "0.9.4")
            self.assertEqual(reference["model"], "RhO_3states")
            self.assertEqual(reference["rho_type"], "ChR2")

    def test_the_reference_used_the_same_parameters_this_package_cites(self) -> None:
        # The comparison is only meaningful if both sides ran the same constants.
        # A drift here would make a passing deviation meaningless.
        from optosim.parameters import CHR2_3STATE

        reference = load(1)
        for name in ("g0", "phi_m", "k_a", "k_r", "p", "q", "Gd", "Gr0", "E", "v0", "v1"):
            self.assertAlmostEqual(
                reference["parameters"][name],
                CHR2_3STATE.value_of(name),
                places=10,
                msg=f"{name} differs between the reference and optosim's cited value",
            )

    def test_the_reference_state_order_matches_this_packages(self) -> None:
        # C, O, D. A different ordering would silently compare the wrong column.
        reference = load(1)
        self.assertEqual(reference["state_order"], list(ChR2ThreeState().state_names))
        self.assertEqual(reference["open_state_index"], 1)


class TestTransientAgreement(unittest.TestCase):
    """The trace-against-trace result, at two irradiances."""

    def _deviation(self, irradiance: int) -> tuple[float, float]:
        reference = load(irradiance)
        mine = [
            state[1]
            for state in analytic_step_response(
                ChR2ThreeState(), reference["photon_flux"], reference["times_ms"]
            )
        ]
        worst = 0.0
        worst_at = 0.0
        for time_ms, got, want in zip(
            reference["times_ms"], mine, reference["open_fraction"], strict=True
        ):
            if abs(got - want) > worst:
                worst, worst_at = abs(got - want), time_ms
        return worst, worst_at

    def test_the_transient_agrees_at_one_milliwatt(self) -> None:
        worst, _ = self._deviation(1)
        self.assertLess(worst, AGREEMENT_TOLERANCE, f"max |dO| = {worst:.3e}")

    def test_the_transient_agrees_at_ten_milliwatts(self) -> None:
        # The tenfold change is the point: Ga and Gr scale differently with flux,
        # so compensating rate errors cannot cancel at both.
        worst, _ = self._deviation(10)
        self.assertLess(worst, AGREEMENT_TOLERANCE, f"max |dO| = {worst:.3e}")

    def test_the_measured_deviations_are_at_the_integrator_tolerance_scale(self) -> None:
        # Pinned so a regression that degrades agreement to merely "within 1e-3"
        # is caught. odeint's default tolerance is ~1.5e-8; agreement at that
        # scale means the equations match and the residual is the reference's own
        # integration error.
        for irradiance in (1, 10):
            worst, _ = self._deviation(irradiance)
            self.assertLess(worst, 1e-6, f"at {irradiance} mW/mm^2: {worst:.3e}")

    def test_the_peak_agrees_not_just_the_plateau(self) -> None:
        # The plateau was already known to agree via the closed form. The peak is
        # the part the steady state cannot see.
        for irradiance in (1, 10):
            reference = load(irradiance)
            mine = max(
                state[1]
                for state in analytic_step_response(
                    ChR2ThreeState(), reference["photon_flux"], reference["times_ms"]
                )
            )
            self.assertAlmostEqual(mine, max(reference["open_fraction"]), places=6)

    def test_agreement_is_not_an_artefact_of_comparing_a_trace_with_itself(self) -> None:
        # A deliberately wrong rate must make the comparison fail. Without this,
        # a bug that loaded optosim's own output as the "reference" would pass
        # every test above.
        from optosim.parameters import CHR2_3STATE

        reference = load(1)
        wrong = ChR2ThreeState(parameters=CHR2_3STATE.with_overrides(Gd=0.2))
        mine = [
            state[1]
            for state in analytic_step_response(
                wrong, reference["photon_flux"], reference["times_ms"]
            )
        ]
        worst = max(
            abs(a - b) for a, b in zip(mine, reference["open_fraction"], strict=True)
        )
        self.assertGreater(
            worst,
            AGREEMENT_TOLERANCE,
            "a wrong Gd still agreed with PyRhO, so the comparison proves nothing",
        )


#: Relative to the peak current. The open-fraction deviations sit at ~5.5e-8,
#: and the current is a smooth function of the open fraction at fixed voltage,
#: so the current deviation is expected at the same scale.
CURRENT_TOLERANCE_RELATIVE = 1e-5


def reference_currents_na(reference: dict) -> list[float]:
    """The reference's photocurrent series, in NANOAMPS.

    PyRhO's ``calcI`` returns nA. References written before that was noticed
    stored the series under ``current_pa``, naming the wrong unit while holding
    the right numbers, so both keys are read as nanoamps.
    """
    series = reference.get("current_na")
    if series is None:
        series = reference.get("current_pa")
    if series is None:
        raise AssertionError(
            "the reference carries no current series, so g0, E, v0 and v1 are "
            "untested against PyRhO"
        )
    return series


class TestPhotocurrentAgreement(unittest.TestCase):
    """The current, which is the only thing that exercises g0, E, v0 and v1.

    Comparing the open fraction alone leaves four of the eleven cited ChR2
    parameters transcribed but unchecked, under a heading that claims agreement
    with PyRhO. This closes that.

    The unit is the trap here and is handled explicitly: PyRhO reports nA,
    optosim reports pA, and the two differ by exactly a thousand. A comparison
    that got this wrong would read as a catastrophic disagreement rather than
    as a unit error, which is precisely how the mislabelled key survived.
    """

    def _deviation(self, irradiance: int) -> tuple[float, float, float]:
        reference = load(irradiance)
        model = ChR2ThreeState()
        voltage = reference["voltage_mv"]
        states = analytic_step_response(
            model, reference["photon_flux"], reference["times_ms"]
        )
        worst = 0.0
        worst_at = 0.0
        peak = 0.0
        for time_ms, state, their_na in zip(
            reference["times_ms"], states, reference_currents_na(reference), strict=True
        ):
            mine_pa = model.current_pa((state[0], state[1], state[2]), voltage)
            theirs_pa = their_na * 1000.0
            if abs(mine_pa - theirs_pa) > worst:
                worst, worst_at = abs(mine_pa - theirs_pa), time_ms
            peak = max(peak, abs(theirs_pa))
        return worst, worst_at, worst / peak

    def test_the_photocurrent_agrees_at_one_milliwatt(self) -> None:
        worst, at_ms, relative = self._deviation(1)
        self.assertLess(
            relative,
            CURRENT_TOLERANCE_RELATIVE,
            f"max |dI| = {worst:.3e} pA at t = {at_ms:g} ms ({relative:.3e} of peak)",
        )

    def test_the_photocurrent_agrees_at_ten_milliwatts(self) -> None:
        worst, at_ms, relative = self._deviation(10)
        self.assertLess(
            relative,
            CURRENT_TOLERANCE_RELATIVE,
            f"max |dI| = {worst:.3e} pA at t = {at_ms:g} ms ({relative:.3e} of peak)",
        )

    def test_the_reference_current_is_nanoamps_not_picoamps(self) -> None:
        # The regression that prompted this. PyRhO's peak at 1 mW/mm^2 is about
        # -3.26, and optosim's is about -3263 pA. If a future reference stored
        # picoamps under either key, this fails rather than letting the factor
        # of 1000 be absorbed as a disagreement.
        reference = load(1)
        peak_na = max(abs(value) for value in reference_currents_na(reference))
        self.assertAlmostEqual(peak_na, 3.2635, places=3, msg=f"peak = {peak_na}")

    def test_a_wrong_conductance_is_caught_by_the_current_and_not_by_the_open_fraction(
        self,
    ) -> None:
        # The point of the whole test class: g0 does not appear in the
        # photocycle, so doubling it leaves every open-fraction test passing.
        # Only the current sees it.
        from optosim.parameters import CHR2_3STATE

        reference = load(1)
        wrong = ChR2ThreeState(parameters=CHR2_3STATE.with_overrides(g0=314000.0))
        states = analytic_step_response(
            wrong, reference["photon_flux"], reference["times_ms"]
        )
        voltage = reference["voltage_mv"]
        worst = 0.0
        peak = 0.0
        for state, their_na in zip(
            states, reference_currents_na(reference), strict=True
        ):
            mine_pa = wrong.current_pa((state[0], state[1], state[2]), voltage)
            theirs_pa = their_na * 1000.0
            worst = max(worst, abs(mine_pa - theirs_pa))
            peak = max(peak, abs(theirs_pa))

        open_fraction_unchanged = max(
            abs(state[1] - want)
            for state, want in zip(states, reference["open_fraction"], strict=True)
        )
        self.assertLess(
            open_fraction_unchanged,
            1e-6,
            "doubling g0 changed the open fraction, so this test is not isolating "
            "what it claims to isolate",
        )
        self.assertGreater(
            worst / peak,
            CURRENT_TOLERANCE_RELATIVE,
            "doubling g0 still agreed with PyRhO's current, so the comparison "
            "does not actually test g0",
        )
