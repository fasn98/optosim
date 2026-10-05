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
