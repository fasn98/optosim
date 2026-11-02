"""The numbers in the figure captions, pinned to the runs that produce them.

This exists because of a specific failure. Figure 2's legend said the leaky
integrate-and-fire baseline fires 35 times while the figure drew 39, and the
paper's Results repeated the 35 without saying over what window. Both numbers
were right -- 35 during the light pulse, 39 over the 120 ms simulated -- but
nothing checked the caption against the run, so the two could disagree
indefinitely and did.

The anti-fabrication guard cannot catch this. It verifies that no module on an
output path can reach a random-number generator; it has nothing to say about a
sentence describing a figure. A statement in prose is only ever checked by
something that recomputes what it claims, which is what these tests do.

If one of these fails, the figure and the caption have come apart. Regenerate
with ``scripts/make_figures.py``, which prints every caption number, and fix
whichever of the two is wrong.
"""

from __future__ import annotations

import unittest

from optosim.light import IrradianceProtocol
from optosim.simulate import simulate_hh, simulate_lif

# Imported from the figure script rather than restated, so there is one
# definition of the protocol and a test cannot drift from the picture.
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

from make_figures import (  # noqa: E402
    COUPLED_DT_MS,
    COUPLED_IRRADIANCE_MW_MM2,
    COUPLED_PULSE_DURATION_MS,
    COUPLED_PULSE_ONSET_MS,
    COUPLED_RUN_DURATION_MS,
)


class TestFigureTwoProtocol(unittest.TestCase):
    """The coupled-response figure: what was run, and what it counts."""

    @classmethod
    def setUpClass(cls) -> None:
        protocol = IrradianceProtocol.single(
            COUPLED_PULSE_ONSET_MS,
            COUPLED_PULSE_DURATION_MS,
            COUPLED_IRRADIANCE_MW_MM2,
        )
        cls.offset_ms = COUPLED_PULSE_ONSET_MS + COUPLED_PULSE_DURATION_MS
        cls.hh = simulate_hh(
            protocol, duration_ms=COUPLED_RUN_DURATION_MS, dt_ms=COUPLED_DT_MS
        )
        cls.lif = simulate_lif(
            protocol, duration_ms=COUPLED_RUN_DURATION_MS, dt_ms=COUPLED_DT_MS
        )

    def test_the_run_extends_past_light_off(self) -> None:
        # If it did not, the during-pulse and whole-window counts would coincide
        # and the ambiguity this file guards against could not arise -- but the
        # figure's x axis would also stop at 100 ms, which is not the figure.
        self.assertGreater(COUPLED_RUN_DURATION_MS, self.offset_ms)

    def test_the_hodgkin_huxley_neuron_fires_once_and_blocks(self) -> None:
        self.assertEqual(len(self.hh.spike_times_ms), 1)

    def test_the_baseline_fires_35_times_during_the_pulse(self) -> None:
        during = [t for t in self.lif.spike_times_ms if t < self.offset_ms]
        self.assertEqual(len(during), 35, "the caption's during-pulse count")

    def test_the_baseline_fires_39_times_over_the_whole_window(self) -> None:
        # The number the figure actually draws. A caption naming only 35 is
        # describing a different window from the one plotted.
        self.assertEqual(len(self.lif.spike_times_ms), 39, "what the figure shows")

    def test_the_four_extra_spikes_all_follow_light_off(self) -> None:
        after = [t for t in self.lif.spike_times_ms if t >= self.offset_ms]
        self.assertEqual(len(after), 4)
        self.assertTrue(all(t > self.offset_ms for t in after))

    def test_a_longer_run_adds_no_further_spikes(self) -> None:
        # The reviewer reported 39 for a 200 ms run. The baseline stops firing
        # about 14 ms after light-off, so any window past ~114 ms gives 39 --
        # which is why the 120 ms window is not an arbitrary truncation.
        protocol = IrradianceProtocol.single(
            COUPLED_PULSE_ONSET_MS,
            COUPLED_PULSE_DURATION_MS,
            COUPLED_IRRADIANCE_MW_MM2,
        )
        longer = simulate_lif(protocol, duration_ms=200.0, dt_ms=COUPLED_DT_MS)
        self.assertEqual(len(longer.spike_times_ms), len(self.lif.spike_times_ms))
        self.assertLess(max(longer.spike_times_ms), COUPLED_RUN_DURATION_MS)


if __name__ == "__main__":
    unittest.main()
