#!/usr/bin/env python3
"""Regenerate both figures in the paper, from the committed code and data.

Why this exists. The figures were produced once, by hand, and the scripts that
made them were never committed -- so the two pictures in the paper could not be
reproduced from this repository, and the protocol behind Figure 2 lived only in
the figure itself. A reviewer reading the paper could not tell what window the
spike count referred to, and neither could the author.

That gap has a cost that was paid: Figure 2's legend says the leaky
integrate-and-fire baseline fires 35 times while the figure draws 39. Both
numbers are right about different windows -- 35 during the light pulse, 39 over
the whole 120 ms simulated -- and nothing in the repository said which one the
caption meant. This script states the protocol as data, prints every number that
appears in a caption, and draws the figures from the same run, so the two cannot
drift apart again.

Usage:
    python scripts/make_figures.py                   # writes figures/
    python scripts/make_figures.py --output-dir .    # next to content.tex
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from optosim.light import IrradianceProtocol  # noqa: E402
from optosim.opsin.chr2 import ChR2ThreeState, analytic_step_response  # noqa: E402
from optosim.simulate import simulate_hh, simulate_lif  # noqa: E402
from optosim.units import irradiance_to_photon_flux  # noqa: E402

# ---------------------------------------------------------------------------
# The protocols, stated once, here, as the single source for both the figures
# and the sentences in the paper that describe them.
# ---------------------------------------------------------------------------

#: Figure 1: a step of constant irradiance held for the whole run.
PHOTOCYCLE_IRRADIANCES_MW_MM2 = (1.0, 10.0)
PHOTOCYCLE_DURATION_MS = 200.0
PHOTOCYCLE_DT_MS = 0.01

#: Figure 2. The pulse starts at t = 0 and the run continues past light-off,
#: which is the whole reason the spike count is ambiguous without this block:
#: the baseline keeps firing for about 14 ms after the light goes out.
COUPLED_PULSE_ONSET_MS = 0.0
COUPLED_PULSE_DURATION_MS = 100.0
COUPLED_IRRADIANCE_MW_MM2 = 5.0
COUPLED_RUN_DURATION_MS = 120.0
COUPLED_DT_MS = 0.01

#: Every 250th sample of the reference is drawn, so the markers are legible at
#: 2.5 ms spacing rather than 20001 overlapping circles.
REFERENCE_MARKER_STRIDE = 250

OPTOSIM_COLOURS = {1.0: "#2e7d8f", 10.0: "#d2601a"}
HH_COLOUR = "#2e7d8f"
LIF_COLOUR = "#9a9a9a"
LIGHT_COLOUR = "#f3e6c8"


def load_reference(irradiance: float, directory: pathlib.Path) -> dict | None:
    path = directory / f"pyrho_reference_{irradiance:g}mw.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    if payload.get("source") != "PyRhO":
        raise ValueError(f"{path} is not a PyRhO reference")
    return payload


def figure_photocycle(output: pathlib.Path, reference_dir: pathlib.Path) -> None:
    """Figure 1: optosim against the PyRhO reference, at two irradiances."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    model = ChR2ThreeState()
    figure, axes = plt.subplots(figsize=(7.2, 3.6))

    print("\nFigure 1 -- three-state photocycle")
    print(f"  step of constant irradiance, {PHOTOCYCLE_DURATION_MS:g} ms, "
          f"dt = {PHOTOCYCLE_DT_MS:g} ms")

    reference_plotted = False
    for irradiance in PHOTOCYCLE_IRRADIANCES_MW_MM2:
        flux = irradiance_to_photon_flux(irradiance)
        steps = int(round(PHOTOCYCLE_DURATION_MS / PHOTOCYCLE_DT_MS))
        times = [index * PHOTOCYCLE_DT_MS for index in range(steps + 1)]
        opens = [state[1] for state in analytic_step_response(model, flux, times)]
        axes.plot(
            times,
            opens,
            color=OPTOSIM_COLOURS[irradiance],
            linewidth=1.6,
            label=f"optosim, {irradiance:g} mW/mm$^2$",
            zorder=3,
        )

        reference = load_reference(irradiance, reference_dir)
        if reference is None:
            print(f"  {irradiance:g} mW/mm^2: NO REFERENCE -- the figure would show "
                  "optosim alone, which is not what its caption claims")
            continue
        marker_times = reference["times_ms"][::REFERENCE_MARKER_STRIDE]
        marker_opens = reference["open_fraction"][::REFERENCE_MARKER_STRIDE]
        axes.plot(
            marker_times,
            marker_opens,
            linestyle="none",
            marker="o",
            markersize=3.4,
            markerfacecolor="none",
            markeredgecolor="#333333",
            markeredgewidth=0.7,
            label="PyRhO reference" if not reference_plotted else None,
            zorder=4,
        )
        reference_plotted = True

        worst = max(
            abs(mine - theirs)
            for mine, theirs in zip(
                [
                    state[1]
                    for state in analytic_step_response(
                        model, reference["photon_flux"], reference["times_ms"]
                    )
                ],
                reference["open_fraction"],
                strict=True,
            )
        )
        peak = max(reference["open_fraction"])
        plateau = reference["open_fraction"][-1]
        print(f"  {irradiance:g} mW/mm^2: peak O = {peak:.4f}, plateau O = {plateau:.4f}, "
              f"peak/plateau = {peak / plateau:.2f}, max |dO| = {worst:.3e}")

    axes.set_xlabel("time (ms)")
    axes.set_ylabel("open fraction $O$")
    axes.set_title("Three-state ChR2 photocycle: peak-then-plateau")
    axes.set_xlim(0, PHOTOCYCLE_DURATION_MS)
    axes.set_ylim(bottom=0.0)
    axes.grid(True, linewidth=0.4, alpha=0.4)
    axes.legend(frameon=True, fontsize=9)
    figure.tight_layout()
    figure.savefig(output, dpi=200)
    plt.close(figure)
    print(f"  written to {output}")


def figure_coupled(output: pathlib.Path) -> dict[str, float]:
    """Figure 2: the coupled models under one pulse. Returns the caption numbers."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    protocol = IrradianceProtocol.single(
        COUPLED_PULSE_ONSET_MS, COUPLED_PULSE_DURATION_MS, COUPLED_IRRADIANCE_MW_MM2
    )
    offset_ms = COUPLED_PULSE_ONSET_MS + COUPLED_PULSE_DURATION_MS

    hh = simulate_hh(protocol, duration_ms=COUPLED_RUN_DURATION_MS, dt_ms=COUPLED_DT_MS)
    lif = simulate_lif(protocol, duration_ms=COUPLED_RUN_DURATION_MS, dt_ms=COUPLED_DT_MS)

    def split(spikes: tuple[float, ...]) -> tuple[int, int]:
        during = sum(1 for t in spikes if COUPLED_PULSE_ONSET_MS <= t < offset_ms)
        return during, len(spikes) - during

    hh_during, hh_after = split(hh.spike_times_ms)
    lif_during, lif_after = split(lif.spike_times_ms)

    print("\nFigure 2 -- coupled response")
    print(f"  pulse: {COUPLED_IRRADIANCE_MW_MM2:g} mW/mm^2 from "
          f"{COUPLED_PULSE_ONSET_MS:g} to {offset_ms:g} ms")
    print(f"  run:   0 to {COUPLED_RUN_DURATION_MS:g} ms, dt = {COUPLED_DT_MS:g} ms")
    print(f"  Hodgkin-Huxley: {len(hh.spike_times_ms)} spike(s) in the window "
          f"({hh_during} during the pulse, {hh_after} after light-off)")
    print(f"  LIF baseline:   {len(lif.spike_times_ms)} spike(s) in the window "
          f"({lif_during} during the pulse, {lif_after} after light-off)")
    if lif.spike_times_ms:
        print(f"  last LIF spike at {max(lif.spike_times_ms):.2f} ms")
    print("  NOTE: a caption that says only '35 spikes' is naming the during-pulse\n"
          "  count while the figure draws the whole window. Say both.")

    figure, axes = plt.subplots(figsize=(7.2, 3.6))
    axes.axvspan(
        COUPLED_PULSE_ONSET_MS,
        offset_ms,
        color=LIGHT_COLOUR,
        zorder=0,
        label=f"light on ({COUPLED_IRRADIANCE_MW_MM2:g} mW/mm$^2$)",
    )
    axes.plot(
        lif.voltage.times_ms,
        lif.voltage.values,
        color=LIF_COLOUR,
        linewidth=0.9,
        label=f"LIF baseline ({lif_during} spikes during, {len(lif.spike_times_ms)} total)",
        zorder=2,
    )
    axes.plot(
        hh.voltage.times_ms,
        hh.voltage.values,
        color=HH_COLOUR,
        linewidth=1.1,
        label=f"Hodgkin-Huxley ({len(hh.spike_times_ms)} spike)",
        zorder=3,
    )
    axes.set_xlabel("time (ms)")
    axes.set_ylabel("membrane potential (mV)")
    axes.set_title(
        f"Coupled response to a {COUPLED_PULSE_DURATION_MS:g} ms light pulse"
    )
    axes.set_xlim(0, COUPLED_RUN_DURATION_MS)
    axes.grid(True, linewidth=0.4, alpha=0.4)
    axes.legend(frameon=True, fontsize=8, loc="upper right")
    figure.tight_layout()
    figure.savefig(output, dpi=200)
    plt.close(figure)
    print(f"  written to {output}")

    return {
        "hh_total": len(hh.spike_times_ms),
        "lif_during_pulse": lif_during,
        "lif_total": len(lif.spike_times_ms),
        "lif_last_spike_ms": max(lif.spike_times_ms) if lif.spike_times_ms else 0.0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="figures")
    parser.add_argument("--reference-dir", default="runs")
    args = parser.parse_args()

    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figure_photocycle(
        output_dir / "fig1_photocycle.png", pathlib.Path(args.reference_dir)
    )
    figure_coupled(output_dir / "fig2_coupled.png")

    print(
        "\nEvery number printed above appears in a caption or in the Results. "
        "They come from\nthis run, so a caption that disagrees with them is the "
        "caption's error, not a mystery."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
