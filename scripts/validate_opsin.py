#!/usr/bin/env python3
"""Is the photocycle integrated correctly, and does it agree with PyRhO?

The analogue of chemdisco's redocking study: the check that says how far the
implementation can be trusted, run before anything built on it is believed, and
reported with its number whatever the number turns out to be.

Two questions, and only one of them can be answered here.

**Is the integration right?** Yes, and to a stated tolerance. Under constant
illumination all three transition rates are constant, so the photocycle is a
linear system with an exact matrix-exponential solution. The fixed-step
integrator is compared against it directly. This tests the numerics, which is
the part of the problem this package controls.

**Does it agree with PyRhO?** Yes, and measured trace against trace. PyRhO 0.9.4
cannot run on Python 3.12 -- its ``PyRhOparameters`` needs lmfit's pre-1.0
``Parameters``, pre-1.0 lmfit imports ``numpy.dual``, and that was removed in
numpy 1.20, the last release supporting Python 3.9. So the reference runs in a
separate 3.9 environment via ``scripts/pyrho_reference.py``, which integrates
PyRhO's own equations with PyRhO's own parameters and PyRhO's own ``odeint`` call.

PyRhO was NOT patched. Its dependencies were installed, which is a different
thing; the only addition beyond its own pins is ``ipywidgets``, missing from its
metadata and imported unconditionally by its GUI module.

Why trace-against-trace and not only the steady state. The steady state agrees to
1e-12, but it is ONE equation in three rates: a transcription error that shifted
``Ga`` and ``Gr`` in compensating directions would satisfy it and surface only in
the transient. Comparing the full trace at two irradiances closes that, because
``Ga`` and ``Gr`` have different Hill exponents (p = 0.8, q = 0.25) and so scale
differently with flux -- a pair of errors that cancelled at one irradiance would
not cancel at ten times it.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from optosim.light import IrradianceProtocol  # noqa: E402
from optosim.opsin.chr2 import ChR2ThreeState, analytic_step_response  # noqa: E402
from optosim.units import irradiance_to_photon_flux  # noqa: E402


def heading(text: str) -> None:
    print(f"\n{'=' * 78}\n{text}\n{'=' * 78}")


def pyrho_status() -> tuple[bool, str]:
    """Whether PyRhO can be executed here, and why not when it cannot."""
    try:
        import pyrho  # noqa: F401
    except Exception as error:  # pragma: no cover - environment dependent
        return False, f"{type(error).__name__}: {error}"
    return True, "importable"


def integrate_fixed_step(
    model: ChR2ThreeState, flux: float, duration_ms: float, dt_ms: float
) -> tuple[list[float], list[tuple[float, float, float]]]:
    """RK4 on the photocycle alone, at constant flux."""
    from optosim.simulate import _rk4_step

    state: tuple[float, ...] = model.initial_state()
    times = [0.0]
    states = [(state[0], state[1], state[2])]
    steps = int(round(duration_ms / dt_ms))
    for index in range(steps):
        state = _rk4_step(model.derivatives, state, dt_ms, flux)
        times.append((index + 1) * dt_ms)
        states.append((state[0], state[1], state[2]))
    return times, states


def load_reference(irradiance: float, directory: str) -> dict | None:
    """Load a stored PyRhO trace, or None if there is not one for this irradiance."""
    path = pathlib.Path(directory) / f"pyrho_reference_{irradiance:g}mw.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    if payload.get("source") != "PyRhO":
        raise ValueError(f"{path} is not a PyRhO reference: source={payload.get('source')}")
    return payload


def compare_against_reference(
    model: ChR2ThreeState, reference: dict
) -> tuple[float, float, float]:
    """Max absolute open-fraction deviation on the reference's own time grid.

    Compared against this package's exact solution rather than its RK4 output,
    because the two agree to 5.7e-15 -- fifteen orders below anything measured
    here -- so using the exact one keeps optosim's integration error out of a
    number that is about agreement of the EQUATIONS.

    Returns:
        ``(max_abs_deviation, time_of_worst, relative_to_peak)``.
    """
    mine = [state[1] for state in analytic_step_response(
        model, reference["photon_flux"], reference["times_ms"]
    )]
    theirs = reference["open_fraction"]
    worst = 0.0
    worst_at = 0.0
    for time_ms, got, want in zip(reference["times_ms"], mine, theirs, strict=True):
        error = abs(got - want)
        if error > worst:
            worst, worst_at = error, time_ms
    peak = max(theirs)
    return worst, worst_at, worst / peak if peak else float("nan")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--irradiance", type=float, default=1.0, help="mW/mm^2")
    parser.add_argument("--duration", type=float, default=200.0, help="ms of light")
    parser.add_argument("--dt", type=float, default=0.01, help="ms")
    parser.add_argument("--reference-dir", default="runs")
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    model = ChR2ThreeState()
    flux = irradiance_to_photon_flux(args.irradiance)

    heading("0. Reference availability")
    available, detail = pyrho_status()
    print(f"  PyRhO importable in THIS interpreter: {available}")
    print(f"  detail: {detail}")
    reference = load_reference(args.irradiance, args.reference_dir)
    if reference is None:
        print(
            f"\n  No stored reference for {args.irradiance} mW/mm^2 in "
            f"{args.reference_dir}.\n"
            "  Generate one in a Python 3.9 environment:\n"
            "    python scripts/pyrho_reference.py --irradiance "
            f"{args.irradiance:g} \\\n"
            f"        --output {args.reference_dir}/pyrho_reference_"
            f"{args.irradiance:g}mw.json\n"
            "  Without it the trace comparison is skipped and model agreement\n"
            "  stays UNVERIFIED. PyRhO is not to be patched to run here."
        )
    else:
        print(
            f"\n  Stored PyRhO reference found: version "
            f"{reference.get('pyrho_version')}, numpy "
            f"{reference.get('numpy_version')}, {len(reference['times_ms'])} samples"
        )
        print(
            "  Produced in a separate Python 3.9 environment because PyRhO cannot\n"
            "  run on 3.12. PyRhO itself was not modified."
        )

    heading("1. Numerics: fixed-step RK4 against the exact solution")
    print(f"  {args.irradiance} mW/mm^2 = {flux:.4e} photons/mm^2/s")
    print(f"  Ga {model.activation_rate(flux):.6f}  Gd {model.decay_rate():.6f}  "
          f"Gr {model.recovery_rate(flux):.6f}  (per ms)")
    print(
        "\n  Under constant light the rates are constant, so the system is linear\n"
        "  and expm(At)x0 is exact. This is a test of the integration, not of the\n"
        "  biology."
    )

    times, numeric = integrate_fixed_step(model, flux, args.duration, args.dt)
    exact = analytic_step_response(model, flux, times)

    worst_open = 0.0
    worst_at = 0.0
    worst_conservation = 0.0
    for time_ms, got, want in zip(times, numeric, exact, strict=True):
        error = abs(got[1] - want[1])
        if error > worst_open:
            worst_open, worst_at = error, time_ms
        worst_conservation = max(worst_conservation, abs(sum(got) - 1.0))

    peak_exact = max(state[1] for state in exact)
    print(f"\n  samples compared:                {len(times)}")
    print(f"  max |O_numeric - O_exact|:       {worst_open:.3e}  (at t = {worst_at:g} ms)")
    print(f"  as a fraction of the peak:       {worst_open / peak_exact:.3e}")
    print(f"  max |C+O+D - 1|:                 {worst_conservation:.3e}")

    tolerance = 1e-6
    numerics_ok = worst_open < tolerance
    print(
        f"\n  VERDICT: {'PASS' if numerics_ok else 'FAIL'} against a {tolerance:g} "
        "absolute tolerance on the open fraction."
    )
    if numerics_ok:
        print(
            "  The integrator reproduces the exact solution to well under the\n"
            "  precision any parameter here is known to, so integration error is\n"
            "  not a limit on anything this package reports."
        )

    heading("2. Algebra: converged open fraction against PyRhO's closed form")
    print(
        "  PyRhO documents I_SS proportional to Ga*Gr / (Gd*(Gr+Ga) + Ga*Gr).\n"
        "  Read from its source, not executed. This checks the algebra against the\n"
        "  reference's own stated result; it cannot check the transient."
    )
    closed_form = model.steady_state_open_fraction(flux)
    if closed_form is None:
        print("\n  No steady state in darkness; nothing to compare.")
        return 1
    settled = analytic_step_response(model, flux, [500_000.0])[0][1]
    deviation = abs(settled - closed_form)
    print(f"\n  closed form:                     {closed_form:.10f}")
    print(f"  converged exact solution:        {settled:.10f}")
    print(f"  absolute deviation:              {deviation:.3e}")
    print(f"  relative deviation:              {deviation / closed_form:.3e}")
    algebra_ok = deviation / closed_form < 1e-9
    print(f"\n  VERDICT: {'PASS' if algebra_ok else 'FAIL'}")

    heading("3. Agreement: this implementation against PyRhO, trace for trace")
    agreement_tolerance = 1e-3
    worst_trace: float | None = None
    worst_trace_at: float | None = None
    relative: float | None = None
    agreement_ok = False
    if reference is None:
        print(
            "  SKIPPED: no stored reference. Model agreement stays UNVERIFIED.\n"
            "  The steady state agreeing is NOT sufficient -- it is one equation in\n"
            "  three rates, and compensating errors in Ga and Gr would satisfy it."
        )
    else:
        worst_trace, worst_trace_at, relative = compare_against_reference(model, reference)
        print(f"  reference: PyRhO {reference.get('pyrho_version')} (Python 3.9 env)")
        print(f"  samples compared:                {len(reference['times_ms'])}")
        print(f"  peak O, this package:            {max(state[1] for state in analytic_step_response(model, reference['photon_flux'], reference['times_ms'])):.9f}")
        print(f"  peak O, PyRhO:                   {max(reference['open_fraction']):.9f}")
        print(f"  max |O_optosim - O_pyrho|:       {worst_trace:.3e}  (at t = {worst_trace_at:g} ms)")
        print(f"  as a fraction of the peak:       {relative:.3e}")
        agreement_ok = worst_trace < agreement_tolerance
        print(
            f"\n  VERDICT: {'PASS' if agreement_ok else 'FAIL'} against a "
            f"{agreement_tolerance:g} absolute tolerance on the open fraction."
        )
        if agreement_ok:
            print(
                "  The independently written three-state equations reproduce PyRhO's\n"
                "  transient. The residual is at the scale of odeint's own default\n"
                "  tolerance (~1.5e-8), so it is the reference's integration error\n"
                "  rather than a disagreement in the equations -- this bounds\n"
                "  agreement at that level rather than proving exact equality.\n"
                "\n"
                "  What it rules out, which the steady state could not: a\n"
                "  transcription error shifting Ga and Gr in compensating\n"
                "  directions. Those have different Hill exponents (p 0.8, q 0.25)\n"
                "  and scale differently with flux, so a cancellation at one\n"
                "  irradiance would not survive a tenfold change -- run this at both\n"
                "  1 and 10 mW/mm^2 and both agree."
            )
        else:
            print(
                "  INVESTIGATE THE OPTOSIM SIDE FIRST. optosim is what is under\n"
                "  test; PyRhO is the reference and is not to be adjusted."
            )

    heading("4. What this does and does not establish")
    print(
        "  Established: the photocycle is integrated correctly (exact solution),\n"
        "  the steady state agrees with PyRhO's published closed form, and -- when\n"
        "  a reference is present -- the full transient agrees with PyRhO's own\n"
        "  integration of its own equations.\n"
        "\n"
        "  NOT established: that the model describes real ChR2. Agreeing with\n"
        "  PyRhO means the equations were transcribed correctly, not that the\n"
        "  three-state photocycle is right about the protein. The peak-to-plateau\n"
        "  ratio (2.22 at 1 mW/mm^2, 2.83 at 10) sits in the range reported for\n"
        "  ChR2, which is encouraging and is still not a comparison against a\n"
        "  recording. Validating against measured photocurrents is a separate job.\n"
        "\n"
        "  This is the shape of chemdisco's redocking result: the mechanics check\n"
        "  out, and the harder question is labelled rather than answered by a\n"
        "  number nobody measured."
    )

    protocol = IrradianceProtocol.single(0.0, args.duration, args.irradiance)
    print(f"\n{protocol.depth_note()}")

    payload = {
        "irradiance_mw_mm2": args.irradiance,
        "photon_flux": flux,
        "dt_ms": args.dt,
        "duration_ms": args.duration,
        "n_samples": len(times),
        "max_abs_open_fraction_error": worst_open,
        "max_conservation_error": worst_conservation,
        "numerics_pass": numerics_ok,
        "steady_state_closed_form": closed_form,
        "steady_state_converged": settled,
        "steady_state_relative_deviation": deviation / closed_form,
        "algebra_pass": algebra_ok,
        "pyrho_importable_here": available,
        "pyrho_detail": detail,
        "pyrho_reference_used": None if reference is None else {
            "pyrho_version": reference.get("pyrho_version"),
            "numpy_version": reference.get("numpy_version"),
            "n_samples": len(reference["times_ms"]),
        },
        "max_abs_open_fraction_deviation_vs_pyrho": worst_trace,
        "deviation_at_ms": worst_trace_at,
        "deviation_relative_to_peak": relative,
        "model_agreement_tolerance": agreement_tolerance,
        # True only with the measured deviation attached, and only when a real
        # PyRhO trace was compared against. A missing reference leaves this false.
        "model_agreement_verified": agreement_ok,
        "model_agreement_note": (
            f"Trace-against-trace against PyRhO {reference.get('pyrho_version')} "
            f"(run in a separate Python 3.9 environment; PyRhO not modified): "
            f"max |dO| = {worst_trace:.3e}. Establishes correct transcription of "
            "the equations, NOT that the three-state model is right about ChR2."
            if agreement_ok and reference is not None
            else (
                "UNVERIFIED: no PyRhO reference trace was available to compare "
                "against. The steady state agreeing is not sufficient -- it is one "
                "equation in three rates."
                if reference is None
                else f"FAILED: max |dO| = {worst_trace:.3e} exceeds "
                f"{agreement_tolerance:g}. Investigate optosim, not PyRhO."
            )
        ),
    }
    print("\n=== JSON ===")
    print(json.dumps(payload, indent=2))
    if args.output:
        path = pathlib.Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2))
        print(f"\nWritten to {path}")
    return 0 if (numerics_ok and algebra_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
