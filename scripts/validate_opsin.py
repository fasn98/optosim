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

**Is the model right about ChR2?** NOT ANSWERED. The intended reference was
PyRhO, the established implementation of these photocycle models, run on the same
protocol. PyRhO 0.9.4 installs but does not import under Python 3.12: its
``PyRhOparameters`` subclasses lmfit's ``Parameters`` and calls
``OrderedDict.__setitem__`` on itself, which stopped working when lmfit moved
``Parameters`` off ``OrderedDict``. Pinning lmfit to the 0.9.x series its own
source names does not help, because that series does not build on 3.12.

PyRhO was NOT patched to make this run. Editing the reference until it agrees is
not a validation, and a deviation measured against a modified reference would
mean nothing.

What is available without executing it: PyRhO's parameter table, read from source
and cited in :mod:`optosim.parameters`, and its documented closed-form steady
state. Comparing this implementation's converged open fraction against that
formula is a real cross-check of the algebra against the reference's own stated
result -- it is simply not a trace-against-trace comparison, and it cannot catch
a disagreement in the transient.

So the honest status is: numerics validated, model agreement UNVERIFIED. That
distinction is the finding, and HANDOFF.md carries it as the top pending item.
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--irradiance", type=float, default=1.0, help="mW/mm^2")
    parser.add_argument("--duration", type=float, default=200.0, help="ms of light")
    parser.add_argument("--dt", type=float, default=0.01, help="ms")
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    model = ChR2ThreeState()
    flux = irradiance_to_photon_flux(args.irradiance)

    heading("0. Reference availability")
    available, detail = pyrho_status()
    print(f"  PyRhO importable: {available}")
    print(f"  detail: {detail}")
    if not available:
        print(
            "\n  The trace-against-trace comparison CANNOT be run. PyRhO was not\n"
            "  patched to make it run: editing the reference until it agrees is not\n"
            "  a validation. Model agreement is therefore UNVERIFIED, and what\n"
            "  follows validates the numerics only."
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

    heading("3. What this does and does not establish")
    print(
        "  Established: the photocycle is integrated correctly, and the steady\n"
        "  state agrees with the formula PyRhO publishes for the same model.\n"
        "\n"
        "  NOT established: that the model describes real ChR2, or that this\n"
        "  implementation agrees with PyRhO's transient. The peak-to-plateau ratio\n"
        "  this model produces (2.22 at 1 mW/mm^2, 2.83 at 10) is in the range\n"
        "  reported for ChR2, which is encouraging and is not a measurement of\n"
        "  agreement with anything.\n"
        "\n"
        "  This is the same shape as chemdisco's redocking result: the mechanics\n"
        "  are right and the harder question is left open and labelled, rather\n"
        "  than answered by a number that was not measured."
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
        "pyrho_importable": available,
        "pyrho_detail": detail,
        "model_agreement_verified": False,
        "model_agreement_note": (
            "PyRhO could not be executed under Python 3.12 and was deliberately "
            "not patched; a deviation against a modified reference would be "
            "meaningless. Trace-against-trace agreement is UNVERIFIED."
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
