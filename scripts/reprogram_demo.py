#!/usr/bin/env python3
"""The conceptual reprogramming branch, exercised and reported.

Three things, in order:

1. **The integration is checked against the exact solution** for the one part of
   the model that has one -- ``P`` under constant induction. Same move as
   ``scripts/validate_opsin.py``: it validates the integrator, and nothing about
   whether the model is right about biology.

2. **A matched-exposure protocol comparison.** Every protocol delivers the same
   total induction time and is followed by the same recovery window, because
   without both of those the comparison measures the bookkeeping instead of the
   protocol. The first version of this script had the cyclic arms ending at the
   edge of the run, so they were scored with no time to recover and looked worse
   than a sustained pulse that had 200 units of quiet afterwards. That was an
   artefact of the schedule, not a property of cyclic induction.

3. **A drive sweep**, which is where the trade-off becomes a picture: stronger
   induction removes more of the age marker and also pushes the cell past the
   commitment threshold.

Nothing here is a measurement, and the time axis is uncalibrated. Run it and
read the caveat block it prints.

    python scripts/reprogram_demo.py --output runs/reprogram_demo.json
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from optosim.report import parameter_report  # noqa: E402
from optosim.reprogram import (  # noqa: E402
    ConceptualReprogramming,
    FactorProtocol,
    ReprogramReport,
    classify_regime,
    regime_flags,
    simulate_reprogramming,
)

#: Model-time after the last withdrawal. Identical for every protocol compared,
#: so recovery is measured rather than accidentally granted to one arm.
RECOVERY_WINDOW = 150.0
ONSET = 10.0


def check_exact_step_response(expression: float = 1.0, dt_ms: float = 0.05) -> dict[str, float]:
    """Compare the integrated ``P`` against its closed form under constant drive.

    With ``f`` held, ``dP/dt = k_on c f (1-P) - k_off P`` is linear with constant
    coefficients, so ``P(t)`` is a single exponential approach to ``P_inf``. The
    integrator has no excuse for disagreeing.
    """
    model = ConceptualReprogramming()
    duration = 300.0
    result = simulate_reprogramming(
        FactorProtocol.single(0.0, duration, expression), duration_ms=duration, dt_ms=dt_ms
    )
    worst = 0.0
    worst_at = 0.0
    for time_ms, numeric in zip(result.pluripotency.times_ms, result.pluripotency.values, strict=True):
        exact = model.pluripotency_step_response(expression, 0.0, time_ms)
        error = abs(numeric - exact)
        if error > worst:
            worst, worst_at = error, time_ms
    p_inf = model.pluripotency_steady_state(expression)
    return {
        "expression": expression,
        "dt": dt_ms,
        "max_abs_error": worst,
        "at_model_time": worst_at,
        "p_inf": p_inf,
        "final_p": result.pluripotency.values[-1],
        "steady_state_residual": abs(result.pluripotency.values[-1] - p_inf),
    }


def _run(protocol: FactorProtocol, span: float) -> object:
    return simulate_reprogramming(
        protocol, duration_ms=ONSET + span + RECOVERY_WINDOW, dt_ms=0.05
    )


def matched_exposure_comparison(exposure: float = 40.0, expression: float = 1.0) -> list[dict]:
    """Sustained against cyclic at equal total induction time and equal recovery."""
    arms: list[tuple[str, FactorProtocol, float]] = [
        ("sustained", FactorProtocol.single(ONSET, exposure, expression), exposure)
    ]
    for n_pulses in (2, 4, 8):
        pulse = exposure / n_pulses
        period = pulse * 3.0  # duty cycle 1/3
        span = period * (n_pulses - 1) + pulse
        arms.append(
            (
                f"cyclic x{n_pulses} (duty 1/3)",
                FactorProtocol.train(n_pulses, ONSET, pulse, period, expression),
                span,
            )
        )

    rows = []
    for name, protocol, span in arms:
        result = _run(protocol, span)
        rows.append(
            {
                "protocol": name,
                "total_induction": exposure,
                "schedule_span": span,
                "rejuvenation": result.rejuvenation().value,
                "peak_pluripotency": result.peak_pluripotency().value,
                "identity_floor": result.identity_floor().value,
                "identity_final": result.identity_recovered().value,
                "risk_index": result.final_risk().value,
                "regime": classify_regime(result),
                "hazards": list(regime_flags(result)),
            }
        )
    return rows


def drive_sweep(exposure: float = 40.0) -> list[dict]:
    """The trade-off, as one variable: expression against outcome."""
    rows = []
    for expression in (0.1, 0.2, 0.3, 0.5, 0.75, 1.0):
        result = _run(FactorProtocol.single(ONSET, exposure, expression), exposure)
        rows.append(
            {
                "expression": expression,
                "rejuvenation": result.rejuvenation().value,
                "peak_pluripotency": result.peak_pluripotency().value,
                "identity_floor": result.identity_floor().value,
                "risk_index": result.final_risk().value,
                "regime": classify_regime(result),
            }
        )
    return rows


def potency_sweep(exposure: float = 40.0, expression: float = 0.3) -> list[dict]:
    """What a MORE POTENT factor does to a protocol that was previously benign.

    The question the news coverage raises and does not answer. Note hard what
    this is not: ``potency`` is a multiplier on the drive in these equations,
    and a published fold-change in marker expression is a measurement of a
    different quantity. The sweep shows what the MODEL says happens if drive
    rises, which is a statement about the equations only.
    """
    rows = []
    for potency in (1.0, 2.0, 5.0, 10.0):
        model = ConceptualReprogramming(
            parameters=ConceptualReprogramming().parameters.with_overrides(potency=potency)
        )
        result = simulate_reprogramming(
            FactorProtocol.single(ONSET, exposure, expression),
            duration_ms=ONSET + exposure + RECOVERY_WINDOW,
            dt_ms=0.05,
            model=model,
        )
        rows.append(
            {
                "potency": potency,
                "expression": expression,
                "rejuvenation": result.rejuvenation().value,
                "peak_pluripotency": result.peak_pluripotency().value,
                "identity_floor": result.identity_floor().value,
                "risk_index": result.final_risk().value,
                "regime": classify_regime(result),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, default=None, help="write results as JSON")
    parser.add_argument("--exposure", type=float, default=40.0, help="total induction, model-time")
    args = parser.parse_args()

    print(parameter_report(ConceptualReprogramming().parameters))

    print("\n" + "=" * 78)
    print("1. The integrator against the exact step response for P")
    print("=" * 78)
    checks = [check_exact_step_response(f) for f in (0.25, 1.0)]
    for check in checks:
        print(
            f"  f={check['expression']:<5} max |P_numeric - P_exact| = "
            f"{check['max_abs_error']:.3e} at t={check['at_model_time']:g}; "
            f"P_inf={check['p_inf']:.6f}, distance from P_inf at the end = "
            f"{check['steady_state_residual']:.3e}"
        )
    print("  The agreement column is the integration error. The last column is NOT an")
    print("  error: it is how far the exponential has genuinely still to run after 300")
    print("  units at that drive, and the exact solution is the same distance away.")
    print("  Validates the integration of P. Says nothing about the biology.")

    print("\n" + "=" * 78)
    print(f"2. Matched exposure: {args.exposure:g} units of induction, "
          f"{RECOVERY_WINDOW:g} units of recovery for every arm")
    print("=" * 78)
    comparison = matched_exposure_comparison(args.exposure)
    print(f"  {'protocol':>22} {'rejuv':>7} {'peakP':>7} {'Imin':>6} {'Ifin':>6} {'risk':>8}  regime")
    for row in comparison:
        print(
            f"  {row['protocol']:>22} {row['rejuvenation']:>7.4f} {row['peak_pluripotency']:>7.4f} "
            f"{row['identity_floor']:>6.3f} {row['identity_final']:>6.3f} "
            f"{row['risk_index']:>8.3f}  {row['regime']}"
        )

    print("\n" + "=" * 78)
    print("3. Drive sweep: the trade-off as one variable")
    print("=" * 78)
    sweep = drive_sweep(args.exposure)
    print(f"  {'f':>5} {'rejuv':>7} {'peakP':>7} {'Imin':>6} {'risk':>8}  regime")
    for row in sweep:
        print(
            f"  {row['expression']:>5} {row['rejuvenation']:>7.4f} {row['peak_pluripotency']:>7.4f} "
            f"{row['identity_floor']:>6.3f} {row['risk_index']:>8.3f}  {row['regime']}"
        )

    print("\n" + "=" * 78)
    print("4. Potency sweep: a protocol that was benign, with the drive multiplied")
    print("=" * 78)
    potency = potency_sweep(args.exposure)
    print(f"  {'potency':>7} {'rejuv':>7} {'peakP':>7} {'Imin':>6} {'risk':>8}  regime")
    for row in potency:
        print(
            f"  {row['potency']:>7} {row['rejuvenation']:>7.4f} {row['peak_pluripotency']:>7.4f} "
            f"{row['identity_floor']:>6.3f} {row['risk_index']:>8.3f}  {row['regime']}"
        )
    print("  'potency' multiplies the drive in these equations. A published fold-change")
    print("  in marker expression measures a different quantity and is not this number.")

    reference = _run(FactorProtocol.single(ONSET, args.exposure, 1.0), args.exposure)
    print("\n" + ReprogramReport(reference).full())

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(
                {
                    "model": reference.model,
                    "integrator": reference.integrator,
                    "uncalibrated_time": True,
                    "all_parameters_uncited": True,
                    "exact_step_response_checks": checks,
                    "matched_exposure_comparison": comparison,
                    "drive_sweep": sweep,
                    "potency_sweep": potency,
                },
                indent=2,
            )
        )
        print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
