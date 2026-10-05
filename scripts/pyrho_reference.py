#!/usr/bin/env python3
"""Generate a ChR2 three-state reference trace with PyRhO itself.

RUNS IN A SEPARATE PYTHON 3.9 ENVIRONMENT, not the optosim venv. PyRhO 0.9.4
cannot run on 3.12: its ``PyRhOparameters`` needs lmfit's pre-1.0 ``Parameters``,
pre-1.0 lmfit imports ``numpy.dual``, and ``numpy.dual`` was removed in numpy
1.20, which is the last version supporting Python 3.9. The chain pins the whole
environment::

    micromamba create -n pyrho39 -c conda-forge python=3.9 'numpy<1.20'
    pip install 'setuptools<81' 'lmfit<1.0' 'numpy<1.20' pyrho ipywidgets

PyRhO was NOT modified. Installing its dependencies is not patching it; the only
thing added beyond its own pins is ``ipywidgets``, which its GUI module imports
unconditionally and which is genuinely missing from its metadata.

What this writes is PyRhO's own integration of its own equations with its own
parameters, so the comparison in ``validate_opsin.py`` is a genuine
trace-against-trace test of optosim's independently written implementation.

Why trace-against-trace and not just the steady state: the steady state already
agrees to 1e-12, but it is ONE equation in three rates. A transcription error
that shifted Ga and Gr in compensating directions would satisfy it and show up
only in the transient. That is the risk this closes.

Usage (in the 3.9 environment):
    python scripts/pyrho_reference.py --irradiance 1 --output runs/pyrho_reference_1mw.json
"""

from __future__ import annotations

import argparse
import json
import pathlib
import warnings

warnings.filterwarnings("ignore")

# Planck * c, and the ChR2 peak wavelength optosim uses. Duplicated here rather
# than imported, because this script runs in a different interpreter and must not
# depend on optosim -- if it imported the code under test, a shared error in the
# unit conversion would cancel out and the comparison would prove nothing.
PLANCK_TIMES_C = 6.62607015e-34 * 299792458.0
WAVELENGTH_NM = 470.0


def irradiance_to_photon_flux(irradiance_mw_mm2: float) -> float:
    """mW/mm^2 -> photons/mm^2/s. Independent restatement, deliberately."""
    energy_per_photon_j = PLANCK_TIMES_C / (WAVELENGTH_NM * 1e-9)
    return (irradiance_mw_mm2 * 1e-3) / energy_per_photon_j


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--irradiance", type=float, default=1.0, help="mW/mm^2")
    parser.add_argument("--duration", type=float, default=200.0, help="ms")
    parser.add_argument("--dt", type=float, default=0.01, help="ms")
    parser.add_argument("--voltage", type=float, default=-70.0, help="mV, for the current")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    import numpy as np
    import pyrho
    from pyrho.models import RhO_3states
    from pyrho.parameters import modelFits

    params = modelFits["3"]["ChR2"]
    model = RhO_3states(params=params, rhoType="ChR2")

    flux = irradiance_to_photon_flux(args.irradiance)

    # Dark-adapted start, then a constant-flux step.
    #
    # `solveStates` is PyRhO's ODE right-hand side, not a solution, so it is
    # integrated exactly as PyRhO's own simulator does it -- simulators.py:210,
    #     odeint(RhO.solveStates, RhO.s0, t, args=(None,), Dfun=RhO.jacobian)
    # -- including PyRhO's analytic Jacobian. The equations, the parameters and
    # the integration call are all PyRhO's; nothing from optosim takes part.
    from scipy.integrate import odeint

    model.initStates(phi=0.0)
    model.setLight(flux)
    times = np.arange(0.0, args.duration + args.dt, args.dt)
    solution = odeint(
        model.solveStates, model.s_0, times, args=(None,), Dfun=model.jacobian
    )

    open_index = model.stateVars.index("O")
    open_fraction = [float(row[open_index]) for row in solution]

    currents: list[float] | None
    try:
        currents = [float(model.calcI(args.voltage, states=row)) for row in solution]
    except Exception as error:  # pragma: no cover - depends on PyRhO internals
        print(f"  photocurrent not retrievable: {type(error).__name__}: {error}")
        currents = None

    payload = {
        "source": "PyRhO",
        "pyrho_version": getattr(pyrho, "__version__", "unknown"),
        "numpy_version": np.__version__,
        "model": "RhO_3states",
        "rho_type": "ChR2",
        "parameters": {name: float(params[name].value) for name in params},
        "state_order": list(getattr(model, "stateVars", ["C", "O", "D"])),
        "open_state_index": open_index,
        "irradiance_mw_mm2": args.irradiance,
        "photon_flux": flux,
        "wavelength_nm": WAVELENGTH_NM,
        "duration_ms": args.duration,
        "dt_ms": args.dt,
        "voltage_mv": args.voltage,
        "times_ms": [float(t) for t in times],
        "open_fraction": open_fraction,
        "current_pa": currents,
    }

    path = pathlib.Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    print(f"PyRhO {payload['pyrho_version']} reference written to {path}")
    print(f"  {len(times)} samples, state order {payload['state_order']}, "
          f"open index {open_index}")
    print(f"  peak O = {max(open_fraction):.6f}, final O = {open_fraction[-1]:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
