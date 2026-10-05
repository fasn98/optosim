# Handoff

What is pending and what is unverified. `README.md` carries the method and the
measurements; this file is only the open work.

Repository state at handoff: core complete, 108 tests passing, agreement with
PyRhO measured trace against trace, and the membrane area derived from a cited
capacitance. The open items are scientific rather than structural.

## Run these first

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -q                                      # expect 108 passed
python scripts/validate_opsin.py --irradiance 1    # numerics, algebra, agreement
python scripts/validate_opsin.py --irradiance 10   # all three again at 10x flux
```

## How to work in this repository

### Detach anything long, and commit each result as it lands

**Anything expected to take more than ~10 minutes runs under `setsid`/`nohup` or
`tmux`. Commit each result the moment its file closes — not at the end of the
block of work.**

This is not a style preference. It was learned twice, the hard way, in a single
afternoon:

- A `discover` run in the sibling chemdisco repo was killed by a session teardown
  at 30 of 90 ligands, having written no output. Forty minutes of docking, gone,
  with nothing to show it had run.
- Mid-edit source changes to that same repo's `discover.py` were lost the same
  way, because they were being held for one tidy commit at the end.

What survived both teardowns was exactly what had been detached and what had been
committed. A later `discover` run, launched under `setsid`, ran for 101 minutes
straight through a session ending and finished normally.

Concretely:

```bash
setsid nohup bash -c 'long thing here' >> runs/thing.console.log 2>&1 < /dev/null &
```

Then commit the artefact as soon as it exists, with its log. Two reasons it has to
be *immediately* rather than at the end:

1. **A local, unpushed commit is already enough.** The work survives in `.git`
   even if nothing reaches the remote. There is no excuse for holding an edit.
2. **A partial artefact must not be committed as though it were complete.** Wait
   for the file to *close*, not for the work block to finish — those are different
   moments, and conflating them is how a truncated measurement gets into the
   record looking like a finished one. If a run is still writing, leave it out and
   say so.

Corollary for this repository specifically: `scripts/pyrho_reference.py` runs in a
separate Python 3.9 environment and takes a while on the 200 ms traces. Its output
is committed (`runs/pyrho_reference_*.json`) precisely so nobody has to re-run it,
and so the agreement test works without PyRhO installed.

## Pending work, most important first

### 1. Validate against MEASURED photocurrents — DONE: agreement with PyRhO

Agreement with PyRhO is now **verified trace against trace**: max |dO| of
**1.739e-08** at 1 mW/mm² and **3.389e-08** at 10 mW/mm², against a 1e-3 bar. The
residual is at `odeint`'s own tolerance scale, so it is the reference's
integration error, not an equation mismatch. Details in the README.

Reference traces are committed in `runs/pyrho_reference_{1,10}mw.json` so tests
run without PyRhO, which cannot be installed alongside optosim — it needs Python
3.9 and `numpy<1.20`. Regenerate with `scripts/pyrho_reference.py` in the 3.9
environment (see its docstring for the exact pins).

**What is still open is the harder question.** Agreeing with PyRhO shows the
equations were transcribed correctly. It says nothing about whether the
three-state photocycle describes ChR2. The remaining validation is against
**measured** photocurrents — digitised published traces, or recordings — and that
has not been attempted.

Note what would make that comparison hard, and plan for it: a recording carries
the cell's own expression level and geometry, so a mismatch in absolute current
amplitude would be uninformative. The shape is the testable part — the
peak-to-plateau ratio, the rise time, and the off-decay. This model's ratios are
2.22 at 1 mW/mm² and 2.83 at 10, which sit in the reported range but have never
been compared with a specific trace.

### 2. The four-state photocycle

`OpsinModel` is the seam and nothing uses it yet. Three states already produce a
peak-then-plateau at a plausible ratio (see README), so the reason to want four
is *not* the peak: it is separate light- and dark-adapted branches, giving
bi-exponential off-kinetics and a wavelength dependence this model has no
representation for.

PyRhO's `modelFits['4']['ChR2']` has the parameters, readable from source the
same way the three-state ones were.

### 3. Tissue optics

Irradiance is the value at the opsin and `depth_note()` says so in every report.
`LightField` is the interface. A depth-attenuation layer would implement
`irradiance_at` and nothing downstream changes.

Until it exists, no number here can be related to a light source at the surface.
That is the single largest gap between this package and an experiment.

### 4. The membrane area — DONE, and the sensitivity is now measured

Derived from `C_whole / C_m` = 100 pF / 1 µF/cm² = 1e-4 cm², both cited.
`compartment_from_capacitance()` takes a measured capacitance with its own
source. The value did not change; its justification did.

The sensitivity is measured and printed in every report, and it is **regime
dependent**: near threshold a 16× area range moves the spike count from 8 to 0,
while under saturating light it does not move at all. See the README table.

What remains, and it is a reporting habit rather than code: a run near threshold
should state the regime it sat in, because there the spike count is a statement
about the assumed cell size and almost nothing else.

### 5. The baseline's uncited parameters

The LIF threshold, reset and refractory period have no source — LIF is
phenomenological and has no canonical set. The 39×-vs-1 gap in the README is
therefore partly a property of those three choices. A sensitivity sweep over them
would say how much, and would turn the gap from an illustration into a
measurement.

### 6. Temperature: no Q10 correction — a real scientific gap

The package runs **squid axon kinetics at squid axon temperature** driving a
**mammalian** opsin fit. The Hodgkin-Huxley rate expressions are used exactly as
published, with no temperature scaling, and the ChR2 fit comes from mammalian
cell recordings at a different temperature again.

Channel kinetics are strongly temperature dependent — rates typically change
severalfold over the gap between squid-axon and mammalian body temperature — so
the spike counts this package reports are not what the same protocol would
produce at 37 °C. This is not a stylistic omission to be noted in passing; it
limits what any number here can be compared with.

A Q10 factor on the HH rates is small work. The harder part is being honest about
what the reference temperatures actually were for each parameter set, and whether
scaling the HH rates while leaving the opsin fit unscaled is coherent at all. Do
not add a Q10 term that silently assumes both sets share a reference temperature.

Deliberately not fixed in the current round of work, and recorded here rather
than left to be noticed.

## Settled, so it does not get re-litigated

- **The three-state model does produce a peak.** It was claimed not to, a test
  asserting monotonicity failed, and the test was wrong. Ratios 2.22 at 1 mW/mm²
  and 2.83 at 10, from slow dark recovery. Measured, in the README.
- **−57.9 mV is a genuine equilibrium at 5 mW/mm².** Hand-verified current
  balance: ionic +15.11 against photocurrent −15.03 µA/cm². A first estimate said
  otherwise and had used the wrong `n_inf`. There is a test.
- **`nS * mV == pA` exactly.** That identity is why the package works in nS and
  converts PyRhO's pS `g0` once.
- **Randomness is not wanted.** `test_no_fabrication.py` fails the build if a
  module reaches an RNG unallowlisted. The allowlist is empty on purpose; two
  identical runs agree bit for bit. If stochastic gating is ever added, it goes on
  the allowlist with a reason, not in quietly.
- **PyRhO needs its own Python 3.9 environment** and is not installable beside
  optosim. The reference traces are committed so nothing routine depends on it.
  Do not patch PyRhO to make it run here; a deviation measured against a modified
  reference would mean nothing.
- **`provenance.py` is a verbatim port.** It will drift from chemdisco's copy.
  When a third project wants it, make it a shared dependency rather than a third
  copy.

## The standard this project holds itself to

Every output is `PREDICTED`, because that is what a simulator produces. A
constant is cited or it is marked uncited and listed. A value that is absent is
`None`, never a default. A validation that could not be run is reported as not
run, with the reason, and the reference is not edited to make it pass.
