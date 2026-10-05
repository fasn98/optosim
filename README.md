# optosim

A biophysical optogenetics simulator: a single-compartment membrane, a
channelrhodopsin photocycle, and a light protocol driving them — built so that
**nothing it emits can be mistaken for a measurement.**

## Why this exists

A simulator produces predictions. That is not a shortcoming to be apologised
for; it is the definition. The failure mode worth engineering against is
different: a simulated trace, plotted on the same axes as a recording and
tabulated with the same number of decimal places, quietly acquires the authority
of data. Six months later nobody remembers which curve was measured.

So every quantity here carries its provenance, and the provenance of a
simulation output is always `PREDICTED`. There is no code path by which a
computed series can claim otherwise — `Trace` raises on any other origin at
construction, and a test asserts it.

The discipline is inherited from [chemdisco](https://github.com/fasn98/chemdisco),
and `optosim/provenance.py` is ported from it verbatim. The emphasis is inverted,
though. chemdisco mostly handled measurements and had to stop a heuristic being
laundered into one. Here essentially *everything* is predicted, so the question
is not "is this number measured?" but "which of the constants behind it is?"

## How this codebase answers that

`optosim/parameters.py` is the analogue of chemdisco's curation report. Every
model constant is either:

- **cited** — `Origin.MEASURED` with a source precise enough to find again, or
- **uncited** — `Origin.HEURISTIC`, with a *required* string saying what is
  missing, listed explicitly in every report.

`value_of()` raises on an unknown name and on a `None` value rather than
returning a plausible number, because a simulation that runs to completion on an
invented constant is the exact failure this package exists to prevent. Overriding
a parameter *degrades* its provenance to `HEURISTIC` — a value can lose its
citation by being changed and can never gain one.

### The honest parameter count

| Set | Parameters | Uncited |
|---|---|---|
| Hodgkin-Huxley (squid axon) | 8 | 0 |
| ChR2 three-state (PyRhO fit) | 11 | 0 |
| Integrate-and-fire baseline | 6 | **3** |
| Compartment geometry | 1 | **1** |

What "cited" does **not** mean: a constant fitted by its authors to *their*
preparation is traceable to published work, not measured in the cell being
simulated. This is squid axon kinetics at squid axon temperature driving a
mammalian opsin fit, with no Q10 correction.

The four uncited values are the LIF threshold, reset and refractory period — LIF
has no canonical parameter set, it is phenomenological — and the **membrane
area**. That last one is the most consequential number in the package: it scales
photocurrent density linearly, so a different cell size changes every spike count
reported here.

## Architecture

```
optosim/provenance.py   Quantity + Origin. Ported verbatim from chemdisco.
optosim/units.py        Strict unit handling. Written in chemdisco's spirit,
                        not ported — its pIC50 conversions mean nothing here.
optosim/parameters.py   Every constant, its provenance, its citation.
optosim/neuron/hh.py    Hodgkin-Huxley, single compartment.
optosim/neuron/iaf.py   Leaky integrate-and-fire: the pessimistic baseline.
optosim/opsin/chr2.py   ChR2 three-state Markov photocycle + exact solution.
optosim/light.py        Irradiance protocols. Tissue optics NOT modelled.
optosim/simulate.py     The coupled system, integrated. Returns Traces.
optosim/report.py       What every number is standing on.
scripts/validate_opsin.py  The validation anchor.
```

### The unit trap, handled explicitly

This package straddles two conventions and that is the likeliest source of a
silent factor of a thousand:

- Hodgkin-Huxley is **per unit area**: `mS/cm^2`, `uA/cm^2`.
- Opsin conductances are published **whole-cell and absolute**: PyRhO's `g0` is
  `1.57e5 pS`.

`nS * mV == pA` exactly, so the package works in nS internally and converts `g0`
from pS once. Crossing into HH's current density requires the membrane area and
happens in one function, `opsin_current_density`, which returns `None` — not
`0.0` — when the area is unusable. A test pins the resulting current magnitude so
that reading pS as nS fails loudly.

## Validation

### The photocycle numerics: validated to machine precision

Under constant illumination all three transition rates are constant, so the
photocycle is a linear system with an exact matrix-exponential solution. The
integrator is compared against it directly — 20001 samples at 1 mW/mm²:

```
max |O_numeric - O_exact|    5.662e-15   (at t = 7.37 ms)
as a fraction of the peak    1.907e-14
max |C+O+D - 1|              5.218e-15
```

Fifteen orders of magnitude below the precision any parameter here is known to.
Integration error limits nothing this package reports.

The converged open fraction also matches the closed form PyRhO publishes for
this model — `Ga·Gr / (Gd·(Gr+Ga) + Ga·Gr)` — to a relative `1.148e-12`.

### Agreement with PyRhO: NOT VERIFIED

This is the part that did not work, and it is reported as a result rather than
omitted.

The intended anchor was PyRhO itself, run on the same light pulse, with the
deviation reported. PyRhO 0.9.4 installs and **does not import** under Python
3.12: its `PyRhOparameters` subclasses lmfit's `Parameters` and calls
`OrderedDict.__setitem__` on itself, which broke when lmfit moved `Parameters`
off `OrderedDict`. Pinning lmfit to the 0.9.x series PyRhO's own source names
does not help — that series does not build on 3.12 — and only 3.12 is available
on this machine.

**PyRhO was not patched.** Editing the reference until it agrees is not a
validation, and a deviation measured against a modified reference would mean
nothing.

So: numerics validated, **model agreement unverified**. What was salvaged without
executing it is PyRhO's parameter table, read from source and cited to file and
line, plus the closed-form steady state used above. That checks the algebra
against the reference's own stated result; it cannot catch a disagreement in the
*transient*, which is where photocycle models are most likely to differ.

Encouraging, and not evidence of agreement: the peak-to-plateau ratio this model
produces — **2.22** at 1 mW/mm², **2.83** at 10 — sits in the range reported for
ChR2.

### A correction this produced

The first version of `chr2.py` claimed three states cannot produce ChR2's
characteristic peak-then-plateau transient, and that the peak requires the
four-state model. A test asserting monotonicity failed, and **the test was what
was wrong**. Measured:

```
 1 mW/mm^2   peak O = 0.297 at 12.3 ms, plateau 0.134, ratio 2.22
10 mW/mm^2   peak O = 0.627 at  4.7 ms, plateau 0.222, ratio 2.83
```

The mechanism is slow dark recovery: `Gd` is 0.104/ms against `Gr` at 0.0002/ms
dark and ~0.1/ms lit, so the closed state empties into the open state faster than
desensitised channels return to it, and O overshoots before settling. No second
open state is needed for a peak. What four states actually buys is separate
light- and dark-adapted branches — bi-exponential off-kinetics and a wavelength
dependence this model cannot represent. Those are the reasons to want it.

### The coupling: checked physically, not just numerically

Under sustained light the HH membrane settles at a depolarised plateau, and there
the ionic currents must exactly cancel the photocurrent density. At 5 mW/mm² the
plateau sits at **−57.9 mV**:

```
I_K   = +23.45 uA/cm^2
I_Na  =  -7.29
I_L   =  -1.05
sum   = +15.11   against a photocurrent density of -15.03   ->  dV/dt = -0.08 mV/ms
```

This was hand-verified before being trusted, because a first estimate said −57.9
mV could not be an equilibrium. The estimate had used `n_inf ≈ 0.35` where the
correct value is 0.43. The code was right and the estimate was wrong — recorded
because that is the sort of check worth doing before publishing a trace, and
because the arithmetic now lives in a test.

## Biophysical detail against the baseline

The LIF model is paired with HH for the reason chemdisco pairs a scaffold split
with a random split: not because it is better — it is strictly less — but because
**the gap between them says how much of a result comes from the biophysics.**

On a 100 ms pulse at 5 mW/mm²:

```
Hodgkin-Huxley          1 spike
integrate-and-fire     39 spikes
```

A 39× gap. HH settles at the depolarised equilibrium above and stops firing; the
baseline has no delayed rectifier and no sodium inactivation, so it integrates to
threshold and resets indefinitely. **A baseline that cannot represent that
mechanism will always overestimate spiking under strong sustained light.**

The HH spike count is also *non-monotonic* in irradiance — 0, 0, 1, 2, 2, 1
across 0.01 to 5 mW/mm² — which is the approach to that equilibrium as more
sodium inactivates.

The gap does **not** say which model is closer to a real neuron, and the
baseline's threshold and refractory period are both uncited, so its magnitude is
partly a property of those choices. The report says so every time it prints the
comparison.

## What this cannot do

Stated here so it does not have to be inferred from a missing function:

- **Tissue optics.** Irradiance is the value *at the opsin*. No absorption, no
  scattering, no depth attenuation. Blue light falls off steeply in brain tissue,
  so these values are not the output of a surface light source. `LightField` is
  the seam; nothing implements it yet.
- **Opsin expression.** One uniform conductance stands in for expression level,
  trafficking and membrane distribution.
- **Spatial structure.** One isopotential compartment. No dendrites, no axon, no
  propagation.
- **Temperature.** HH rates at the source's temperature, no Q10 correction.
- **Everything other than three HH conductances plus the photocurrent.** No
  calcium, no chloride, no pumps, no adaptation.
- **Four-state photocycle.** One open state, so off-kinetics are
  single-exponential and there is no wavelength dependence.
- **Stochasticity.** Deliberately none. `test_no_fabrication.py` fails the build
  if any module reaches a random number generator without being allowlisted with
  a reason. The allowlist is empty, so two identical runs agree bit for bit.

## Installation and running

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -q
python scripts/validate_opsin.py --output runs/opsin_validation.json
```

## Licence

MIT.
