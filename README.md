# optosim

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23177396.svg)](https://doi.org/10.5281/zenodo.23177396)

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
| Compartment geometry | 1 | 0 — *derived, see below* |

What "cited" does **not** mean: a constant fitted by its authors to *their*
preparation is traceable to published work, not measured in the cell being
simulated. This is squid axon kinetics at squid axon temperature driving a
mammalian opsin fit, with no Q10 correction.

The three remaining uncited values are the LIF threshold, reset and refractory
period — LIF is phenomenological and has no canonical parameter set.

### The membrane area, and what it controls

This was the most consequential uncited number in the package, because
photocurrent density is exactly linear in 1/area. It is now **derived** from two
cited quantities rather than chosen:

```
area = C_whole / C_m = 100 pF / 1 uF/cm^2 = 1e-4 cm^2
```

The inversion works because specific membrane capacitance is close to 1 µF/cm²
across cell types, while whole-cell capacitance is measured routinely in any
patch-clamp recording. `compartment_from_capacitance()` accepts a measured
capacitance with its own source; the default states that 100 pF is
*representative* of a cortical pyramidal neuron and that reported values span
roughly 50–300 pF. **The value did not change — its justification did.**

The sensitivity is the more useful half, and it is regime dependent. Spike counts
across a 16× area range (0.25× to 4× the default):

| Irradiance | 0.25× | 0.5× | 1× | 2× | 4× | |
|---|---|---|---|---|---|---|
| 0.05 mW/mm² | 4 | 1 | 0 | 0 | 0 | the result **is** the assumption |
| 0.1 mW/mm² | 8 | 2 | 1 | 0 | 0 | |
| 0.3 mW/mm² | 6 | 3 | 2 | 1 | 0 | |
| 1.0 mW/mm² | 1 | 2 | 2 | 1 | 1 | non-monotonic |
| 5.0 mW/mm² | 1 | 1 | 1 | 1 | 1 | insensitive |

Near threshold the spike count is a statement about the assumed cell size and
almost nothing else. Under saturating light it is entirely robust to it. So
neither "the area matters" nor "it doesn't" is true on its own, and **a run
reported without saying which regime it sat in is not interpretable.** Every run
report prints this.

Note which quantity the linearity belongs to: the *density* is exactly linear in
1/area — that is just the unit conversion. The *spike count* is not, and is not
even monotonic, because a stronger sustained current settles the membrane at a
more depolarised equilibrium with more sodium inactivated. A result cannot be
rescaled to another cell size by arithmetic; it has to be re-run.

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

### Agreement with PyRhO: verified, trace against trace

| Irradiance | max \|O_optosim − O_pyrho\| | at |
|---|---|---|
| 1 mW/mm² | **1.739e-08** | 0.70 ms |
| 10 mW/mm² | **3.389e-08** | 1.76 ms |

Against a 1e-3 bar for correct transcription — five orders of margin. The
residual sits at the scale of `odeint`'s own default tolerance (~1.5e-8), so it
is the reference's integration error rather than a disagreement in the equations.
That bounds agreement at that level; it does not prove exact equality.

**Why this was needed when the steady state already agreed to 1e-12.** The steady
state is *one* equation in three rates. A transcription error shifting `Ga` and
`Gr` in compensating directions would satisfy it and surface only in the
transient. Both irradiances are checked because `Ga` and `Gr` have different Hill
exponents (p = 0.8, q = 0.25) and scale differently with flux, so a cancellation
at one would not survive a tenfold change. Both agree.

**Getting PyRhO to run took a pinned chain**, diagnosed rather than guessed:
`PyRhOparameters` needs lmfit's pre-1.0 `Parameters` → pre-1.0 lmfit imports
`numpy.dual` → removed in numpy 1.20 → numpy 1.19 is the last release building on
Python 3.9. So 3.10 fails too (the first attempt died there on `numpy.dual`). The
reference runs under micromamba + Python 3.9 + `numpy<1.20` + `lmfit<1.0` +
`pyrho`, plus `ipywidgets`, which PyRhO's GUI module imports unconditionally and
its metadata omits.

**PyRhO was not patched.** Installing a package's dependencies is not editing it.
Nothing from optosim takes part in producing the reference: `scripts/pyrho_reference.py`
uses PyRhO's parameter table, its `solveStates` right-hand side, its analytic
Jacobian, and the same `odeint` call its own simulator makes
(`simulators.py:210`). It even restates the irradiance→flux conversion locally
rather than importing optosim's, so a shared error there cannot cancel out and
leave the comparison looking clean.

**What this does not establish:** that the three-state photocycle is right about
*ChR2*. Agreeing with PyRhO means the equations were transcribed correctly.
Validating against measured photocurrents is a separate job and remains undone.
The peak-to-plateau ratio — **2.22** at 1 mW/mm², **2.83** at 10 — sits in the
range reported for ChR2, which is encouraging and is not a comparison with a
recording.

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
- **Temperature — a real scientific gap, not a stylistic one.** The HH rates are
  at the source's temperature with no Q10 correction, so this runs *squid axon
  kinetics at squid axon temperature* driving a *mammalian* opsin fit. Channel
  kinetics are strongly temperature dependent, so the spike counts here are not
  what the same protocol would produce at 37 °C. Pending in HANDOFF.
- **Everything other than three HH conductances plus the photocurrent.** No
  calcium, no chloride, no pumps, no adaptation.
- **Four-state photocycle.** One open state, so off-kinetics are
  single-exponential and there is no wavelength dependence.
- **Stochasticity.** Deliberately none. `test_no_fabrication.py` fails the build
  if any module reaches a random number generator without being allowlisted with
  a reason. The allowlist is empty, so two identical runs agree bit for bit.

## The reprogramming branch (`optosim.reprogram`)

A second model in the same discipline, on the `reprogram-sim` branch: partial
cellular reprogramming as a coupled dynamical system. Where the parent package
drives a membrane with light, this drives a cell's regulatory state with
reprogramming-factor expression (OSK/OSKM), and asks the question that decides
whether partial reprogramming is a therapy or a tumour — **can the epigenetic-age
marker be lowered without the cell losing its identity or committing to
pluripotency?**

Four lumped variables, factor expression `f(t)` in `[0, 1]` from the protocol:

```
dP/dt = k_on * potency * f * (1 - P)  -  k_off * P      pluripotency programme
dI/dt = -k_dediff * P * I  +  k_rediff * (1 - P) * (1 - I)   somatic identity
dA/dt = -k_rejuv * P * A                                 epigenetic-age marker
dR/dt =  k_risk * ( relu(P - theta_commit) + relu(theta_identity - I) )   risk index
```

The structure carries the whole argument: **`A` falls only while `P` is active,
and `P` is also what erodes `I`.** The benefit and the hazard share a cause, so
there is no setting at which rejuvenation is free. That is the one thing this
model is for.

### What it is, stated before the numbers

**Every parameter is uncited, and that is the correct state, not a gap to close
later.** There is no published parameter set for a four-variable
phenomenological model of reprogramming; attaching a citation to one of these
coefficients would launder a modelling choice into a result, which is the single
failure this package is built to prevent. `parameter_report` lists all nine,
every run, and a test asserts the set contains no cited parameter.

**The time axis is uncalibrated.** Times are carried in `ms` only because that
is what the unit registry recognises. One unit is arbitrary model-time — not a
millisecond, and the model makes no claim about hours, days or cell cycles. Only
trajectory *shape* and outcome *ordering* mean anything; the clock does not.

The inherited machinery still applies: every trace is `PREDICTED` and cannot be
constructed otherwise, the integration grid lands exactly on every induction
edge, and `test_no_fabrication.py` fails the build if any module reaches an RNG.

### Verified: the integrator, against a closed form

Under constant `f` the `P` equation is linear with constant coefficients, so it
has an exact solution, and `pluripotency_step_response` returns it:

| drive | max \|P_numeric − P_exact\| |
|---|---|
| f = 0.25 | **4.791e-14** |
| f = 1.0 | **2.050e-12** |

Same move as the opsin validation, and it establishes the same narrow thing: the
arithmetic is right. It says nothing about whether the model describes a cell.
A potency change moves the analytic anchor too, and a test checks the integrator
against it there as well — otherwise `potency` could have been wired into the
integration and not the reference.

### Measured: cyclic induction beats sustained at matched exposure

Every arm delivers **40 units of induction** and gets the **same 150-unit
recovery window** afterwards:

| protocol | rejuvenation | peak P | min I | risk |
|---|---|---|---|---|
| sustained | 0.843 | 0.718 | 0.316 | 2.334 |
| cyclic ×2 (duty 1/3) | 0.894 | 0.668 | 0.331 | 1.143 |
| cyclic ×4 | 0.914 | **0.593** | 0.334 | **1.033** |
| cyclic ×8 | **0.920** | 0.537 | 0.334 | 1.202 |

Cyclic induction reaches **more** rejuvenation at a **lower** peak pluripotency
and **less** than half the risk index. The mechanism is in the `I` equation:
identity restores in the gaps, while `A` has no recovery term and so ratchets
down. There is a test on each of those three inequalities, because they are the
claim this section makes in prose.

**Matching the recovery window is what makes the comparison valid, and getting
it wrong is easy.** The first version of the demo script ran each protocol to a
fixed end time, so the cyclic arms — whose schedules span further — were scored
almost immediately after their last pulse while the sustained arm had 200 units
of quiet to recover in. Cyclic looked worse on final identity. That was the
bookkeeping, not the biology, and the comparison now matches both exposure and
recovery.

### Measured: potency is not safety

The sweep the news coverage invites. A protocol at `f = 0.3` that crosses no
threshold at the reference drive, with the drive multiplied:

| potency | rejuvenation | peak P | min I | risk | regime |
|---|---|---|---|---|---|
| 1× | 0.609 | 0.393 | 0.587 | **0.000** | no threshold crossed |
| 2× | 0.767 | 0.588 | 0.420 | 0.000 | identity not recovered |
| 5× | 0.880 | 0.798 | 0.254 | 5.499 | committed; identity lost |
| 10× | 0.915 | 0.889 | 0.185 | 10.220 | committed; identity lost |

So in this model a more potent factor does **not** widen the safe window — it
moves the same protocol out of it, and the protocol has to be re-tuned to get
back in. That is a property of these equations and nothing else.

**What it is emphatically not.** `potency` is a multiplier on the drive term. A
published fold-change in pluripotency-*marker expression* — such as the >50×
reported for the AI-designed RetroSOX/RetroKLF variants (OpenAI/Retro
Biosciences, 2025) — is an in-vitro measurement of a different quantity, and
entering it here as a value for `potency` would be exactly the category error
this package refuses. The sweep asks "what do the equations do if drive rises",
not "what does that variant do in a cell".

### What the model cannot see

Listed because these are the mechanisms that would matter most, and their
absence bounds every number above:

- **Cell cycle, proliferation and clonal expansion.** A reprogrammed cell that
  divides is how a damaged clone grows out. There is no division here.
- **Pre-existing oncogenic mutations.** The reported interaction in which
  reprogramming dramatically accelerates cancer in an already-mutated cell has
  no representation, so this model is **structurally blind** to that pathway —
  the one a safety argument would most need to address.
- **Cell-type dependence.** Liver and intestine are where in vivo reprogramming
  is reported to go wrong first; there is one generic cell here.
- **Network bistability.** Commitment is a scalar threshold crossing, where in
  reality it is a self-sustaining network state. A threshold cannot be
  irreversible; a network state can.
- Delivery, mosaic expression, population heterogeneity, immune response,
  senescence and tissue context.

The regime label reports **every** hazard condition a run met rather than the
worst one. That is a correction: an ordered worst-wins label printed
`identity lost` for a cyclic arm that was better than a sustained arm it printed
`committed` for, because the conditions are not mutually exclusive and ranking
them hid the comparison. There is a regression test.

```bash
python scripts/reprogram_demo.py --output runs/reprogram_demo.json
```

## Installation and running

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -q
python scripts/validate_opsin.py --output runs/opsin_validation.json
```

## How to cite

Citation metadata is in [`CITATION.cff`](CITATION.cff), which GitHub renders as a
"Cite this repository" button and most reference managers read directly.

A Zenodo DOI will be minted from the first archived release and added to
`CITATION.cff` and `.zenodo.json`; until then, cite the repository and the commit
or tag you used.

**What this is honest about being.** The claim is *verified implementation and
provenance enforced by construction*, not a validated model of ChR2. Specifically:
the photocycle integrator matches the exact solution to 5.7e-15, the steady state
matches PyRhO's published closed form to 1.1e-12, and the full transient matches
PyRhO's own integration to 1.7e-08. That establishes the equations were
transcribed correctly — **it does not establish that the three-state photocycle
describes the real protein**, which would need comparison against measured
photocurrents and has not been attempted.

Every output this package emits is labelled `PREDICTED`, because that is what a
simulator produces. Cite it as an implementation with a stated validation
boundary, not as evidence about channelrhodopsin.

### Data and dependencies

All code is MIT (see [`LICENSE`](LICENSE)). The ChR2 parameters are read from
PyRhO's published fit table and cited in `optosim/parameters.py`; PyRhO itself is
not vendored. The committed reference traces under `runs/` were generated by
PyRhO 0.9.4 and are included so the agreement test runs without it — PyRhO
cannot be installed alongside this package, since it needs Python 3.9 and
`numpy<1.20`.

## Licence

MIT.
