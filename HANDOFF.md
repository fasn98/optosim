# Handoff

What is pending and what is unverified. `README.md` carries the method and the
measurements; this file is only the open work.

Repository state at handoff: core complete, 91 tests passing, one validation
anchor measured and one explicitly not.

## Run these first

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -q                                      # expect 91 passed
python scripts/validate_opsin.py               # numerics PASS, agreement UNVERIFIED
```

## Pending work, most important first

### 1. Agreement with PyRhO is UNVERIFIED

**The top open item, and the one most likely to be mistaken for done** — the
numerics validation passes to machine precision, which looks like a green tick
and is answering a different question.

What is validated: the integrator reproduces the exact matrix-exponential
solution (5.66e-15), and the converged open fraction matches PyRhO's published
closed-form steady state (relative 1.148e-12). Both concern the *algebra and the
numerics*, which are the parts this package controls.

What is not: whether the model's **transient** agrees with PyRhO's. That is where
photocycle implementations actually diverge, and the closed form cannot see it.

Why it is not done: PyRhO 0.9.4 installs but does not import under Python 3.12.
`PyRhOparameters` subclasses lmfit's `Parameters` and calls
`OrderedDict.__setitem__` on itself; lmfit moved `Parameters` off `OrderedDict`.
Pinning lmfit to the 0.9.x series PyRhO's own source names fails because that
series does not build on 3.12, and only 3.12 is on this machine.

**PyRhO was deliberately not patched.** A deviation measured against a modified
reference means nothing. Options, in the order I would try them:

1. A Python 3.9/3.10 environment (conda, docker, pyenv) where PyRhO's pinned
   dependency set installs. Cleanest, and leaves the reference untouched.
2. Digitised published ChR2 photocurrent traces as the reference instead. Also
   leaves PyRhO alone, and is arguably a better anchor — it compares against a
   recording rather than another simulator.
3. Reimplement PyRhO's integration independently and compare. Weakest: two
   implementations of the same equations agreeing says little.

Report the deviation whatever it is. If it diverges, investigate before adjusting
anything — the README records one case where a confident expectation about this
model was simply the wrong sign.

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

### 4. The membrane area is uncited and matters most

`DEFAULT_COMPARTMENT["area"]` is `1e-4 cm^2`, chosen so PyRhO's whole-cell
conductance gives photocurrent densities of a plausible size against HH's
`uA/cm^2`. It scales photocurrent density **linearly**, so every spike count in
the README moves with it.

Either source it from a real preparation, or make the sweep explicit and report
spike counts as a function of it rather than at one arbitrary value.

### 5. The baseline's uncited parameters

The LIF threshold, reset and refractory period have no source — LIF is
phenomenological and has no canonical set. The 39×-vs-1 gap in the README is
therefore partly a property of those three choices. A sensitivity sweep over them
would say how much, and would turn the gap from an illustration into a
measurement.

### 6. Temperature

HH rates are at the source's temperature with no Q10 correction, so the package
currently runs squid axon kinetics against a mammalian opsin fit. A Q10 term is
small work and would remove an obvious objection.

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
- **`provenance.py` is a verbatim port.** It will drift from chemdisco's copy.
  When a third project wants it, make it a shared dependency rather than a third
  copy.

## The standard this project holds itself to

Every output is `PREDICTED`, because that is what a simulator produces. A
constant is cited or it is marked uncited and listed. A value that is absent is
`None`, never a default. A validation that could not be run is reported as not
run, with the reason, and the reference is not edited to make it pass.
