"""Partial cellular reprogramming as a coupled dynamical system.

A branch of optosim's discipline applied to a different problem. Where the
parent package drives a membrane with light, this one drives a cell's
regulatory state with reprogramming-factor expression (the Yamanaka OSK/OSKM
programme), and asks the question that governs whether partial reprogramming is
a therapy or a tumour: **can epigenetic age be lowered without the cell losing
its identity or committing to pluripotency?**

What this is, stated before anything else can be mistaken for more
---------------------------------------------------------------------
This is a **conceptual, phenomenological** model. It is not fitted to any
dataset, it does not represent any specific gene-regulatory network, and its
time axis is **uncalibrated** -- the clock carries no claim about hours, days or
cell cycles (see :mod:`optosim.reprogram.factors`). Every one of its parameters
is ``Origin.HEURISTIC`` and listed as such, because there is no published
parameter set for a model at this level of abstraction and inventing one with a
citation attached is precisely the failure optosim exists to prevent.

What it *is* good for is the thing the news coverage skips: making the
**trade-off** explicit and visible. Stronger or longer factor expression lowers
the epigenetic-age marker faster -- and also pushes the cell past the
pluripotency commitment point, where identity collapses and the teratoma-risk
index runs away. The "safe window" everyone refers to is a region in this
trade-off, not a slogan, and the model's job is to let you see it move as you
change the protocol.

Read :meth:`ReprogramReport.caveats` output before quoting any number from a run.
"""

from __future__ import annotations

from .factors import FactorProtocol, FactorPulse
from .model import ConceptualReprogramming
from .parameters import REPROGRAM_CONCEPTUAL
from .report import HAZARDS, ReprogramReport, classify_regime, regime_flags
from .simulate import ReprogramResult, simulate_reprogramming

__all__ = [
    "HAZARDS",
    "REPROGRAM_CONCEPTUAL",
    "ConceptualReprogramming",
    "FactorProtocol",
    "FactorPulse",
    "ReprogramReport",
    "ReprogramResult",
    "classify_regime",
    "regime_flags",
    "simulate_reprogramming",
]
