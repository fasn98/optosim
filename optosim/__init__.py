"""optosim: a biophysical optogenetics simulator with honest provenance.

A simulator produces predictions, by definition. The value of this package is
not that its numbers are right -- it cannot know that -- but that nothing it
emits can be mistaken for a measurement.

See :mod:`optosim.provenance` for the mechanism and :mod:`optosim.parameters`
for which constants anyone is actually on record for.
"""

from .provenance import Origin, ProvenanceError, Quantity, weakest

__all__ = ["Origin", "ProvenanceError", "Quantity", "weakest"]
