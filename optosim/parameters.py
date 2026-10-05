"""Model parameters, each carrying where it came from.

The analogue of chemdisco's curation report. There, the question was what the
dataset discarded; here it is which of the numbers driving a simulation anyone
actually measured.

Three categories, and the middle one is the one that matters:

* **Cited** -- a value read from a paper or a published table. ``Origin.MEASURED``
  with a source precise enough to find it again. Note the limit of that claim: a
  constant fitted by its authors to *their* preparation is traceable to published
  work, not measured in the cell being simulated.
* **Uncited** -- a value of the right order of magnitude taken from nobody in
  particular, needed to make a model run. ``Origin.HEURISTIC``. Legitimate, and
  listed explicitly by :func:`optosim.report.parameter_report` every time.
* **Absent** -- ``value=None``. Never a silent default.

Nothing here is allowed to be a bare float. A parameter set is queried through
:meth:`ParameterSet.value_of`, which raises on an unknown name rather than
returning a plausible number, because the failure this package most needs to
avoid is a simulation that runs to completion on an invented constant.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field

from .provenance import Origin, ProvenanceError, Quantity
from .units import normalise_unit

#: Where the ChR2 three-state numbers in this file were read from. PyRhO is the
#: reference implementation for opsin photocycle models; this is the direct
#: source for the values, quoted to the file and table so it can be checked.
PYRHO_SOURCE = (
    "PyRhO 0.9.4, parameters.py modelFits['3']['ChR2'] (lines 106-117); "
    "package: Evans et al. (2016) Front Neuroinform 10:8. The three-state "
    "photocycle itself is from Nikolic et al. (2009) Photochem Photobiol "
    "85:400-411 -- named as model provenance, not as a table this package "
    "verified independently"
)

#: Where the Hodgkin-Huxley numbers came from.
HH_SOURCE = (
    "Hodgkin & Huxley (1952) J Physiol 117:500-544, in the modern "
    "resting-potential convention (V_rest ~ -65 mV) as tabulated in standard "
    "treatments, e.g. Dayan & Abbott (2001) Theoretical Neuroscience ch. 5"
)


@dataclass(frozen=True, slots=True)
class Parameter:
    """One model constant, its provenance, and what it means.

    Attributes:
        name: Identifier used in the equations, e.g. ``"g_Na"``.
        quantity: The value and its :class:`~optosim.provenance.Origin`.
        description: What it is, in words, for the report.
    """

    name: str
    quantity: Quantity
    description: str

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ProvenanceError("a parameter must have a name")
        normalise_unit(self.quantity.unit)

    @property
    def is_cited(self) -> bool:
        """Whether anyone is on record for this number."""
        return self.quantity.origin is Origin.MEASURED

    def describe(self) -> str:
        line = f"  {self.name:<12} {self.quantity.label(digits=4):<34} {self.description}"
        if not self.is_cited:
            line += "\n" + " " * 16 + f"UNCITED ({self.quantity.origin.value}): {self.quantity.source}"
        elif self.quantity.source:
            line += "\n" + " " * 16 + f"source: {self.quantity.source}"
        return line


def cited(
    name: str, value: float, unit: str, description: str, source: str, *, notes: tuple[str, ...] = ()
) -> Parameter:
    return Parameter(
        name=name,
        quantity=Quantity.measured(value, normalise_unit(unit), source, notes=notes),
        description=description,
    )


def uncited(
    name: str, value: float, unit: str, description: str, why: str, *, notes: tuple[str, ...] = ()
) -> Parameter:
    """A parameter with no citation. ``why`` must say what is missing."""
    if not why.strip():
        raise ProvenanceError(
            "an uncited parameter must say what is missing; that string is the "
            "only thing standing between it and being read as sourced"
        )
    return Parameter(
        name=name,
        quantity=Quantity.heuristic(value, normalise_unit(unit), why, notes=notes),
        description=description,
    )


@dataclass(frozen=True, slots=True)
class ParameterSet(Mapping[str, Parameter]):
    """A named, immutable collection of parameters that reports its own gaps."""

    name: str
    model: str
    parameters: tuple[Parameter, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        seen: set[str] = set()
        for parameter in self.parameters:
            if parameter.name in seen:
                raise ProvenanceError(f"duplicate parameter {parameter.name!r}")
            seen.add(parameter.name)

    def __getitem__(self, key: str) -> Parameter:
        for parameter in self.parameters:
            if parameter.name == key:
                return parameter
        raise KeyError(
            f"no parameter {key!r} in {self.name}. Available: "
            f"{', '.join(p.name for p in self.parameters)}"
        )

    def __iter__(self) -> Iterator[str]:
        return (p.name for p in self.parameters)

    def __len__(self) -> int:
        return len(self.parameters)

    def value_of(self, key: str) -> float:
        """The number, or an exception. Never a substitute.

        Equations call this. An unknown name is a programming error and a
        parameter whose value is ``None`` is an unanswered question; both raise,
        because the alternative is a simulation that completes on a number
        nobody chose.
        """
        return self[key].quantity.require()

    @property
    def uncited(self) -> tuple[Parameter, ...]:
        """Every parameter nobody is on record for."""
        return tuple(p for p in self.parameters if not p.is_cited)

    def with_overrides(self, **values: float) -> ParameterSet:
        """Replace values, degrading provenance to HEURISTIC.

        An overridden parameter is no longer the published one, so it stops
        claiming the citation. This is the provenance-degradation rule from
        chemdisco applied to parameters: a value can lose its source by being
        changed, and can never gain one.
        """
        updated: list[Parameter] = []
        for parameter in self.parameters:
            if parameter.name not in values:
                updated.append(parameter)
                continue
            updated.append(
                Parameter(
                    name=parameter.name,
                    quantity=Quantity.heuristic(
                        values[parameter.name],
                        parameter.quantity.unit,
                        f"overridden by caller; original was "
                        f"{parameter.quantity.value} from: {parameter.quantity.source}",
                    ),
                    description=parameter.description,
                )
            )
        return ParameterSet(name=f"{self.name} (overridden)", model=self.model, parameters=tuple(updated))


# -- Hodgkin-Huxley, squid giant axon --------------------------------------

HH_SQUID_AXON = ParameterSet(
    name="Hodgkin-Huxley squid giant axon",
    model="hh",
    parameters=(
        cited("C_m", 1.0, "uF/cm^2", "membrane capacitance per unit area", HH_SOURCE),
        cited("g_Na", 120.0, "mS/cm^2", "peak sodium conductance density", HH_SOURCE),
        cited("g_K", 36.0, "mS/cm^2", "peak potassium conductance density", HH_SOURCE),
        cited("g_L", 0.3, "mS/cm^2", "leak conductance density", HH_SOURCE),
        cited("E_Na", 50.0, "mV", "sodium reversal potential", HH_SOURCE),
        cited("E_K", -77.0, "mV", "potassium reversal potential", HH_SOURCE),
        cited(
            "E_L",
            -54.387,
            "mV",
            "leak reversal potential",
            HH_SOURCE,
            notes=(
                "chosen in the source convention so that the resting potential "
                "comes out at -65 mV; it is a derived bookkeeping value rather "
                "than an independently measured reversal potential",
            ),
        ),
        cited("V_rest", -65.0, "mV", "resting potential this parameter set settles at", HH_SOURCE),
    ),
)


# -- Leaky integrate-and-fire, the pessimistic baseline --------------------
#
# Paired with HH the way chemdisco pairs a scaffold split with a random split:
# not because it is better, but because the gap between them says how much of
# the result is coming from the biophysical detail rather than from the setup.
#
# Almost every value here is uncited, and that is the honest state of it. There
# is no canonical LIF parameter set the way there is a canonical HH one -- LIF is
# a phenomenological model whose parameters are chosen per application.

LIF_BASELINE = ParameterSet(
    name="Leaky integrate-and-fire baseline",
    model="lif",
    parameters=(
        cited(
            "C_m",
            1.0,
            "uF/cm^2",
            "membrane capacitance per unit area",
            HH_SOURCE,
            notes=("taken from the HH set so the two models share a capacitance",),
        ),
        cited(
            "g_L",
            0.3,
            "mS/cm^2",
            "leak conductance density",
            HH_SOURCE,
            notes=("taken from the HH set so the comparison differs in spiking mechanism, not in leak",),
        ),
        cited(
            "E_L",
            -54.387,
            "mV",
            "leak reversal potential",
            HH_SOURCE,
            notes=("shared with the HH set for the same reason",),
        ),
        uncited(
            "V_threshold",
            -50.0,
            "mV",
            "voltage at which a spike is declared",
            "no source: LIF has no canonical threshold, and this one is chosen "
            "to sit above the HH set's -65 mV rest without being tuned to "
            "reproduce HH's rheobase. Any comparison of spike COUNTS between the "
            "two models is therefore a comparison against an arbitrary choice",
        ),
        uncited(
            "V_reset",
            -65.0,
            "mV",
            "voltage the membrane is reset to after a spike",
            "no source: set to the HH resting potential for want of a better "
            "reason",
        ),
        uncited(
            "t_refractory",
            2.0,
            "ms",
            "absolute refractory period",
            "no source: 2 ms is the order of magnitude usually quoted for "
            "cortical neurons, from nobody in particular. It caps the model's "
            "firing rate at 500 Hz, which is an artefact of this choice and not "
            "a property of anything measured",
        ),
    ),
)


# -- ChR2, three-state photocycle -----------------------------------------
#
# Every value read from PyRhO's own fitted table, quoted to the line. The model
# is C -> O -> D -> C with light-driven activation and light-sensitive recovery.

CHR2_3STATE = ParameterSet(
    name="ChR2 three-state photocycle (PyRhO fit)",
    model="chr2_3state",
    parameters=(
        cited("g0", 1.57e5, "", "maximum conductance, pS (absolute, whole-cell)", PYRHO_SOURCE,
              notes=("unit is pS; carried dimensionless here because the package's unit "
                     "registry works in nS, and chr2.py converts once, explicitly",)),
        cited("phi_m", 5e17, "photons/mm^2/s", "Hill constant for photoactivation", PYRHO_SOURCE),
        cited("k_a", 5.0, "1/ms", "maximum photoactivation rate C->O", PYRHO_SOURCE),
        cited("k_r", 0.1, "1/ms", "maximum light-dependent recovery rate D->C", PYRHO_SOURCE),
        cited("p", 0.8, "", "Hill coefficient for photoactivation", PYRHO_SOURCE),
        cited("q", 0.25, "", "Hill coefficient for light-dependent recovery", PYRHO_SOURCE),
        cited("Gd", 0.104, "1/ms", "open-state decay rate O->D (light independent)", PYRHO_SOURCE),
        cited("Gr0", 2e-4, "1/ms", "dark recovery rate D->C", PYRHO_SOURCE),
        cited("E", 0.0, "mV", "ChR2 reversal potential", PYRHO_SOURCE),
        cited("v0", 43.0, "mV", "voltage rectification scale", PYRHO_SOURCE),
        cited("v1", 17.1, "mV", "voltage rectification amplitude", PYRHO_SOURCE),
    ),
)

# -- Compartment geometry -------------------------------------------------
#
# The bridge between the opsin's absolute conductance (nS) and HH's per-area
# currents (uA/cm^2). It scales photocurrent density LINEARLY, so every spike
# count in this package moves with it, which is why it is derived here rather
# than chosen.
#
# The derivation. Whole-cell capacitance is one of the most routinely reported
# numbers in patch-clamp electrophysiology, and specific membrane capacitance is
# close to 1 uF/cm^2 across cell types -- that near-constancy is itself a
# long-standing experimental result and is what makes the inversion possible:
#
#     area = C_whole / C_m
#
# With C_m = 1 uF/cm^2 from the Hodgkin-Huxley set and a whole-cell capacitance
# of 100 pF, a value squarely in the range reported for cortical pyramidal
# neurons, the area is 1e-4 cm^2 = 100 um^2 * 100 = 1e4 um^2.
#
# Worth being precise about what that last step means: 100 pF / (1 uF/cm^2)
# = 1e-10 F / 1e-6 F/cm^2 = 1e-4 cm^2. The arithmetic is exact; the input is a
# representative value, not a measurement of any particular cell.
#
# So the number is unchanged -- the previous uncited value happened to be right --
# but it is now DERIVED from two cited quantities rather than chosen to make the
# output look plausible. That distinction is the whole point: the value did not
# move, its justification did.

#: Specific membrane capacitance, the quantity whose near-constancy across cell
#: types is what lets a whole-cell capacitance be inverted into an area.
SPECIFIC_CAPACITANCE_SOURCE = (
    "specific membrane capacitance ~1 uF/cm^2, the value used by the "
    "Hodgkin-Huxley set and approximately conserved across cell types; see "
    "Gentet, Stuart & Clements (2000) Biophys J 79:314-320 for a direct "
    "measurement in neurons. Named as the basis for the inversion"
)

#: Whole-cell capacitance the default area is derived from.
WHOLE_CELL_CAPACITANCE_SOURCE = (
    "whole-cell capacitance 100 pF, a representative value for a cortical "
    "pyramidal neuron as routinely reported in patch-clamp recordings. A "
    "REPRESENTATIVE value, not a measurement of the cell being simulated: "
    "reported capacitances span roughly 50-300 pF and the area scales with it "
    "directly"
)

#: Whole-cell capacitance in pF, and the specific capacitance in uF/cm^2, that
#: :data:`DEFAULT_COMPARTMENT` inverts. Exposed so the sensitivity line in the
#: report can restate the derivation.
DEFAULT_WHOLE_CELL_CAPACITANCE_PF = 100.0
DEFAULT_SPECIFIC_CAPACITANCE_UF_CM2 = 1.0


def area_from_capacitance(
    whole_cell_pf: float, specific_uf_cm2: float = DEFAULT_SPECIFIC_CAPACITANCE_UF_CM2
) -> float:
    """Membrane area in cm^2 from a whole-cell capacitance in pF.

    ``area = C_whole / C_m``. pF -> F is 1e-12 and uF/cm^2 -> F/cm^2 is 1e-6, so
    the ratio carries a factor of 1e-6.

    This is the honest way to set the compartment size: whole-cell capacitance is
    measured as a matter of routine in any patch-clamp recording, and specific
    capacitance is close to 1 uF/cm^2 across cell types, so the inversion turns a
    reported number into a geometry instead of requiring one to be invented.

    Raises:
        ValueError: on a non-positive input, which has no geometric reading.
    """
    if whole_cell_pf <= 0:
        raise ValueError(f"whole-cell capacitance must be positive, got {whole_cell_pf} pF")
    if specific_uf_cm2 <= 0:
        raise ValueError(f"specific capacitance must be positive, got {specific_uf_cm2}")
    return (whole_cell_pf * 1e-12) / (specific_uf_cm2 * 1e-6)


def compartment_from_capacitance(
    whole_cell_pf: float = DEFAULT_WHOLE_CELL_CAPACITANCE_PF,
    *,
    source: str = "",
) -> ParameterSet:
    """A geometry parameter set derived from a whole-cell capacitance.

    Args:
        whole_cell_pf: Measured or representative whole-cell capacitance.
        source: Where that capacitance came from. Supply it when using a real
            measurement; the default names the representative value instead, and
            says that is what it is.
    """
    area = area_from_capacitance(whole_cell_pf)
    citation = source.strip() or WHOLE_CELL_CAPACITANCE_SOURCE
    return ParameterSet(
        name=f"Single compartment, area from {whole_cell_pf:g} pF whole-cell capacitance",
        model="geometry",
        parameters=(
            cited(
                "area",
                area,
                "cm^2",
                "membrane area, derived as C_whole / C_m",
                f"derived: {whole_cell_pf:g} pF / "
                f"{DEFAULT_SPECIFIC_CAPACITANCE_UF_CM2:g} uF/cm^2. "
                f"Capacitance: {citation}. Specific capacitance: "
                f"{SPECIFIC_CAPACITANCE_SOURCE}",
                notes=(
                    "DERIVED from two cited quantities, not measured directly. The "
                    "whole-cell capacitance is representative unless a source was "
                    "supplied, and the area scales with it linearly",
                    "photocurrent density scales linearly with 1/area, so every "
                    "spike count in a run depends on this value",
                ),
            ),
        ),
    )


#: The default compartment, derived rather than chosen. See above for the
#: arithmetic and for what "derived" is and is not claiming.
DEFAULT_COMPARTMENT = compartment_from_capacitance()
