"""What every number in a run is standing on.

chemdisco reports what each stage of curation discarded. The equivalent question
for a simulator is not what it threw away -- it threw away nothing -- but which
of the constants driving it anyone ever measured, and that is what this prints.

The report is deliberately not a summary of results. It is a list of the inputs
with their provenance, the uncited ones called out, and the modelling choices
that are not represented at all. A reader who sees a voltage trace and wants to
know how much to believe it should be able to get the answer from here.
"""

from __future__ import annotations

from .light import IrradianceProtocol
from .parameters import ParameterSet
from .provenance import Origin


def parameter_report(*sets: ParameterSet) -> str:
    """Every parameter, its provenance, and an explicit list of the uncited ones."""
    lines = ["=" * 78, "Parameters, and who is on record for them", "=" * 78]
    total = 0
    uncited: list[tuple[str, str]] = []
    for parameter_set in sets:
        lines.append(f"\n{parameter_set.name}  [model: {parameter_set.model}]")
        for name in parameter_set:
            parameter = parameter_set[name]
            lines.append(parameter.describe())
            total += 1
            if not parameter.is_cited:
                uncited.append((parameter_set.name, parameter.name))

    lines.append("\n" + "-" * 78)
    if not uncited:
        lines.append(
            f"All {total} parameters carry a citation. Note what that does and does "
            "not mean:\na constant fitted by its authors to their own preparation is "
            "traceable to\npublished work, not measured in the cell being simulated."
        )
    else:
        lines.append(
            f"{len(uncited)} of {total} parameters are UNCITED. These are chosen "
            "values of\nplausible magnitude, not measurements, and every result "
            "below depends on them:"
        )
        for set_name, parameter_name in uncited:
            lines.append(f"    {parameter_name}  ({set_name})")
        lines.append(
            "\nAn uncited parameter is not a defect to be hidden -- a model needs "
            "numbers\nbefore anyone has measured them. It is a defect to present a "
            "result as though\nthey had been."
        )
    return "\n".join(lines)


def protocol_report(protocol: IrradianceProtocol) -> str:
    """The stimulus, and what the optics does not model."""
    return "\n".join(
        [
            "=" * 78,
            "Optical stimulus",
            "=" * 78,
            protocol.describe(),
        ]
    )


def unmodelled_report() -> str:
    """Effects absent from the model, listed so they are not inferred as handled.

    chemdisco's standard: where a limitation is known it is written down rather
    than left to be worked out from an absence. These are the things a reader
    would most reasonably assume were included.
    """
    return "\n".join(
        [
            "=" * 78,
            "Not modelled -- stated rather than left to be inferred",
            "=" * 78,
            "  Tissue optics. Irradiance is the value AT THE OPSIN. No absorption,",
            "  scattering or depth attenuation. Relating these numbers to a light",
            "  source at the surface is a layer that does not exist yet.",
            "",
            "  Opsin expression. A single uniform conductance stands in for",
            "  expression level, trafficking, and membrane distribution. The",
            "  compartment area and g0 together set the photocurrent scale, and the",
            "  area is uncited.",
            "",
            "  Spatial structure. One compartment: no dendrites, no axon, no",
            "  propagation. Every current is assumed to act on one isopotential",
            "  patch.",
            "",
            "  Temperature. The Hodgkin-Huxley rates are at the source's",
            "  temperature with no Q10 correction, so this is squid axon kinetics",
            "  at squid axon temperature driving a mammalian opsin fit.",
            "",
            "  Calcium, chloride, pumps, and any current other than the three HH",
            "  conductances plus the photocurrent. No adaptation, no",
            "  after-hyperpolarisation beyond what the HH potassium current gives.",
            "",
            "  Photobleaching, and the four-state photocycle's light- and",
            "  dark-adapted branches: the three-state model has one open state, so",
            "  off-kinetics are single-exponential and wavelength dependence is",
            "  absent.",
        ]
    )


def run_report(result, *, include_parameters: bool = True) -> str:
    """The full provenance report for one :class:`~optosim.simulate.SimulationResult`."""
    blocks = [
        "=" * 78,
        f"Simulated: {result.model}",
        "=" * 78,
        f"  integrator: {result.integrator}",
        f"  window:     {result.voltage.times_ms[0]:g} to "
        f"{result.voltage.times_ms[-1]:g} ms, {len(result.voltage)} samples",
        "",
        "  Outputs, every one of them PREDICTED:",
        f"    V peak            {result.voltage.peak().label(digits=2)}",
        f"    V trough          {result.voltage.trough().label(digits=2)}",
        f"    O peak            {result.open_fraction.peak().label(digits=4)}",
        f"    photocurrent min  {result.photocurrent_pa.trough().label(digits=1)}",
        f"    spikes            {result.n_spikes.label(digits=0)}",
        f"    mean rate         {result.mean_rate_hz().label(digits=2)} Hz",
        "",
        "  A simulated spike count is not a recording. There is no measurement",
        "  anywhere in this report; the only question it answers is what the model",
        "  does given the parameters listed below.",
        "",
        protocol_report(result.protocol),
    ]
    if include_parameters:
        blocks.append(parameter_report(*result.parameter_sets))
    blocks.append(unmodelled_report())
    return "\n".join(blocks)


def model_comparison_report(hh_result, lif_result) -> str:
    """Hodgkin-Huxley against the integrate-and-fire baseline.

    The analogue of chemdisco's scaffold split beside its random split: the
    pairing exists so the gap can be read. A protocol that drives both to the
    same answer has demonstrated nothing about channel dynamics; a large gap
    says the biophysics is carrying the result -- it does not say which model is
    right, and the baseline's threshold is uncited, so the gap's sign is not on
    its own a finding.
    """
    hh_spikes = len(hh_result.spike_times_ms)
    lif_spikes = len(lif_result.spike_times_ms)
    lines = [
        "=" * 78,
        "Biophysical detail against the pessimistic baseline",
        "=" * 78,
        f"  Hodgkin-Huxley        {hh_spikes:>4} spikes",
        f"  integrate-and-fire    {lif_spikes:>4} spikes",
    ]
    if hh_spikes == lif_spikes:
        lines.append(
            "\n  The two agree. This protocol does not discriminate between them, so\n"
            "  nothing here is evidence about channel dynamics."
        )
    else:
        ratio = (
            f"{max(hh_spikes, lif_spikes) / min(hh_spikes, lif_spikes):.1f}x"
            if min(hh_spikes, lif_spikes) > 0
            else "one of them produced no spikes at all"
        )
        louder = "integrate-and-fire" if lif_spikes > hh_spikes else "Hodgkin-Huxley"
        lines.append(
            f"\n  They disagree by {ratio}, with {louder} firing more. The gap is the\n"
            "  measurement this pairing exists to produce: it says how much of the\n"
            "  result comes from the conductance dynamics rather than from the\n"
            "  stimulus. It does NOT say which model is closer to a real neuron, and\n"
            "  the baseline's threshold and refractory period are both uncited, so\n"
            "  the magnitude is partly a property of those choices."
        )
    if lif_spikes > hh_spikes:
        lines.append(
            "\n  Why this direction is expected here: the baseline has no delayed\n"
            "  rectifier and no sodium inactivation, so a sustained inward current\n"
            "  makes it integrate to threshold and reset indefinitely. The\n"
            "  Hodgkin-Huxley membrane instead settles at a depolarised equilibrium\n"
            "  where the potassium current balances the photocurrent, and stops\n"
            "  firing. A baseline that cannot represent that mechanism will always\n"
            "  overestimate spiking under strong sustained light."
        )
    return "\n".join(lines)


def provenance_audit(result) -> tuple[bool, str]:
    """Check that nothing in a result claims to be measured.

    The structural guarantee, checked rather than asserted. A simulator whose
    outputs could be labelled MEASURED would defeat the whole package, so this
    is verified on the way out and covered by a test.
    """
    problems: list[str] = []
    for trace in (result.voltage, result.open_fraction, result.photocurrent_pa):
        if trace.origin is not Origin.PREDICTED:
            problems.append(f"trace {trace.name} claims origin {trace.origin.value}")
    for quantity in (result.n_spikes, result.mean_rate_hz()):
        if quantity.origin is not Origin.PREDICTED:
            problems.append(f"summary claims origin {quantity.origin.value}")
    if problems:
        return False, "; ".join(problems)
    return True, "every output is PREDICTED, as a simulator's outputs must be"
