"""The optical stimulus: irradiance as a function of time.

Scope, stated plainly because the omission matters
--------------------------------------------------
Irradiance here is the irradiance **at the opsin**. This module does not model
light propagation through tissue: no absorption, no scattering, no
depth-dependent attenuation, no geometry of the fibre or the illuminated volume.

That is not an oversight to be inferred from the absence of a function. A real
preparation illuminated at the surface with 10 mW/mm^2 delivers very much less
than that a millimetre down -- brain tissue attenuates blue light steeply -- so
an irradiance set here is what a cell *receives*, and relating it to what a light
source *emits* is a layer that does not exist yet.

:class:`LightField` is the seam for it. A tissue model would implement the same
``irradiance_at`` and take depth into account; nothing downstream changes. Until
one does, :meth:`IrradianceProtocol.depth_note` is what the report prints, so
every run says the attenuation is unmodelled rather than letting a reader assume
it was handled.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

from .units import CHR2_PEAK_WAVELENGTH_NM, irradiance_to_photon_flux


class LightField(Protocol):
    """What the opsin needs from the optics: irradiance at a time."""

    def irradiance_at(self, time_ms: float) -> float: ...


@dataclass(frozen=True, slots=True)
class LightPulse:
    """A rectangular pulse of constant irradiance.

    Args:
        onset_ms: When the light comes on.
        duration_ms: How long it stays on. Must be positive; a zero-length pulse
            is not a stimulus and is refused rather than silently ignored.
        irradiance_mw_mm2: Irradiance at the opsin while on.
    """

    onset_ms: float
    duration_ms: float
    irradiance_mw_mm2: float

    def __post_init__(self) -> None:
        if self.duration_ms <= 0:
            raise ValueError(
                f"pulse duration must be positive, got {self.duration_ms} ms. A "
                "zero-length pulse is refused rather than silently delivering "
                "nothing."
            )
        if self.irradiance_mw_mm2 < 0:
            raise ValueError("irradiance cannot be negative; darkness is 0.0")
        if self.onset_ms < 0:
            raise ValueError("pulse onset cannot be before the start of the run")

    @property
    def offset_ms(self) -> float:
        return self.onset_ms + self.duration_ms

    def contains(self, time_ms: float) -> bool:
        """Half-open ``[onset, offset)``, so abutting pulses do not double up."""
        return self.onset_ms <= time_ms < self.offset_ms

    def describe(self) -> str:
        return (
            f"{self.irradiance_mw_mm2:g} mW/mm^2 from {self.onset_ms:g} to "
            f"{self.offset_ms:g} ms ({self.duration_ms:g} ms)"
        )


@dataclass(frozen=True, slots=True)
class IrradianceProtocol:
    """A train of pulses, and the wavelength they are delivered at."""

    pulses: tuple[LightPulse, ...] = field(default_factory=tuple)
    wavelength_nm: float = CHR2_PEAK_WAVELENGTH_NM

    def __post_init__(self) -> None:
        if self.wavelength_nm <= 0:
            raise ValueError("wavelength must be positive")
        ordered = sorted(self.pulses, key=lambda p: p.onset_ms)
        for earlier, later in zip(ordered, ordered[1:], strict=False):
            if later.onset_ms < earlier.offset_ms:
                raise ValueError(
                    f"pulses overlap: {earlier.describe()} and {later.describe()}. "
                    "Overlapping pulses would sum to an irradiance no light "
                    "source in the protocol actually delivered, so they are "
                    "refused rather than added."
                )

    @classmethod
    def single(
        cls, onset_ms: float, duration_ms: float, irradiance_mw_mm2: float, **kwargs: float
    ) -> IrradianceProtocol:
        return cls(pulses=(LightPulse(onset_ms, duration_ms, irradiance_mw_mm2),), **kwargs)

    @classmethod
    def train(
        cls,
        n_pulses: int,
        onset_ms: float,
        duration_ms: float,
        period_ms: float,
        irradiance_mw_mm2: float,
        **kwargs: float,
    ) -> IrradianceProtocol:
        """A regular pulse train, the standard optogenetic driving protocol."""
        if n_pulses < 1:
            raise ValueError("a train needs at least one pulse")
        if period_ms < duration_ms:
            raise ValueError(
                f"period {period_ms} ms is shorter than the pulse duration "
                f"{duration_ms} ms, which would overlap"
            )
        return cls(
            pulses=tuple(
                LightPulse(onset_ms + index * period_ms, duration_ms, irradiance_mw_mm2)
                for index in range(n_pulses)
            ),
            **kwargs,
        )

    def irradiance_at(self, time_ms: float) -> float:
        for pulse in self.pulses:
            if pulse.contains(time_ms):
                return pulse.irradiance_mw_mm2
        return 0.0

    def photon_flux_at(self, time_ms: float) -> float:
        """What the opsin model actually consumes, in photons/mm^2/s."""
        return irradiance_to_photon_flux(
            self.irradiance_at(time_ms), wavelength_nm=self.wavelength_nm
        )

    @property
    def switch_times_ms(self) -> tuple[float, ...]:
        """Every instant the irradiance changes.

        An adaptive integrator must be told about these. Light steps are
        discontinuities, and a step-size controller that straddles one will
        either smear the onset or spend a great many rejected steps finding it.
        """
        edges: set[float] = set()
        for pulse in self.pulses:
            edges.add(pulse.onset_ms)
            edges.add(pulse.offset_ms)
        return tuple(sorted(edges))

    def depth_note(self) -> str:
        """The caveat every report carries while tissue optics is unmodelled."""
        return (
            "Irradiance is specified AT THE OPSIN. Propagation through tissue is "
            "not modelled: no absorption, scattering or depth attenuation. Blue "
            "light falls off steeply in brain tissue, so these values are not "
            "the output of a light source at the surface and must not be read as "
            "one. Tissue optics is a future layer, not an effect already included."
        )

    def describe(self) -> str:
        if not self.pulses:
            return "  no light: the whole run is in darkness"
        lines = [f"  {len(self.pulses)} pulse(s) at {self.wavelength_nm:g} nm:"]
        lines.extend(f"    {pulse.describe()}" for pulse in self.pulses)
        lines.append(f"  {self.depth_note()}")
        return "\n".join(lines)
