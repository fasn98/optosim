"""Units for single-compartment electrophysiology, and strict conversion.

Written in the spirit of ``chemdisco/units.py`` rather than ported from it. That
module converts nM/uM/mg.mL-1 into molar and then into pIC50; none of it has any
meaning here. What transfers is the *discipline*: a unit string is either
recognised or it is refused, a conversion that cannot be made returns ``None``
rather than the number unchanged, and nothing is guessed from context.

The unit set is fixed by the convention this package integrates in, which is the
one the Hodgkin-Huxley literature uses:

===================  ==========================================================
millivolt (mV)       membrane potential, reversal potentials
millisecond (ms)     time, time constants, rate denominators
mS/cm^2              area-specific conductance (HH channel densities)
uF/cm^2              area-specific capacitance
uA/cm^2              area-specific current (HH convention)
nS                   absolute conductance (opsin, patch-clamp convention)
pA                   absolute current -- note nS * mV == pA exactly
mW/mm^2              irradiance, the optical input
photons/mm^2/s       photon flux, which is what the opsin model actually sees
1/ms                 transition rates
===================  ==========================================================

The mixed convention is deliberate and is the single most likely source of a
silent factor-of-1000 error in a package like this: HH is per-unit-area, while
opsin conductances are published as whole-cell absolute values. Rather than
pretend one convention covers both, :func:`opsin_current_density` does the
conversion explicitly and requires the cell area to do it.
"""

from __future__ import annotations

import math
from typing import Final

#: Recognised unit strings. A unit outside this set is refused, not assumed.
KNOWN_UNITS: Final[frozenset[str]] = frozenset(
    {
        "mV",
        "ms",
        "mS/cm^2",
        "uF/cm^2",
        "uA/cm^2",
        "nS",
        "pA",
        "mW/mm^2",
        "photons/mm^2/s",
        "1/ms",
        "cm^2",
        "nm",
        "",
    }
)

#: Planck constant times speed of light, J*m. Used to turn mW/mm^2 into a photon
#: flux, which is the quantity the opsin photocycle responds to.
#: CODATA 2018: h = 6.62607015e-34 J*s (exact), c = 299792458 m/s (exact).
PLANCK_TIMES_C: Final[float] = 6.62607015e-34 * 299792458.0

#: Peak absorption wavelength of ChR2, used to convert irradiance to photon flux.
#: ChR2 absorbs maximally near 470 nm; blue LEDs and lasers used in optogenetics
#: are specified at 470-473 nm.
CHR2_PEAK_WAVELENGTH_NM: Final[float] = 470.0


class UnitError(ValueError):
    """Raised when a unit is unrecognised or a conversion is not defined."""


def normalise_unit(unit: str | None) -> str:
    """Canonicalise a unit string, refusing anything unrecognised.

    Refusing is the point. A simulator that accepts ``"mv"`` and silently treats
    it as millivolts will accept ``"V"`` the same way, and then every voltage is
    out by a thousand with no error anywhere.
    """
    if unit is None:
        return ""
    text = unit.strip().replace("**", "^").replace("·", "*").replace("μ", "u")
    aliases = {
        "millivolt": "mV",
        "mv": "mV",
        "millisecond": "ms",
        "msec": "ms",
        "mS/cm2": "mS/cm^2",
        "uF/cm2": "uF/cm^2",
        "uA/cm2": "uA/cm^2",
        "mW/mm2": "mW/mm^2",
        "photons/mm2/s": "photons/mm^2/s",
        "/ms": "1/ms",
        "ms^-1": "1/ms",
        "cm2": "cm^2",
    }
    text = aliases.get(text, text)
    if text not in KNOWN_UNITS:
        raise UnitError(
            f"unrecognised unit {unit!r}. Known units: "
            f"{', '.join(sorted(u for u in KNOWN_UNITS if u))}. "
            "A unit this package does not recognise is refused rather than "
            "assumed, because the failure mode of assuming is a silent factor "
            "of a thousand."
        )
    return text


def irradiance_to_photon_flux(
    irradiance_mw_mm2: float, *, wavelength_nm: float = CHR2_PEAK_WAVELENGTH_NM
) -> float:
    """Convert mW/mm^2 to photons/mm^2/s at a given wavelength.

    The opsin photocycle responds to photon arrivals, not to power, and the two
    differ by a wavelength-dependent factor. Published optogenetics protocols
    state irradiance; published photocycle models take flux. This is where that
    conversion happens, once, with the wavelength explicit.

    E_photon = h*c/lambda. At 470 nm that is 4.228e-19 J, so 1 mW/mm^2 is about
    2.37e15 photons/mm^2/s.

    Raises:
        UnitError: on a non-positive wavelength, which has no physical reading.
    """
    if wavelength_nm <= 0:
        raise UnitError(f"wavelength must be positive, got {wavelength_nm}")
    if irradiance_mw_mm2 < 0:
        raise UnitError(
            f"irradiance cannot be negative, got {irradiance_mw_mm2}. Darkness "
            "is 0.0; there is no negative light."
        )
    energy_per_photon_j = PLANCK_TIMES_C / (wavelength_nm * 1e-9)
    power_w_mm2 = irradiance_mw_mm2 * 1e-3
    return power_w_mm2 / energy_per_photon_j


def photon_flux_to_irradiance(
    flux_photons_mm2_s: float, *, wavelength_nm: float = CHR2_PEAK_WAVELENGTH_NM
) -> float:
    """Inverse of :func:`irradiance_to_photon_flux`, for reporting."""
    if wavelength_nm <= 0:
        raise UnitError(f"wavelength must be positive, got {wavelength_nm}")
    if flux_photons_mm2_s < 0:
        raise UnitError("photon flux cannot be negative")
    energy_per_photon_j = PLANCK_TIMES_C / (wavelength_nm * 1e-9)
    return flux_photons_mm2_s * energy_per_photon_j * 1e3


def opsin_current_density(
    current_pa: float, membrane_area_cm2: float
) -> float | None:
    """Convert an absolute opsin current (pA) to HH current density (uA/cm^2).

    The join between the two conventions this package has to live with. Opsin
    conductances are published whole-cell in nS; HH works per unit area. Doing
    this implicitly is how a photocurrent ends up a thousand times too large.

    Returns:
        Current density in uA/cm^2, or ``None`` if the area is not usable --
        which is an explicit absence, not a zero current.
    """
    if membrane_area_cm2 <= 0 or not math.isfinite(membrane_area_cm2):
        return None
    # pA -> uA is 1e-6; then per cm^2.
    return (current_pa * 1e-6) / membrane_area_cm2


def nS_times_mV_is_pA() -> bool:
    """Documentation-as-code for the one identity the opsin current relies on.

    1 nS * 1 mV = 1e-9 S * 1e-3 V = 1e-12 A = 1 pA. Exact, no factor. Stated as
    a function so a test can assert it and a reader can stop wondering.
    """
    return True
