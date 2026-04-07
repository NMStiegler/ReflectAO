# Acronyms: AO — adaptive optics; DM — deformable mirror; FITS — Flexible Image Transport System;
# FWHM — full width at half maximum; HDU — Header/Data Unit;
# LBWFS — low-bandwidth wavefront sensor; LGRMSWF — LGS RMS wavefront residual (example FITS keyword);
# LGS — laser guide star; MJD — Modified Julian Date; ms — millisecond; nm — nanometre;
# RMS — root mean square; STRAP — Keck low-order AO sensor; STINTTIM — STRAP integration time (example FITS keyword);
# UTC — Coordinated Universal Time; WFE — wavefront error; WFS — wavefront sensor.

"""FITS header parsing: keyword maps, adapters, and :func:`read_fits_header`.

This module is the main integration point between :class:`astropy.io.fits.Header`
and :class:`~reflectao.observation.frame_metadata.FrameMetadata`.

Developers extending instrument support subclass :class:`DefaultHeaderAdapter` and
construct a :class:`HeaderKeywordMap`; end users reading files should call
:func:`read_fits_header` with a supported ``instrument`` name.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Literal, Protocol, TypeAlias, runtime_checkable

import astropy.units as u
from astropy.io import fits
from astropy.io.fits.header import Header
from astropy.time import Time
from astropy.units import Quantity

from reflectao.observation.frame_metadata import FrameMetadata

# -----------------------------------------------------------------------------
# Type aliases
# -----------------------------------------------------------------------------

InstrumentName: TypeAlias = Literal["osiris", "nirc2"]
"""Supported instrument identifiers for :func:`read_fits_header` (lowercase)."""


# -----------------------------------------------------------------------------
# Unit and quantity helpers
# -----------------------------------------------------------------------------


def _parse_input_unit(spec: str) -> u.Unit:
    """Map a short unit string from adapter configuration to an :class:`astropy.units.Unit`.

    Parameters
    ----------
    spec : str
        Human-readable unit label (e.g. ``\"nm\"``, ``\"micron\"``, ``\"ms\"``).

    Returns
    -------
    astropy.units.Unit
        The corresponding Astropy unit.

    Raises
    ------
    ValueError
        If ``spec`` is not a recognized alias.
    """

    s = spec.lower().strip()
    if s in ("nm", "nanometer", "nanometers"):
        return u.nm
    if s in ("a", "ang", "angstrom", "angstroms", "\u212b"):
        return u.Angstrom
    if s in ("um", "micron", "microns", "micrometer", "micrometers"):
        return u.um
    if s in ("s", "sec", "second", "seconds"):
        return u.s
    if s in ("ms", "millisecond", "milliseconds"):
        return u.ms
    if s in ("hz", "hertz"):
        return u.Hz
    if s in ("deg", "degree", "degrees"):
        return u.deg
    if s in ("arcsec", "arcsecond", "arcseconds"):
        return u.arcsec
    if s in ("m", "meter", "meters", "metre", "metres"):
        return u.m
    if s in ("k", "kelvin"):
        return u.K
    raise ValueError(f"Unsupported unit string for FITS adapter: {spec!r}")


def _quantity_from_card(
    header: Header,
    key: str,
    *,
    value_unit: u.Unit,
    to_unit: u.Unit,
) -> Quantity | None:
    """Read a numeric FITS card and return a :class:`~astropy.units.Quantity`.

    Parameters
    ----------
    header : astropy.io.fits.Header
        Header to read from.
    key : str
        FITS keyword name.
    value_unit : astropy.units.Unit
        Physical unit of the raw numeric value on the card.
    to_unit : astropy.units.Unit
        Desired output unit after conversion.

    Returns
    -------
    astropy.units.Quantity or None
        Converted quantity, or ``None`` if the key is missing or not numeric.
    """

    if key not in header:
        return None
    try:
        raw = float(header[key])
    except (TypeError, ValueError):
        return None
    return (raw * value_unit).to(to_unit)


def _quantity_optional(
    header: Header,
    key: str | None,
    *,
    value_unit: u.Unit,
    to_unit: u.Unit,
) -> Quantity | None:
    """Like :func:`_quantity_from_card` but accepts a missing keyword name (returns ``None``).

    Parameters
    ----------
    header : astropy.io.fits.Header
        Header to read from.
    key : str or None
        FITS keyword, or ``None`` to skip.
    value_unit, to_unit : astropy.units.Unit
        Same as :func:`_quantity_from_card`.

    Returns
    -------
    astropy.units.Quantity or None
    """

    if key is None:
        return None
    return _quantity_from_card(header, key, value_unit=value_unit, to_unit=to_unit)


# -----------------------------------------------------------------------------
# Header card primitives
# -----------------------------------------------------------------------------


def _card_float(header: Header, key: str) -> float | None:
    """Parse a FITS card as ``float``, or ``None`` if missing or invalid.

    Parameters
    ----------
    header : astropy.io.fits.Header
    key : str

    Returns
    -------
    float or None
    """

    if key not in header:
        return None
    try:
        return float(header[key])
    except (TypeError, ValueError):
        return None


def _card_str(header: Header, key: str) -> str | None:
    """Return stripped string value for a FITS card, or ``None`` if empty/missing.

    Parameters
    ----------
    header : astropy.io.fits.Header
    key : str

    Returns
    -------
    str or None
    """

    if key not in header:
        return None
    raw = header[key]
    if raw is None:
        return None
    s = str(raw).strip()
    return s or None


# -----------------------------------------------------------------------------
# Keyword resolution policy and HeaderKeywordMap
# -----------------------------------------------------------------------------


class KeywordResolutionPolicy(str, Enum):
    """How to combine multiple candidate FITS keywords for one logical slot.

    New policies can be added here and handled in :func:`_resolve_float`.
    """

    FIRST_NON_NULL = "first_non_null"
    """Use the first candidate keyword whose value parses as a finite float."""


def _resolve_float(
    header: Header,
    candidate_keys: tuple[str, ...],
    policy: KeywordResolutionPolicy,
) -> float | None:
    """Resolve a float from an ordered list of FITS keywords.

    Parameters
    ----------
    header : astropy.io.fits.Header
    candidate_keys : tuple of str
        Keywords to try in order (e.g. ``(\"ELAPTIME\", \"EXPTIME\")``).
    policy : KeywordResolutionPolicy
        Resolution strategy.

    Returns
    -------
    float or None
        Parsed value, or ``None`` if no candidate yields a float.

    Raises
    ------
    NotImplementedError
        If ``policy`` is not implemented.
    """

    if policy is KeywordResolutionPolicy.FIRST_NON_NULL:
        for key in candidate_keys:
            v = _card_float(header, key)
            if v is not None:
                return v
        return None
    raise NotImplementedError(f"Unsupported keyword resolution policy: {policy!r}")


@dataclass(frozen=True, slots=True)
class HeaderKeywordMap:
    """Maps logical observation fields to FITS keywords and on-disk units.

    Exposure duration may be stored under several alternative keywords; use
    :attr:`exptime_keys` (ordered) with :attr:`keyword_policy` to select one.

    Attributes
    ----------
    All fields correspond to FITS keyword names unless noted. Optional fields
    use ``None`` to disable reading that quantity. See instrument modules for
    concrete keyword sets (e.g. :mod:`reflectao.observation.instruments.osiris`).
    """

    date_obs: str = "DATE-OBS"
    time_obs: str = "TIME-OBS"
    exptime_keys: tuple[str, ...] = ("EXPTIME", "EXPOSURE")
    """Ordered FITS keywords for on-sky exposure time in seconds (see :attr:`keyword_policy`)."""

    keyword_policy: KeywordResolutionPolicy = KeywordResolutionPolicy.FIRST_NON_NULL
    """Policy for resolving :attr:`exptime_keys` and any future multi-candidate fields."""

    airmass: str = "AIRMASS"
    filter_name: str = "FILTER"
    wavelength_nm: str | None = "WAVELEN"
    wavelength_input_unit: str = "nm"
    """Unit of the wavelength keyword on the FITS card (converted to nm internally)."""

    ao_wfs_frame_rate_hz: str | None = None
    ao_strap_integration_time_ms: str | None = None
    strap_integration_input_unit: str = "ms"
    """Unit of ``STINTTIM``-like keywords (converted to seconds internally)."""

    mjd_obs: str | None = None
    telescope_name: str | None = None
    telescope_el_deg: str | None = None
    telescope_az_deg: str | None = None
    parallactic_deg: str | None = None
    target_ra_deg: str | None = None
    target_dec_deg: str | None = None
    tube_temp_c: str | None = None
    """Tube temperature on disk in degrees Celsius (same unit in :attr:`FrameMetadata.tube_temperature`)."""

    ao_lbwfs_fwhm_arcsec: str | None = None
    lgs_rms_wf_residual: str | None = None
    """On-disk unit assumed nanometres for ``LGRMSWF``-style metrics."""

    lgs_rms_input_unit: str = "nm"
    sodium_layer_altitude_m: str | None = None
    dm_reconstructor_file: str | None = None

    extra_keys: tuple[str, ...] = ()
    """Additional header keys copied verbatim into :attr:`FrameMetadata.extra`."""


# -----------------------------------------------------------------------------
# HeaderAdapter protocol and read_fits_header
# -----------------------------------------------------------------------------


@runtime_checkable
class HeaderAdapter(Protocol):
    """Protocol for objects that build :class:`~reflectao.observation.frame_metadata.FrameMetadata` from a header."""

    def extract(self, header: Header, *, source_path: Path | None = None) -> FrameMetadata:
        """Parse ``header`` into a :class:`FrameMetadata` instance.

        Parameters
        ----------
        header : astropy.io.fits.Header
            Header to parse (any HDU).
        source_path : pathlib.Path, optional
            Original file path for provenance on the result.

        Returns
        -------
        FrameMetadata
        """
        ...


_adapters_by_instrument: dict[str, HeaderAdapter] | None = None


def _instrument_adapters() -> dict[str, HeaderAdapter]:
    """Lazily build the mapping from instrument name to adapter (avoids import cycles)."""

    global _adapters_by_instrument
    if _adapters_by_instrument is None:
        from reflectao.observation.instruments import Nirc2HeaderAdapter, OsirisHeaderAdapter

        _adapters_by_instrument = {
            "osiris": OsirisHeaderAdapter(),
            "nirc2": Nirc2HeaderAdapter(),
        }
    return _adapters_by_instrument


def read_fits_header(
    path: str | Path,
    instrument: InstrumentName,
    *,
    hdu: int | tuple[int, int] = 0,
) -> FrameMetadata:
    """Load a FITS header from disk and parse it with a built-in instrument adapter.

    Only the header is read; image pixels are not loaded.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the FITS file.
    instrument : {\"osiris\", \"nirc2\"}
        Which instrument keyword map and adapter to use (lowercase).
    hdu : int or tuple of int, optional
        Extension index or ``(group, number)`` passed to :func:`astropy.io.fits.getheader`.

    Returns
    -------
    FrameMetadata
        Parsed metadata.

    Raises
    ------
    ValueError
        If ``instrument`` is not one of the supported names.

    Notes
    -----
    For in-memory headers or custom keyword maps, use a concrete
    :class:`HeaderAdapter` implementation (e.g. :class:`OsirisHeaderAdapter`) and
    call :meth:`~HeaderAdapter.extract` directly. :class:`DefaultHeaderAdapter` is
    intended for developer use inside this package and tests.
    """

    adapters = _instrument_adapters()
    if instrument not in adapters:
        supported = ", ".join(sorted(adapters))
        raise ValueError(f"Unknown instrument {instrument!r}; expected one of: {supported}")

    path = Path(path)
    header = fits.getheader(path, ext=hdu)
    adapter = adapters[instrument]
    return adapter.extract(header, source_path=path)


# -----------------------------------------------------------------------------
# Field extractors (HeaderKeywordMap → scalar helpers)
# -----------------------------------------------------------------------------


def _parse_exposure_start(header: Header, keys: HeaderKeywordMap) -> datetime | None:
    """Combine date and time keywords into an aware UTC :class:`datetime.datetime`."""

    date_raw = _card_str(header, keys.date_obs)
    if date_raw is None:
        return None
    time_raw = _card_str(header, keys.time_obs)
    if time_raw is not None and "T" not in date_raw:
        iso_guess = f"{date_raw}T{time_raw}"
    else:
        iso_guess = date_raw
    try:
        t = Time(iso_guess, format="isot", scale="utc")
    except Exception:
        try:
            if time_raw is None:
                t = Time(date_raw, scale="utc")
            else:
                t = Time(f"{date_raw} {time_raw}", scale="utc")
        except Exception:
            return None
    return t.to_datetime(timezone=timezone.utc)


def _exposure_time(header: Header, keys: HeaderKeywordMap) -> Quantity | None:
    """On-sky exposure duration in seconds from :attr:`HeaderKeywordMap.exptime_keys`."""

    v = _resolve_float(header, keys.exptime_keys, keys.keyword_policy)
    if v is None:
        return None
    return (v * u.s).to(u.s)


def _wavelength(header: Header, keys: HeaderKeywordMap) -> Quantity | None:
    """Band-centre wavelength in nanometres."""

    if keys.wavelength_nm is None or keys.wavelength_nm not in header:
        return None
    src = _parse_input_unit(keys.wavelength_input_unit)
    return _quantity_from_card(
        header,
        keys.wavelength_nm,
        value_unit=src,
        to_unit=u.nm,
    )


def _airmass(header: Header, keys: HeaderKeywordMap) -> Quantity | None:
    """Dimensionless airmass."""

    return _quantity_from_card(
        header,
        keys.airmass,
        value_unit=u.dimensionless_unscaled,
        to_unit=u.dimensionless_unscaled,
    )


def _ao_wfs_rate(header: Header, keys: HeaderKeywordMap) -> Quantity | None:
    """High-order WFS rate in hertz."""

    if keys.ao_wfs_frame_rate_hz is None:
        return None
    return _quantity_from_card(
        header,
        keys.ao_wfs_frame_rate_hz,
        value_unit=u.Hz,
        to_unit=u.Hz,
    )


def _ao_strap_time(header: Header, keys: HeaderKeywordMap) -> Quantity | None:
    """STRAP integration time in seconds."""

    if keys.ao_strap_integration_time_ms is None:
        return None
    src = _parse_input_unit(keys.strap_integration_input_unit)
    return _quantity_from_card(
        header,
        keys.ao_strap_integration_time_ms,
        value_unit=src,
        to_unit=u.s,
    )


def _optional_float(header: Header, key: str | None) -> float | None:
    """Float from a single optional keyword name."""

    if key is None:
        return None
    return _card_float(header, key)


def _optional_str(header: Header, key: str | None) -> str | None:
    """String from a single optional keyword name."""

    if key is None:
        return None
    return _card_str(header, key)


def _tube_temperature_c(header: Header, keys: HeaderKeywordMap) -> Quantity | None:
    """Tube temperature in degrees Celsius from the FITS card."""

    if keys.tube_temp_c is None:
        return None
    v = _card_float(header, keys.tube_temp_c)
    if v is None:
        return None
    return Quantity(v, u.deg_C)


def _lgs_rms(header: Header, keys: HeaderKeywordMap) -> Quantity | None:
    """LGS RMS-style metric in nanometres."""

    if keys.lgs_rms_wf_residual is None:
        return None
    src = _parse_input_unit(keys.lgs_rms_input_unit)
    return _quantity_from_card(
        header,
        keys.lgs_rms_wf_residual,
        value_unit=src,
        to_unit=u.nm,
    )


def _gather_extra(header: Header, keys: HeaderKeywordMap) -> dict[str, Any]:
    """Copy selected raw cards into the ``extra`` mapping."""

    out: dict[str, Any] = {}
    for name in keys.extra_keys:
        if name in header:
            out[name] = header[name]
    return out


# -----------------------------------------------------------------------------
# DefaultHeaderAdapter (developer / internal baseline)
# -----------------------------------------------------------------------------


class DefaultHeaderAdapter:
    """Baseline header adapter using :class:`HeaderKeywordMap` conventions.

    Subclass or instantiate with a custom :class:`HeaderKeywordMap` for new
    instruments. This class is not re-exported from the top-level ``reflectao``
    package; use :func:`read_fits_header` for supported instruments or import
    from ``reflectao.observation.fits_header`` when developing adapters.

    Parameters
    ----------
    keywords : HeaderKeywordMap, optional
        Keyword map; defaults to generic FITS-like names.

    Attributes
    ----------
    keywords : HeaderKeywordMap
        The map used by :meth:`extract`.
    """

    def __init__(self, keywords: HeaderKeywordMap | None = None) -> None:
        self._keywords = keywords or HeaderKeywordMap()

    @property
    def keywords(self) -> HeaderKeywordMap:
        """Keyword configuration used for extraction."""

        return self._keywords

    def extract(self, header: Header, *, source_path: Path | None = None) -> FrameMetadata:
        """Build :class:`FrameMetadata` from a FITS header using :attr:`keywords`.

        Parameters
        ----------
        header : astropy.io.fits.Header
        source_path : pathlib.Path, optional

        Returns
        -------
        FrameMetadata
        """

        keys = self._keywords
        extra = _gather_extra(header, keys)
        return FrameMetadata(
            source_path=source_path,
            exposure_start=_parse_exposure_start(header, keys),
            exposure_time=_exposure_time(header, keys),
            wavelength=_wavelength(header, keys),
            filter_name=_card_str(header, keys.filter_name),
            airmass=_airmass(header, keys),
            ao_wfs_frame_rate=_ao_wfs_rate(header, keys),
            ao_strap_integration_time=_ao_strap_time(header, keys),
            mjd_obs=_optional_float(header, keys.mjd_obs),
            telescope_name=_optional_str(header, keys.telescope_name),
            telescope_elevation=_quantity_optional(
                header, keys.telescope_el_deg, value_unit=u.deg, to_unit=u.deg
            ),
            telescope_azimuth=_quantity_optional(
                header, keys.telescope_az_deg, value_unit=u.deg, to_unit=u.deg
            ),
            parallactic_angle=_quantity_optional(
                header, keys.parallactic_deg, value_unit=u.deg, to_unit=u.deg
            ),
            target_ra=_quantity_optional(header, keys.target_ra_deg, value_unit=u.deg, to_unit=u.deg),
            target_dec=_quantity_optional(
                header, keys.target_dec_deg, value_unit=u.deg, to_unit=u.deg
            ),
            tube_temperature=_tube_temperature_c(header, keys),
            ao_lbwfs_fwhm=_quantity_optional(
                header, keys.ao_lbwfs_fwhm_arcsec, value_unit=u.arcsec, to_unit=u.arcsec
            ),
            lgs_rms_wf_residual=_lgs_rms(header, keys),
            sodium_layer_altitude=_quantity_optional(
                header, keys.sodium_layer_altitude_m, value_unit=u.m, to_unit=u.m
            ),
            dm_reconstructor_file=_optional_str(header, keys.dm_reconstructor_file),
            extra=extra,
        )
