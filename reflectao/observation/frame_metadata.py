"""Logical observation metadata extracted from FITS headers (no pixel data).

This module defines :class:`FrameMetadata`, the canonical in-memory representation
of per-frame quantities used to drive AO simulations and comparisons. Values are
normalized to documented physical units at read time; see ``docs/frame_metadata.md``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

from astropy.units import Quantity


@dataclass(frozen=True)
class FrameMetadata:
    """Header-derived quantities for a single FITS exposure (simulation-oriented).

    Physical fields use :class:`astropy.units.Quantity` with **canonical** storage
    units (nanometres, seconds, degrees, hertz, etc.). Instrument-specific
    :class:`~reflectao.observation.fits_header.HeaderAdapter` implementations
    interpret on-disk FITS values and units. Pixel arrays are intentionally omitted.

    Attributes
    ----------
    source_path : pathlib.Path or None
        Path of the file the header was read from, if applicable.
    exposure_start : datetime.datetime or None
        UTC start time when parseable from ``DATE-OBS`` / time keywords.
    exposure_time : astropy.units.Quantity or None
        On-sky exposure duration in seconds.
    wavelength : astropy.units.Quantity or None
        Effective or band-centre wavelength in nanometres.
    filter_name : str or None
        Filter name string from the header.
    airmass : astropy.units.Quantity or None
        Dimensionless airmass.
    ao_wfs_frame_rate : astropy.units.Quantity or None
        High-order WFS sampling rate in hertz.
    ao_strap_integration_time : astropy.units.Quantity or None
        STRAP integration time in seconds.
    mjd_obs : float or None
        Modified Julian Date (dimensionless) when present.
    telescope_name : str or None
        Telescope identifier string.
    telescope_elevation : astropy.units.Quantity or None
        Telescope elevation in degrees.
    telescope_azimuth : astropy.units.Quantity or None
        Telescope azimuth in degrees.
    parallactic_angle : astropy.units.Quantity or None
        Parallactic angle in degrees.
    target_ra : astropy.units.Quantity or None
        Target right ascension in degrees (header convention).
    target_dec : astropy.units.Quantity or None
        Target declination in degrees (header convention).
    tube_temperature : astropy.units.Quantity or None
        Tube temperature in kelvin (converted from Celsius on disk when applicable).
    ao_lbwfs_fwhm : astropy.units.Quantity or None
        LBWFS spot FWHM in arcseconds.
    lgs_rms_wf_residual : astropy.units.Quantity or None
        LGS-related RMS metric in nanometres (exact telemetry meaning is pipeline-specific).
    sodium_layer_altitude : astropy.units.Quantity or None
        Sodium / LGS layer altitude context in metres.
    dm_reconstructor_file : str or None
        Deformable-mirror reconstructor or control filename when logged.
    extra : collections.abc.Mapping
        Raw header values for keys listed in
        :attr:`HeaderKeywordMap.extra_keys
        <reflectao.observation.fits_header.HeaderKeywordMap.extra_keys>`
        that are not mapped to typed fields above.
    """

    source_path: Path | None = None
    exposure_start: datetime | None = None
    """UTC start (or best available reference from headers)."""

    exposure_time: Quantity | None = None
    """On-sky exposure duration (canonical unit: second)."""

    wavelength: Quantity | None = None
    """Effective or band-centre wavelength (canonical unit: nanometre)."""

    filter_name: str | None = None
    airmass: Quantity | None = None
    """Dimensionless airmass."""

    ao_wfs_frame_rate: Quantity | None = None
    """High-order WFS sampling rate (canonical unit: hertz)."""

    ao_strap_integration_time: Quantity | None = None
    """STRAP integration time (canonical unit: second; FITS often uses ms)."""

    mjd_obs: float | None = None
    """Modified Julian date when provided (dimensionless epoch)."""

    telescope_name: str | None = None
    telescope_elevation: Quantity | None = None
    """Telescope elevation (canonical unit: degree)."""

    telescope_azimuth: Quantity | None = None
    """Telescope azimuth (canonical unit: degree)."""

    parallactic_angle: Quantity | None = None
    """Parallactic angle (canonical unit: degree)."""

    target_ra: Quantity | None = None
    """Target right ascension (canonical unit: degree)."""

    target_dec: Quantity | None = None
    """Target declination (canonical unit: degree)."""

    tube_temperature: Quantity | None = None
    """Tube / relevant telescope temperature (canonical unit: kelvin)."""

    ao_lbwfs_fwhm: Quantity | None = None
    """LBWFS spot FWHM (canonical unit: arcsecond)."""

    lgs_rms_wf_residual: Quantity | None = None
    """LGS RMS WFE-style metric from telemetry (canonical unit: nanometre; verify for your pipeline)."""

    sodium_layer_altitude: Quantity | None = None
    """Sodium / LGS altitude context (canonical unit: metre)."""

    dm_reconstructor_file: str | None = None
    """DM reconstructor or control filename when logged."""

    extra: Mapping[str, Any] = field(default_factory=dict)
    """Raw header values not mapped to typed quantities."""
