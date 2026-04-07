# Acronyms: FITS — Flexible Image Transport System; HDU — Header/Data Unit; UTC — Coordinated Universal Time.

"""Tests for :class:`~reflectao.observation.FrameMetadata` and :class:`DefaultHeaderAdapter`.

:class:`~reflectao.observation.fits_header.DefaultHeaderAdapter` is imported from
``reflectao.observation.fits_header`` (developer entry point), not the package root.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits

from reflectao.observation import FrameMetadata, read_fits_header
from reflectao.observation.fits_header import DefaultHeaderAdapter, HeaderKeywordMap


def test_extract_from_in_memory_header() -> None:
    """Custom :class:`HeaderKeywordMap` populates :class:`FrameMetadata` as expected."""

    h = fits.Header()
    h["DATE-OBS"] = "2024-01-15"
    h["TIME-OBS"] = "10:30:00"
    h["EXPTIME"] = 45.0
    h["AIRMASS"] = 1.15
    h["FILTER"] = "Kbb"
    h["WAVELEN"] = 2200.0
    h["NOTE"] = "ignored unless listed in extra_keys"

    adapter = DefaultHeaderAdapter(
        HeaderKeywordMap(extra_keys=("NOTE",)),
    )
    out = adapter.extract(h, source_path=Path("/tmp/example.fits"))

    assert out.source_path == Path("/tmp/example.fits")
    assert u.isclose(out.exposure_time, 45.0 * u.s)
    assert u.isclose(out.airmass, 1.15 * u.dimensionless_unscaled)
    assert out.filter_name == "Kbb"
    assert u.isclose(out.wavelength, 2200.0 * u.nm)
    assert out.exposure_start == datetime(2024, 1, 15, 10, 30, 0, tzinfo=timezone.utc)
    assert out.extra["NOTE"] == "ignored unless listed in extra_keys"


def test_wavelength_angstrom_conversion() -> None:
    """Wavelength input unit string is applied when reading the card."""

    h = fits.Header()
    h["DATE-OBS"] = "2024-01-15T10:30:00"
    h["EXPTIME"] = 1.0
    h["WAVELEN"] = 22000.0
    adapter = DefaultHeaderAdapter(
        HeaderKeywordMap(wavelength_input_unit="angstrom"),
    )
    out = adapter.extract(h)
    assert u.isclose(out.wavelength, 2200.0 * u.nm)


def test_exposure_fallback_keyword(tmp_path: Path) -> None:
    """First missing ``exptime_keys`` candidate falls through to the next."""

    hdu = fits.PrimaryHDU(data=np.zeros((1, 1), dtype=np.float32))
    hdu.header["DATE-OBS"] = "2024-06-01"
    hdu.header["EXPOSURE"] = 12.5
    path = tmp_path / "mini.fits"
    hdu.writeto(path)

    adapter = DefaultHeaderAdapter(HeaderKeywordMap(exptime_keys=("MISSING", "EXPOSURE")))
    header = fits.getheader(path)
    out = adapter.extract(header, source_path=path)
    assert u.isclose(out.exposure_time, 12.5 * u.s)
    assert out.source_path == path


def test_read_fits_header_unknown_instrument(tmp_path: Path) -> None:
    """Unsupported ``instrument`` raises :class:`ValueError` before opening the file."""

    path = tmp_path / "empty.fits"
    fits.PrimaryHDU(data=np.zeros((1, 1), dtype=np.float32)).writeto(path)

    with pytest.raises(ValueError, match="Unknown instrument"):
        read_fits_header(path, "bogus")  # type: ignore[arg-type]


def test_frame_metadata_frozen() -> None:
    """``FrameMetadata`` is an immutable dataclass."""

    o = FrameMetadata(exposure_time=1.0 * u.s)
    with pytest.raises(AttributeError):
        o.exposure_time = 2.0 * u.s  # type: ignore[misc]
