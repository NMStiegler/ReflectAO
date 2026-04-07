"""NIRC2 header mapping vs reference ``maos_utils.run_maos_comp_to_sky_sim`` card usage.

Covers :class:`~reflectao.observation.Nirc2HeaderAdapter`, :func:`~reflectao.observation.read_fits_header`,
and comparison helpers that mirror archived PAARTI conventions.
"""

from __future__ import annotations

import math
import os
from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import pytest
from astropy.io import fits

from reflectao.observation import FrameMetadata, Nirc2HeaderAdapter, read_fits_header


def _reference_targwave_microns_to_m(hdr: fits.Header) -> float:
    """Reference: ``float(hdr['TARGWAVE']) * 1.0e-6`` (microns → metres)."""

    return float(hdr["TARGWAVE"]) * 1.0e-6


def _reference_strap_rate_hz(hdr: fits.Header) -> float:
    """Reference: ``(1.0 / float(STINTTIM)) * 1000.0`` with STINTTIM in ms."""

    return (1.0 / float(hdr["STINTTIM"])) * 1000.0


def _reference_max_wfs_rate_hz(hdr: fits.Header) -> float:
    sh = float(hdr["WSFRRT"])
    st = _reference_strap_rate_hz(hdr)
    return max(sh, st)


def test_nirc2_synthetic_matches_reference_maos_header_usage() -> None:
    h = fits.Header()
    h["DATE-OBS"] = "2017-08-23"
    h["TIME-OBS"] = "14:30:00.00"
    h["EXPTIME"] = 30.0
    h["AIRMASS"] = 1.08
    h["FILTER"] = "Kp"
    h["TARGWAVE"] = 2.12
    h["STINTTIM"] = 2.0
    h["WSFRRT"] = 472.0
    h["AOLBFWHM"] = 1.25
    h["MJD-OBS"] = 57990.5
    h["TELESCOP"] = "Keck II"
    h["TUBETEMP"] = -1.5
    h["LGRMSWF"] = 350.0

    out = Nirc2HeaderAdapter().extract(h, source_path=Path("/tmp/nirc2.fits"))

    assert out.exposure_start == datetime(2017, 8, 23, 14, 30, 0, tzinfo=timezone.utc)
    assert u.isclose(out.exposure_time, 30.0 * u.s)
    assert u.isclose(out.airmass, 1.08 * u.dimensionless_unscaled)
    assert out.filter_name == "Kp"
    assert u.isclose(out.wavelength, 2120.0 * u.nm)
    assert u.isclose(out.ao_wfs_frame_rate, 472.0 * u.Hz)
    assert u.isclose(out.ao_strap_integration_time, 0.002 * u.s)
    assert u.isclose(out.ao_lbwfs_fwhm, 1.25 * u.arcsec)
    assert out.mjd_obs == pytest.approx(57990.5)
    assert out.telescope_name == "Keck II"
    assert u.isclose(out.tube_temperature, -1.5 * u.deg_C)
    assert u.isclose(out.lgs_rms_wf_residual, 350.0 * u.nm)

    wvl_m_reflectao = out.wavelength.to(u.m).value
    wvl_m_ref = _reference_targwave_microns_to_m(h)
    assert wvl_m_reflectao == pytest.approx(wvl_m_ref)

    z_ref = math.degrees(math.acos(1.0 / float(h["AIRMASS"])))
    z_from_out = math.degrees(math.acos(1.0 / out.airmass.value))
    assert z_from_out == pytest.approx(z_ref)

    assert _reference_max_wfs_rate_hz(h) == pytest.approx(500.0)


def test_nirc2_reference_zenith_and_dtrat_ingredients() -> None:
    h = fits.Header()
    h["DATE-OBS"] = "2017-08-23T14:30:00"
    h["AIRMASS"] = 1.25
    h["TARGWAVE"] = 2.2
    h["STINTTIM"] = 2.5
    h["WSFRRT"] = 800.0
    h["FILTER"] = "Kp"
    h["EXPTIME"] = 10.0

    out = Nirc2HeaderAdapter().extract(h)
    hdr_shwfs_frame_rate = out.ao_wfs_frame_rate.to_value(u.Hz)
    hdr_strap_int_time_ms = out.ao_strap_integration_time.to_value(u.ms)
    hdr_strap_frame_rate_hz = (1.0 / hdr_strap_int_time_ms) * 1000.0
    max_frame_rate = max(hdr_shwfs_frame_rate, hdr_strap_frame_rate_hz)
    sim_dt_s = 1.0 / max_frame_rate
    hdr_shwfs_int_time_ms = (1.0 / hdr_shwfs_frame_rate) * 1000.0

    assert sim_dt_s == pytest.approx(1.0 / 800.0)
    assert hdr_shwfs_int_time_ms == pytest.approx(1.25)
    assert hdr_strap_frame_rate_hz == pytest.approx(400.0)


# Default: PAARTI-style tree on department storage (override with REFLECTAO_NIRC2_PSF).
NIRC2_EXAMPLE_PSF = Path("/u/bdigia/work/ao/airopa_input/20170823nirc2_kp/c2061_psf.fits")


def _nirc2_psf_fixture_path() -> Path:
    env = os.environ.get("REFLECTAO_NIRC2_PSF", "").strip()
    return Path(env) if env else NIRC2_EXAMPLE_PSF


_NIRC2_REAL = _nirc2_psf_fixture_path()


@pytest.mark.skipif(
    not _NIRC2_REAL.is_file(),
    reason="NIRC2 example FITS missing: mount airopa_input or set REFLECTAO_NIRC2_PSF",
)
def test_nirc2_real_file_c2061_example() -> None:
    """End-to-end read of ``c2061_psf.fits`` (2017-08-23 night); header values checked 2026-04."""

    out = read_fits_header(_NIRC2_REAL, "nirc2")
    assert out.source_path is not None and out.source_path.samefile(_NIRC2_REAL)
    assert u.isclose(out.wavelength, 2124.0 * u.nm)
    assert out.filter_name == "Kp + clear"
    assert u.isclose(out.airmass, 1.53333 * u.dimensionless_unscaled)
    assert u.isclose(out.ao_wfs_frame_rate, 1000.0 * u.Hz)
    assert u.isclose(out.ao_strap_integration_time, 0.001 * u.s)
    assert out.mjd_obs == pytest.approx(57988.27265)
    assert out.telescope_name == "Keck II"
    assert out.telescope_elevation is not None
    assert out.telescope_azimuth is not None
    assert out.ao_lbwfs_fwhm is not None
    assert out.lgs_rms_wf_residual is not None
    assert out.tube_temperature is not None
    assert "EXPSTART" in out.extra
    assert "EXPSTOP" in out.extra


def test_frame_metadata_accepts_quantity_strap() -> None:
    o = FrameMetadata(ao_strap_integration_time=3e-3 * u.s)
    assert u.isclose(o.ao_strap_integration_time, 0.003 * u.s)
