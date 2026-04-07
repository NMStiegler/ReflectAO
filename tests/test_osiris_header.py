# Acronyms: DM — deformable mirror; FITS — Flexible Image Transport System; FWHM — full width at half maximum;
# LGS — laser guide star; OSIRIS — OH-Suppressing Infra-Red Imaging Spectrograph; UTC — Coordinated Universal Time.

"""Integration and unit tests for :class:`~reflectao.observation.OsirisHeaderAdapter` and OSIRIS ``read_fits_header``."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import pytest
from astropy.io import fits

from reflectao.observation import OsirisHeaderAdapter, read_fits_header

OSIRIS_FIXTURE_PATH = Path("/g3/data/kapa/2026feb26/raw/i260226_a010002.fits")


def test_osiris_adapter_synthetic_header() -> None:
    h = fits.Header()
    h["DATE-OBS"] = "2026-02-26"
    h["UTC"] = "07:58:04.70"
    h["ELAPTIME"] = 8.85
    h["AIRMASS"] = 1.06095471
    h["FILTER"] = "Kp"
    h["WAVECNTR"] = 2114
    h["O1FPS"] = "500"
    h["STINTTIM"] = "2"
    h["TARGWAVE"] = 1.65e-06
    h["MJD-OBS"] = 61097.332
    h["TELESCOP"] = "Keck I"
    h["EL"] = 70.0
    h["AZ"] = -12.0
    h["PARANG"] = 165.0
    h["TARGRA"] = 114.5
    h["TARGDEC"] = 38.88
    h["TUBETEMP"] = -2.0
    h["AOLBFWHM"] = "1.0"
    h["LGRMSWF"] = 380.0
    h["AOFCSALT"] = 90000.0
    h["DMMRFN"] = "test.mr"
    h["GUIDFWHM"] = 2.5

    out = OsirisHeaderAdapter().extract(h, source_path=Path("/tmp/x.fits"))
    assert out.exposure_start == datetime(2026, 2, 26, 7, 58, 4, 700000, tzinfo=timezone.utc)
    assert u.isclose(out.exposure_time, 8.85 * u.s)
    assert u.isclose(out.airmass, 1.06095471 * u.dimensionless_unscaled)
    assert out.filter_name == "Kp"
    assert u.isclose(out.wavelength, 2114.0 * u.nm)
    assert u.isclose(out.ao_wfs_frame_rate, 500.0 * u.Hz)
    assert u.isclose(out.ao_strap_integration_time, 0.002 * u.s)
    assert out.extra["TARGWAVE"] == 1.65e-06
    assert out.mjd_obs == pytest.approx(61097.332)
    assert out.telescope_name == "Keck I"
    assert u.isclose(out.telescope_elevation, 70.0 * u.deg)
    assert u.isclose(out.telescope_azimuth, -12.0 * u.deg)
    assert u.isclose(out.parallactic_angle, 165.0 * u.deg)
    assert u.isclose(out.target_ra, 114.5 * u.deg)
    assert u.isclose(out.target_dec, 38.88 * u.deg)
    assert u.isclose(out.tube_temperature, -2.0 * u.deg_C)
    assert u.isclose(out.ao_lbwfs_fwhm, 1.0 * u.arcsec)
    assert u.isclose(out.lgs_rms_wf_residual, 380.0 * u.nm)
    assert u.isclose(out.sodium_layer_altitude, 90000.0 * u.m)
    assert out.dm_reconstructor_file == "test.mr"
    assert out.extra["GUIDFWHM"] == 2.5


@pytest.mark.skipif(not OSIRIS_FIXTURE_PATH.is_file(), reason="OSIRIS fixture file not on this system")
def test_osiris_real_file_integration() -> None:
    out = read_fits_header(OSIRIS_FIXTURE_PATH, "osiris")
    assert out.source_path == OSIRIS_FIXTURE_PATH
    assert out.exposure_start == datetime(2026, 2, 26, 7, 58, 4, 700000, tzinfo=timezone.utc)
    assert u.isclose(out.exposure_time, 8.85 * u.s)
    assert out.filter_name == "Kp"
    assert u.isclose(out.wavelength, 2114.0 * u.nm)
    assert u.isclose(out.airmass, 1.06095471 * u.dimensionless_unscaled)
    assert u.isclose(out.ao_wfs_frame_rate, 500.0 * u.Hz)
    assert u.isclose(out.ao_strap_integration_time, 0.002 * u.s)
    assert "TARGWAVE" in out.extra
    assert out.mjd_obs == pytest.approx(61097.33199892)
    assert out.telescope_name == "Keck I"
    assert u.isclose(out.telescope_elevation, 70.46442749 * u.deg)
    assert u.isclose(out.telescope_azimuth, -12.33360863 * u.deg)
    assert u.isclose(out.parallactic_angle, 164.87581857 * u.deg)
    assert u.isclose(out.target_ra, 114.53545833 * u.deg)
    assert u.isclose(out.target_dec, 38.88191667 * u.deg)
    assert u.isclose(out.tube_temperature, -2.03654206 * u.deg_C)
    assert u.isclose(out.ao_lbwfs_fwhm, 1.0 * u.arcsec)
    assert u.isclose(out.lgs_rms_wf_residual, 379.9 * u.nm)
    assert u.isclose(out.sodium_layer_altitude, 85868.0 * u.m)
    assert out.dm_reconstructor_file == "26Feb0031.mr"
    assert "GUIDFWHM" in out.extra
    assert "WXPRESS" in out.extra
