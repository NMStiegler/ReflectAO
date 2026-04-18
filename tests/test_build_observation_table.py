from pathlib import Path

import pytest
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
import numpy as np

from reflectao.build_observation_table import build_observation_table
from reflectao.schema import new_empty_observation_table


def test_build_observation_table_one_row_with_data_from_fits_header(tmp_path: Path):
    TEST_PATH = "/g3/data/kapa/2026feb26/raw/i260226_a010002.fits"

    tbl = build_observation_table(TEST_PATH, instrument="OSIRIS")
    assert len(tbl) == 1
    assert tbl["image_path"][0] == TEST_PATH
    assert tbl["telescope_name"][0] == "Keck I"
    assert tbl["telescope_elevation"][0] == 70.46442749 * u.deg
    assert tbl["telescope_azimuth"][0] == -12.33360863  * u.deg
    assert tbl["plate_scale"][0] == 0.0099576 * u.arcsec / u.pixel
    assert tbl["frame_number"][0] == "2"
    assert tbl["set_number"][0] == "10"
    assert tbl["set_name"][0] == "ngc2419"
    assert tbl["target_name"][0] == "NGC_2419"
    assert tbl["object_name"][0] == "ngc2419"
    assert tbl["target_ra"][0] == 114.53545833  * u.deg
    assert tbl["target_dec"][0] == 38.88191667  * u.deg
    assert tbl["epoch"][0] == '2000'
    assert tbl["t_exposure_start"][0] == Time("2026-02-26 07:58:04.700")
    assert tbl["t_exposure_duration"][0] == 8.85 * u.s
    assert tbl["t_int"][0] == 4.425 * u.s
    assert tbl["num_coadds"][0] == 2
    assert tbl["wavelength"][0] == 2120 * u.nm
    assert tbl["filter_name"][0] == "Kp"
    assert tbl["airmass"][0] == 1.06095471
    assert tbl["zenith_angle"][0] == (np.arccos(1.0 / float(tbl["airmass"][0])) * u.rad).to(u.deg) 
    assert tbl["aborted"][0] == False
    assert tbl["dither_name"][0] == "Box4"
    assert tbl["lgs_wfs_rate"][0] == 500 * u.Hz
    assert tbl["tt_wfs_rate"][0] == 500 * u.Hz
    assert tbl["tt_wfs_centroid_gain"][0] == 0.711
    assert tbl["lgs_rms_wfe"][0] == 379.9 * u.nm
    assert tbl["lgs_layer_alt"][0] == 85868 * u.m
    assert tbl["OSIRIS_tt_sensor"][0] == "STRAP"
    assert tbl["ngs_fwhm"][0] == 2.834221576850309 * u.arcsec
    assert tbl["ngs_wavelength"][0] == (6.5E-07 * u.m).to(u.m)
    assert tbl["reconstructor_name"][0] == "26Feb0031.mr"
    assert tbl["dm_gain"][0] == 0.4
    assert tbl["lgs_wfs_detector_gain"][0] == 600
    assert tbl["system_gain"][0] is np.ma.masked # Not sure if it's listed in this file. The header for it is the detector gain
    assert tbl["utt_gain"][0] == 0.1
    assert tbl["dtt_gain"][0] == 0.1
    assert tbl["ao_mode"][0] == "3" # <-- What is this?
    assert tbl["humidity_dome"][0] == 2.7
    assert tbl["T_dome_air"][0] == -0.5 * u.deg_C
    assert tbl["humidity_outside"][0] == 4.0
    assert tbl["T_outside_air"][0] == -0.7 * u.deg_C
    assert tbl["P_barometric"][0] == (614.72 * u.hPa).to(u.Pa)
    assert tbl["wind_direction"][0] == 212.5 * u.deg
    assert tbl["wind_speed"][0] == 0.41 * u.m / u.s
    assert tbl["t_weather_sample"][0] == "07:58:05.34"
    assert tbl["T_tube"][0] == -2.03654206 * u.deg_C

def test_build_observation_table_two_files(tmp_path: Path):
    TEST_PATH = ["/g3/data/kapa/2026feb26/raw/i260226_a010002.fits", "/g3/data/kapa/2026feb26/raw/i260226_a010002.fits"]

    tbl = build_observation_table(TEST_PATH, instrument="OSIRIS")
    assert len(tbl) == 2
    
    index = 0
    assert tbl["image_path"][index] == TEST_PATH[index]
    assert tbl["telescope_name"][index] == "Keck I"
    assert tbl["telescope_elevation"][index] == 70.46442749 * u.deg
    assert tbl["telescope_azimuth"][index] == -12.33360863  * u.deg
    assert tbl["plate_scale"][index] == 0.0099576 * u.arcsec / u.pixel
    assert tbl["frame_number"][index] == "2"
    assert tbl["set_number"][index] == "10"
    assert tbl["set_name"][index] == "ngc2419"
    assert tbl["target_name"][index] == "NGC_2419"
    assert tbl["object_name"][index] == "ngc2419"
    assert tbl["target_ra"][index] == 114.53545833  * u.deg
    assert tbl["target_dec"][index] == 38.88191667  * u.deg
    assert tbl["epoch"][index] == '2000'
    assert tbl["t_exposure_start"][index] == Time("2026-02-26 07:58:04.700")
    assert tbl["t_exposure_duration"][index] == 8.85 * u.s
    assert tbl["t_int"][index] == 4.425 * u.s
    assert tbl["num_coadds"][index] == 2
    assert tbl["wavelength"][index] == 2120 * u.nm
    assert tbl["filter_name"][index] == "Kp"
    assert tbl["airmass"][index] == 1.06095471
    assert tbl["zenith_angle"][index] == (np.arccos(1.0 / float(tbl["airmass"][index])) * u.rad).to(u.deg)
    assert tbl["aborted"][index] == False
    assert tbl["dither_name"][index] == "Box4"
    assert tbl["lgs_wfs_rate"][index] == 500 * u.Hz
    assert tbl["tt_wfs_rate"][index] == 500 * u.Hz
    assert tbl["tt_wfs_centroid_gain"][index] == 0.711
    assert tbl["lgs_rms_wfe"][index] == 379.9 * u.nm
    assert tbl["lgs_layer_alt"][index] == 85868 * u.m
    assert tbl["OSIRIS_tt_sensor"][index] == "STRAP"
    assert tbl["ngs_fwhm"][index] == 2.834221576850309 * u.arcsec
    assert tbl["ngs_wavelength"][index] == (6.5E-07 * u.m).to(u.m)
    assert tbl["reconstructor_name"][index] == "26Feb0031.mr"
    assert tbl["dm_gain"][index] == 0.4
    assert tbl["lgs_wfs_detector_gain"][index] == 600 # <-- Check this it's in counts
    assert tbl["system_gain"][index] is np.ma.masked
    assert tbl["utt_gain"][index] == 0.1
    assert tbl["dtt_gain"][index] == 0.1
    assert tbl["ao_mode"][index] == "3" # <-- What is this?
    assert tbl["humidity_dome"][index] == 2.7
    assert tbl["T_dome_air"][index] == -0.5 * u.deg_C
    assert tbl["humidity_outside"][index] == 4.0
    assert tbl["T_outside_air"][index] == -0.7 * u.deg_C
    assert tbl["P_barometric"][index] == (614.72 * u.hPa).to(u.Pa)
    assert tbl["wind_direction"][index] == 212.5 * u.deg
    assert tbl["wind_speed"][index] == 0.41 * u.m / u.s
    assert tbl["t_weather_sample"][index] == "07:58:05.34"
    assert tbl["T_tube"][index] == -2.03654206 * u.deg_C

    index = 1
    assert tbl["image_path"][index] == TEST_PATH[index]
    assert tbl["telescope_name"][index] == "Keck I"
    assert tbl["telescope_elevation"][index] == 70.46442749 * u.deg
    assert tbl["telescope_azimuth"][index] == -12.33360863  * u.deg
    assert tbl["plate_scale"][index] == 0.0099576 * u.arcsec / u.pixel
    assert tbl["frame_number"][index] == "2"
    assert tbl["set_number"][index] == "10"
    assert tbl["set_name"][index] == "ngc2419"
    assert tbl["target_name"][index] == "NGC_2419"
    assert tbl["object_name"][index] == "ngc2419"
    assert tbl["target_ra"][index] == 114.53545833  * u.deg
    assert tbl["target_dec"][index] == 38.88191667  * u.deg
    assert tbl["epoch"][index] == '2000'
    assert tbl["t_exposure_start"][index] == Time("2026-02-26 07:58:04.700")
    assert tbl["t_exposure_duration"][index] == 8.85 * u.s
    assert tbl["t_int"][index] == 4.425 * u.s
    assert tbl["num_coadds"][index] == 2
    assert tbl["wavelength"][index] == 2120 * u.nm
    assert tbl["filter_name"][index] == "Kp"
    assert tbl["airmass"][index] == 1.06095471
    assert tbl["zenith_angle"][index] == (np.arccos(1.0 / float(tbl["airmass"][index])) * u.rad).to(u.deg)
    assert tbl["aborted"][index] == False
    assert tbl["dither_name"][index] == "Box4"
    assert tbl["lgs_wfs_rate"][index] == 500 * u.Hz
    assert tbl["tt_wfs_rate"][index] == 500 * u.Hz
    assert tbl["tt_wfs_centroid_gain"][index] == 0.711
    assert tbl["lgs_rms_wfe"][index] == 379.9 * u.nm
    assert tbl["lgs_layer_alt"][index] == 85868 * u.m
    assert tbl["OSIRIS_tt_sensor"][index] == "STRAP"
    assert tbl["ngs_fwhm"][index] == 2.834221576850309 * u.arcsec
    assert tbl["ngs_wavelength"][index] == (6.5E-07 * u.m).to(u.m)
    assert tbl["reconstructor_name"][index] == "26Feb0031.mr"
    assert tbl["dm_gain"][index] == 0.4
    assert tbl["lgs_wfs_detector_gain"][index] == 600 # <-- Check this it's in counts
    assert tbl["system_gain"][index] is np.ma.masked
    assert tbl["utt_gain"][index] == 0.1
    assert tbl["dtt_gain"][index] == 0.1
    assert tbl["ao_mode"][index] == "3" # <-- What is this?
    assert tbl["humidity_dome"][index] == 2.7
    assert tbl["T_dome_air"][index] == -0.5 * u.deg_C
    assert tbl["humidity_outside"][index] == 4.0
    assert tbl["T_outside_air"][index] == -0.7 * u.deg_C
    assert tbl["P_barometric"][index] == (614.72 * u.hPa).to(u.Pa)
    assert tbl["wind_direction"][index] == 212.5 * u.deg
    assert tbl["wind_speed"][index] == 0.41 * u.m / u.s
    assert tbl["t_weather_sample"][index] == "07:58:05.34"
    assert tbl["T_tube"][index] == -2.03654206 * u.deg_C


def test_build_observation_table_extending_table(tmp_path: Path):
    TEST_PATH = "/g3/data/kapa/2026feb26/raw/i260226_a010002.fits"

    tbl0 = build_observation_table(TEST_PATH, instrument="OSIRIS")
    assert len(tbl0) == 1
    
    tbl = build_observation_table(TEST_PATH, instrument="OSIRIS", table=tbl0)
    assert len(tbl) == 2
    
    index = 0
    assert tbl["image_path"][index] == TEST_PATH
    assert tbl["telescope_name"][index] == "Keck I"
    assert tbl["telescope_elevation"][index] == 70.46442749 * u.deg
    assert tbl["telescope_azimuth"][index] == -12.33360863  * u.deg
    assert tbl["plate_scale"][index] == 0.0099576 * u.arcsec / u.pixel
    assert tbl["frame_number"][index] == "2"
    assert tbl["set_number"][index] == "10"
    assert tbl["set_name"][index] == "ngc2419"
    assert tbl["target_name"][index] == "NGC_2419"
    assert tbl["object_name"][index] == "ngc2419"
    assert tbl["target_ra"][index] == 114.53545833  * u.deg
    assert tbl["target_dec"][index] == 38.88191667  * u.deg
    assert tbl["epoch"][index] == '2000'
    assert tbl["t_exposure_start"][index] == Time("2026-02-26 07:58:04.700")
    assert tbl["t_exposure_duration"][index] == 8.85 * u.s
    assert tbl["t_int"][index] == 4.425 * u.s
    assert tbl["num_coadds"][index] == 2
    assert tbl["wavelength"][index] == 2120 * u.nm
    assert tbl["filter_name"][index] == "Kp"
    assert tbl["airmass"][index] == 1.06095471
    assert tbl["zenith_angle"][index] == (np.arccos(1.0 / float(tbl["airmass"][index])) * u.rad).to(u.deg)
    assert tbl["aborted"][index] == False
    assert tbl["dither_name"][index] == "Box4"
    assert tbl["lgs_wfs_rate"][index] == 500 * u.Hz
    assert tbl["tt_wfs_rate"][index] == 500 * u.Hz
    assert tbl["tt_wfs_centroid_gain"][index] == 0.711
    assert tbl["lgs_rms_wfe"][index] == 379.9 * u.nm
    assert tbl["lgs_layer_alt"][index] == 85868 * u.m
    assert tbl["OSIRIS_tt_sensor"][index] == "STRAP"
    assert tbl["ngs_fwhm"][index] == 2.834221576850309 * u.arcsec
    assert tbl["ngs_wavelength"][index] == (6.5E-07 * u.m).to(u.m)
    assert tbl["reconstructor_name"][index] == "26Feb0031.mr"
    assert tbl["dm_gain"][index] == 0.4
    assert tbl["lgs_wfs_detector_gain"][index] == 600 # <-- Check this it's in counts
    assert tbl["system_gain"][index] is np.ma.masked
    assert tbl["utt_gain"][index] == 0.1
    assert tbl["dtt_gain"][index] == 0.1
    assert tbl["ao_mode"][index] == "3" # <-- What is this?
    assert tbl["humidity_dome"][index] == 2.7
    assert tbl["T_dome_air"][index] == -0.5 * u.deg_C
    assert tbl["humidity_outside"][index] == 4.0
    assert tbl["T_outside_air"][index] == -0.7 * u.deg_C
    assert tbl["P_barometric"][index] == (614.72 * u.hPa).to(u.Pa)
    assert tbl["wind_direction"][index] == 212.5 * u.deg
    assert tbl["wind_speed"][index] == 0.41 * u.m / u.s
    assert tbl["t_weather_sample"][index] == "07:58:05.34"
    assert tbl["T_tube"][index] == -2.03654206 * u.deg_C

    index = 1
    assert tbl["image_path"][index] == TEST_PATH
    assert tbl["telescope_name"][index] == "Keck I"
    assert tbl["telescope_elevation"][index] == 70.46442749 * u.deg
    assert tbl["telescope_azimuth"][index] == -12.33360863  * u.deg
    assert tbl["plate_scale"][index] == 0.0099576 * u.arcsec / u.pixel
    assert tbl["frame_number"][index] == "2"
    assert tbl["set_number"][index] == "10"
    assert tbl["set_name"][index] == "ngc2419"
    assert tbl["target_name"][index] == "NGC_2419"
    assert tbl["object_name"][index] == "ngc2419"
    assert tbl["target_ra"][index] == 114.53545833  * u.deg
    assert tbl["target_dec"][index] == 38.88191667  * u.deg
    assert tbl["epoch"][index] == '2000'
    assert tbl["t_exposure_start"][index] == Time("2026-02-26 07:58:04.700")
    assert tbl["t_exposure_duration"][index] == 8.85 * u.s
    assert tbl["t_int"][index] == 4.425 * u.s
    assert tbl["num_coadds"][index] == 2
    assert tbl["wavelength"][index] == 2120 * u.nm
    assert tbl["filter_name"][index] == "Kp"
    assert tbl["airmass"][index] == 1.06095471
    assert tbl["zenith_angle"][index] == (np.arccos(1.0 / float(tbl["airmass"][index])) * u.rad).to(u.deg)
    assert tbl["aborted"][index] == False
    assert tbl["dither_name"][index] == "Box4"
    assert tbl["lgs_wfs_rate"][index] == 500 * u.Hz
    assert tbl["tt_wfs_rate"][index] == 500 * u.Hz
    assert tbl["tt_wfs_centroid_gain"][index] == 0.711
    assert tbl["lgs_rms_wfe"][index] == 379.9 * u.nm
    assert tbl["lgs_layer_alt"][index] == 85868 * u.m
    assert tbl["OSIRIS_tt_sensor"][index] == "STRAP"
    assert tbl["ngs_fwhm"][index] == 2.834221576850309 * u.arcsec
    assert tbl["ngs_wavelength"][index] == (6.5E-07 * u.m).to(u.m)
    assert tbl["reconstructor_name"][index] == "26Feb0031.mr"
    assert tbl["dm_gain"][index] == 0.4
    assert tbl["lgs_wfs_detector_gain"][index] == 600 # <-- Check this it's in counts
    assert tbl["system_gain"][index] is np.ma.masked
    assert tbl["utt_gain"][index] == 0.1
    assert tbl["dtt_gain"][index] == 0.1
    assert tbl["ao_mode"][index] == "3" # <-- What is this?
    assert tbl["humidity_dome"][index] == 2.7
    assert tbl["T_dome_air"][index] == -0.5 * u.deg_C
    assert tbl["humidity_outside"][index] == 4.0
    assert tbl["T_outside_air"][index] == -0.7 * u.deg_C
    assert tbl["P_barometric"][index] == (614.72 * u.hPa).to(u.Pa)
    assert tbl["wind_direction"][index] == 212.5 * u.deg
    assert tbl["wind_speed"][index] == 0.41 * u.m / u.s
    assert tbl["t_weather_sample"][index] == "07:58:05.34"
    assert tbl["T_tube"][index] == -2.03654206 * u.deg_C


