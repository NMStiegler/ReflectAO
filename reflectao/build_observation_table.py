"""
Build the canonical per-frame observation table.

Key project invariant:
- every schema field has a column in the internal table
- missing values are represented as **masked** entries (not ``None``)

Instrument-specific header interpretation should come from KAI. If ReflectAO
needs a new piece of metadata and KAI doesn't provide it, we should add it to
KAI rather than duplicating logic here.
"""

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.table import vstack, QTable

import paarti.utils.maos_utils as mu

from .schema import new_empty_observation_table, validate_table_has_schema, schema_column_names, SCHEMA
from .telemetry_utils import *

def _normalize_paths(fits_paths):
    """
    Normalize a single path or list of paths into a list of `Path`.

    :param fits_paths: A single FITS path or a list/tuple of FITS paths.
    :type fits_paths: str or pathlib.Path or list
    :return: List of `Path` objects.
    :rtype: list[pathlib.Path]
    """
    if isinstance(fits_paths, (str, Path)):
        return [Path(fits_paths)]
    return [Path(p) for p in fits_paths]


def _get_kai_instrument(instrument_name):
    """
    Return a KAI instrument object from an instrument name string.

    :param instrument_name: Instrument name, e.g. "OSIRIS".
    :type instrument_name: str
    :raises NotImplementedError: If the instrument is not supported yet.
    :return: KAI instrument instance.
    :rtype: kai.instruments.Instrument
    """
    # Import inside the function so importing reflectao doesn't require KAI
    # until the user actually builds a table.
    from kai.instruments import OSIRIS

    if str(instrument_name).upper() == "OSIRIS":
        return OSIRIS()

    raise NotImplementedError(
        "Instrument not supported yet: {0}. For now only OSIRIS is implemented.".format(
            instrument_name
        )
    )


def build_observation_table(
    fits_paths,
    instrument,
    telemetry_paths=None,
    table=None,
    verbose=False
):
    """
    Build or extend the canonical per-frame observation table.

    :param fits_paths: Path to one FITS file or a list of FITS files.
        One path produces one row; a list produces one row per file.
    :type fits_paths: str or pathlib.Path or list
    :param instrument: Instrument identifier string (e.g. "OSIRIS").
    :type instrument: str
    :param telemetry_paths: Optional paths to already-downloaded telemetry. When not provided, the code will attempt to auto-discover the telemetry paths.
    :type telemetry_paths: str or pathlib.Path, optional
    :param table: Optional existing observation table to extend.
        If provided, it must already contain all schema columns.
    :type table: astropy.table.Table, optional
    :param verbose: If True, print additional information during processing.
    :type verbose: bool, optional
    :raises ValueError: If `table` is provided but does not match the schema.
    :raises NotImplementedError: If the instrument is not supported yet.
    :return: Observation table with one row per FITS file and all schema columns.
    :rtype: astropy.table.Table
    """

    # Use KAI to read the data we want from our FITS headers
    inst = _get_kai_instrument(instrument)

    # Make sure we end up with a list of paths regardless of what
    # was passed in
    fits_paths = _normalize_paths(fits_paths)
    telemetry_paths = _normalize_paths(telemetry_paths) if telemetry_paths is not None else []
    
    # Setup the table to add to. If no table is provided, create a new one.
    if table is not None:
        validate_table_has_schema(table)

    # Helper to clean up the try/except blocks. 
    # Fetches the value, optionally applies a unit, and returns masked on failure.
    def _safe_get(func, hdr, unit=None):
        try:
            val = func(hdr)
            return val * unit if unit is not None else val
        except KeyError:
            return np.ma.masked
    
    # Accumulate all new rows in a list of dictionaries first
    new_rows = []
    expected_cols = schema_column_names(SCHEMA) # <-- Grab the full list of schema columns

    # Add a row to the table for each image, pointed to by a path in our list
    for index, p in enumerate(fits_paths):
        telemetry_file = telemetry_paths[index] if index < len(telemetry_paths) else np.ma.masked
        if telemetry_file is None or telemetry_file is np.ma.masked:
            try:
                telemetry_file = get_telemetry_path_from_image_path(p)
            except Exception as e:
                if verbose:
                    print(f"Error occurred while fetching telemetry path for {p}: {e}")
                telemetry_file = None

        # Get the fits header for the current image
        hdr = fits.getheader(p)
        row_vals = {col: np.ma.masked for col in expected_cols}

        ############# Add metadata from this row's FITS file #############
        ### System information ###
        row_vals["image_path"] = str(p) # Path to the FITS file
        row_vals['telemetry_file_path'] = str(telemetry_file) if (telemetry_file is not np.ma.masked and telemetry_file is not None) else np.ma.masked
        
        ### Telescope / instrument / set information ###
        row_vals['telescope_name'] = inst.get_telescope_name(hdr)
        row_vals['instrument_name'] = inst.get_instrument_name(hdr)
        row_vals['OSIRIS_imaging_mode'] = inst.get_imaging_mode(hdr) # For OSIRIS, get the imaging mode (imag or spec)
        row_vals['telescope_elevation'] = inst.get_telescope_elevation(hdr) * u.deg
        row_vals['telescope_azimuth'] = inst.get_telescope_azimuth(hdr) * u.deg
        row_vals['plate_scale'] = inst.get_plate_scale(hdr) * u.arcsec / u.pixel
        row_vals['frame_number'] = inst.get_frame_number(hdr)
        row_vals['set_number'] = inst.get_set_number(hdr)
        row_vals['set_name'] = inst.get_set_name(hdr)

        ### Target information ###
        row_vals['target_name'] = inst.get_target_name(hdr)
        row_vals['object_name'] = inst.get_object_name(hdr)
        row_vals['target_ra'] = inst.get_target_ra(hdr) * u.deg
        row_vals['target_dec'] = inst.get_target_dec(hdr) * u.deg
        row_vals['epoch'] = inst.get_epoch(hdr)
        
        ### Exposure information ###
        row_vals['t_exposure_start'] = Time(f'{inst.get_exposure_start_date(hdr)} {inst.get_exposure_start_time(hdr)}', scale='utc')        
        row_vals['t_exposure_duration'] = inst.get_exposure_duration(hdr) * u.s
        row_vals['t_int'] = inst.get_integration_time_per_coadd(hdr) * u.s
        row_vals['num_coadds'] = inst.get_number_of_coadds(hdr)
        # Check if t_exposure_duration matches t_int * num_coadds
        if row_vals['t_exposure_duration'] - row_vals['t_int'] * row_vals['num_coadds'] >= 0.1 * u.s:
            if verbose: print(f"Warning: for file {p}, t_exposure_duration ({row_vals['t_exposure_duration']}) does not match t_int * num_coadds ({row_vals['t_int'] * row_vals['num_coadds']})")
        row_vals['wavelength'] = (inst.get_central_wavelength(hdr) * u.um).to(u.nm) if inst.get_central_wavelength(hdr) is not None else np.ma.masked
        row_vals['filter_name'] = inst.get_filter_name(hdr)
        row_vals['airmass'] = inst.get_airmass(hdr)
        row_vals['zenith_angle'] = (np.arccos(1.0 / float(hdr['AIRMASS'])) * u.rad).to(u.deg)
        row_vals['aborted'] = not 'f' in inst.was_exposure_aborted(hdr).strip(" ").lower()
        row_vals['dither_name'] = inst.get_dither_name(hdr)

        ### AO system information ###
        row_vals['ao_closed'] = inst.was_waiting_for_DM_lock(hdr) and inst.was_DM_closed_loop(hdr)
        row_vals['lgs_wfs_rate'] = inst.get_lgs_wfs_rate(hdr) * u.Hz

        # Initialize TT variables to masked to ensure they always exist in the schema
        row_vals['OSIRIS_tt_sensor'] = np.ma.masked
        row_vals['tt_wfs_rate'] = np.ma.masked
        row_vals['tt_wfs_centroid_gain'] = np.ma.masked

        if inst.get_instrument_name(hdr) == "OSIRIS": # Get tip tilt / NGS sensor information, which can differ on OSIRIS/Keck I
            TT_sensor = inst.get_tip_tilt_wfs_name(hdr) # Either STRAP or NIRTTS (TRICK)
            row_vals['OSIRIS_tt_sensor'] = TT_sensor
            if TT_sensor == "STRAP":
                tt_wfs_integration_time = inst.get_STRAP_integration_time(hdr) * u.ms
                row_vals['tt_wfs_rate'] = (1.0 / tt_wfs_integration_time).to(u.Hz)
                row_vals['tt_wfs_centroid_gain'] = inst.get_STRAP_centroid_gain(hdr)
            elif TT_sensor == "NIRTTS":
                # Not sure how to get the frame rate for TRICK
                row_vals['tt_wfs_centroid_gain'] = inst.get_TRICK_centroid_gain(hdr)
        
        row_vals['lgs_rms_wfe'] = inst.get_lgs_rms_wfe(hdr) * u.nm
        row_vals['lgs_layer_alt'] = inst.get_lgs_layer_altitude(hdr) * u.m
        row_vals['lbwfs_fwhm'] = _safe_get(inst.get_lbwfs_fwhm, hdr, u.arcsec)
        row_vals['ngs_wavelength'] = (inst.get_ngs_wavelength(hdr) * u.m).to(u.nm)
        row_vals['reconstructor_name'] = inst.get_reconstructor_name(hdr)
        row_vals['dm_gain'] = inst.get_dm_gain(hdr)
        row_vals['utt_gain'] = inst.get_utt_gain(hdr)
        row_vals['dtt_gain'] = inst.get_dtt_gain(hdr)
        row_vals['ao_mode'] = inst.get_ao_mode(hdr)
        row_vals['ao_hatch_open'] = inst.was_AO_hatch_open(hdr)

        ### Detector information ###
        row_vals['lgs_wfs_detector_gain'] = inst.get_lgs_wfs_detector_gain(hdr)

        ### Weather information ###
        # Early telemetry might not have all these parameters, we can 
        # get them from the weather files hopefully
        row_vals['humidity_dome'] = _safe_get(inst.get_dome_humidity, hdr)
        row_vals['T_dome_air'] = _safe_get(inst.get_dome_temperature, hdr, u.deg_C)
        row_vals['humidity_outside'] = _safe_get(inst.get_outside_humidity, hdr)
        row_vals['T_outside_air'] = _safe_get(inst.get_outside_temperature, hdr, u.deg_C)
        row_vals['P_barometric'] = _safe_get(inst.get_barometric_pressure, hdr, u.hPa)
        row_vals['P_barometric'] = row_vals['P_barometric'].to(u.Pa) if not np.ma.is_masked(row_vals['P_barometric']) else np.ma.masked
       
        row_vals['wind_direction'] = _safe_get(inst.get_wind_direction, hdr, u.deg)
        row_vals['wind_speed'] = _safe_get(inst.get_wind_speed, hdr, u.m / u.s)
        row_vals['t_weather_sample'] = _safe_get(inst.get_weather_sample_timestamp_string, hdr)
        row_vals['T_tube'] = _safe_get(inst.get_tube_temperature, hdr, u.deg_C)

        # Telemetry handling
        telem_files = read_image_telemetry(telemetry_file, verbose=True)
        ocam2k = load_ocam2k_data(telem_files)
        row_vals['num_lgs_wfs'] = 4 if has_four_lgs_data(ocam2k) else 1 # Either in LGS with 1 WFS or KAPA/LTAO with 4

        # Weather handling
        on_sky_fits_file_path = str(row_vals['image_path'])
        folder_with_on_sky_data = "/u/nstieg/work/ao/keck/massdimm/"
        on_sky_conditions = mu.estimate_on_sky_conditions(on_sky_fits_file_path, folder_with_on_sky_data)
        r0, turbulence_profile, wind_speed_profile, wind_direction_profile, _, _, _, _, tau0, theta0, _ = on_sky_conditions

        # Put weather data in row
        row_vals['r0'] = r0
        row_vals['turbulence_profile'] = turbulence_profile
        row_vals['wind_speed_profile'] = wind_speed_profile * (u.m / u.s)
        row_vals['wind_direction_profile'] = wind_direction_profile * u.deg
        row_vals['tau0'] = tau0 * u.s
        row_vals['theta0'] = theta0 * u.arcsec

        # Put in the table
        new_rows.append(row_vals)

    ### Make the table now from the list of dictionaries
    if not new_rows:
        return table if table is not None else new_empty_observation_table(n_rows=0)

    # Build an intermediate table from the new rows. By forcing masked=True, 
    # Astropy properly creates MaskedColumns for any np.ma.masked values.
    new_table = QTable(rows=new_rows, masked=True)
    
    # If a table was passed in, stack them vertically. 
    if table is not None:
        out = vstack([table, new_table], join_type='exact')
    else:
        out = new_table

    return out