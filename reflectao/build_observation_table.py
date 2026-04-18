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

from .schema import new_empty_observation_table, validate_table_has_schema


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
    telemetry_path=None,
    weather_path=None,
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
    :param telemetry_path: Optional path to already-downloaded telemetry.
        TODO: auto-discovery / downloading is not implemented in this slice.
    :type telemetry_path: str or pathlib.Path, optional
    :param weather_path: Optional path to already-downloaded weather data.
        TODO: auto-discovery / downloading is not implemented in this slice.
    :type weather_path: str or pathlib.Path, optional
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

    # TODO: Implement loading telemetry and weather data 
    _ = telemetry_path
    _ = weather_path

    # Use KAI to read the data we want from our FITS headers
    inst = _get_kai_instrument(instrument)

    # Make sure we end up with a list of paths regardless of what
    # was passed in
    paths = _normalize_paths(fits_paths)
    
    # Setup the table to add to. If no table is provided, create a new one.
    if table is None:
        out = new_empty_observation_table(n_rows=0)
    else:
        validate_table_has_schema(table)
        out = table

    # Add a row to the table for each image, pointed to by a path in our list
    for p in paths:
        # Get the fits header for the current image
        hdr = fits.getheader(p)

       # Set up a dictionary to hold the values we'll add for this row. Start
       # everything as a masked value just to be safe
        row_vals = {name: np.ma.masked for name in out.colnames}

        ############# Add metadata from this row's FITS file #############
        # System information
        row_vals["image_path"] = str(p) # Path to the FITS file
        
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
        row_vals['wavelength'] = (inst.get_central_wavelength(hdr) * u.um).to(u.nm)
        row_vals['filter_name'] = inst.get_filter_name(hdr)
        row_vals['airmass'] = inst.get_airmass(hdr)
        row_vals['zenith_angle'] = (np.arccos(1.0 / float(hdr['AIRMASS'])) * u.rad).to(u.deg)
        row_vals['aborted'] = not 'f' in inst.was_exposure_aborted(hdr).strip(" ").lower()
        row_vals['dither_name'] = inst.get_dither_name(hdr)

        ### AO system information ###
        row_vals['lgs_wfs_rate'] = inst.get_lgs_wfs_rate(hdr) * u.Hz
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
        try: row_vals['ngs_fwhm'] = inst.get_ngs_fwhm(hdr) * u.arcsec # Early telemetry might not have this
        except KeyError: pass
        row_vals['ngs_wavelength'] = (inst.get_ngs_wavelength(hdr) * u.m).to(u.nm)
        row_vals['reconstructor_name'] = inst.get_reconstructor_name(hdr)
        row_vals['dm_gain'] = inst.get_dm_gain(hdr)
        row_vals['utt_gain'] = inst.get_utt_gain(hdr)
        row_vals['dtt_gain'] = inst.get_dtt_gain(hdr)
        row_vals['ao_mode'] = inst.get_ao_mode(hdr)

        ### Detector information ###
        row_vals['lgs_wfs_detector_gain'] = inst.get_lgs_wfs_detector_gain(hdr)

        ### Weather information ###
        # Early telemetry might not have all these parameters, we can 
        # get them from the weather files hopefully
        try: row_vals['humidity_dome'] = inst.get_dome_humidity(hdr)
        except KeyError: pass
        try: row_vals['T_dome_air'] = inst.get_dome_temperature(hdr) * u.deg_C
        except KeyError: pass
        try: row_vals['humidity_outside'] = inst.get_outside_humidity(hdr)
        except KeyError: pass
        try: row_vals['T_outside_air'] = inst.get_outside_temperature(hdr) * u.deg_C
        except KeyError: pass
        try: row_vals['P_barometric'] = (inst.get_barometric_pressure(hdr) * u.hPa).to(u.Pa)
        except KeyError: pass
        try: row_vals['wind_direction'] = inst.get_wind_direction(hdr) * u.deg # <-- What reference angle?
        except KeyError: pass
        try: row_vals['wind_speed'] = inst.get_wind_speed(hdr) * u.m / u.s # <-- What unit is this logged as?
        except KeyError: pass
        try: row_vals['t_weather_sample'] = inst.get_weather_sample_timestamp_string(hdr)
        except KeyError: pass
        try: row_vals['T_tube'] = inst.get_tube_temperature(hdr) * u.deg_C
        except KeyError: pass

        # Calculate atm parameters
        # fried, turbpro, windspds, winddrcts, _, _, _, _, _, _ = estimate_on_sky_conditions(sky_file, 
        #                                                                                    sky_folder)

        # Put in the table
        out.add_row(row_vals)

    return out