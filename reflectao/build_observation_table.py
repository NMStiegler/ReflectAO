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

        # Check this is an instrument we support
        keck_ao_instruments = ["OSIRIS", "NIRC2"] # We can add more instruments here as we support them
        if hdr["CURRINST"] not in keck_ao_instruments:
            print(f"Warning: Instrument {hdr['CURRINST']} not in supported instrument list {keck_ao_instruments}. Skipping file {p}.")
            row_vals['instrument_name'] = hdr["CURRINST"] # We can still fill in the instrument name column for unsupported instruments, which is useful for debugging and expanding support in the future.

        ############# Add metadata from this row's FITS file #############
        ### System information ###
        row_vals["image_path"] = str(p) # Path to the FITS file
        row_vals['telemetry_file_path'] = str(telemetry_file) if (telemetry_file is not np.ma.masked and telemetry_file is not None) else np.ma.masked
        
        ### Telescope / instrument / set information ###
        row_vals['telescope_name'] = inst.get_telescope_name(hdr)
        if inst.get_instrument_name(hdr) == hdr["CURRINST"]:
            row_vals['instrument_name'] = inst.get_instrument_name(hdr)
        else:
            print(f"Warning: Instrument name from KAI ({inst.get_instrument_name(hdr)}) does not match CURRINST in header ({hdr['CURRINST']}). Using CURRINST value.")
            row_vals['instrument_name'] = hdr["CURRINST"]
        row_vals['instrument_angle'] = inst.get_instrument_angle(hdr) * u.deg
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
        row_vals['propagating'] = inst.was_laser_propagating(hdr)
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
            elif TT_sensor == "NIRTTS + STRAP":
                # There's some way to tell which one is being used for tip-tilt
                # but I'm not sure how
                pass

        # Find where the tip-tilt star is on the detector using the aotsx/y values, which define 
        # where the tip-tilt stage is located and thus where the tip-tilt star is located since the tip-tilt stage steers the tip-tilt star onto the TT sensor.
        # Based on: kai/reduce/kai_util.py: aotsxy2pix
        ref_aotsx = 1.820 * u.mm  # hard-code
        ref_aotsy = -11.130 * u.mm # hard-code
        tt_stage_scale = 0.727 * (u.mm / u.arcsec) # maybe for NIRC2 only? Should be 0.725? https://www2.keck.hawaii.edu/observing/kecktelgde/ktelinstupdate.pdf
        row_vals['aotsx'] = _safe_get(inst.get_aotsx, hdr, u.mm)
        row_vals['aotsy'] = _safe_get(inst.get_aotsy, hdr, u.mm)
        dx = (row_vals['aotsx'] - ref_aotsx) / tt_stage_scale
        dy = (row_vals['aotsy'] - ref_aotsy) / tt_stage_scale
        row_vals['tip_tilt_radial_offset'] = np.sqrt(dx**2 + dy**2)
        cosa = np.cos(np.radians(-row_vals['instrument_angle']))
        sina = np.sin(np.radians(-row_vals['instrument_angle']))
        rot_matrix = np.array([[cosa, sina], [-sina, cosa]])
        unrotated_position_vector = np.array([dx.to_value(u.arcsec), dy.to_value(u.arcsec)]) * u.arcsec
        rotated_position_vector = rot_matrix.dot(unrotated_position_vector)
        row_vals['tt_star_offset_x'] = rotated_position_vector[0]
        row_vals['tt_star_offset_y'] = rotated_position_vector[1]

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
        row_vals['tt_gs_r_mag'] = get_tt_guide_star_r_mag(p)

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
        if has_four_lgs_data(ocam2k):
            row_vals['num_lgs_wfs'] = 4
        elif has_only_single_lgs_data(ocam2k):
            row_vals['num_lgs_wfs'] = 1
        else:
            row_vals['num_lgs_wfs'] = np.ma.masked # We don't have telemetry for this image, so we can't determine the number of LGS WFSs. This could happen if the telemetry is missing or corrupted.

        # Calculate signal and background levels for each LGS WFS based on telemetry, if possible. If not, leave as masked.
        if has_four_lgs_data(ocam2k):
            # NOTE: Using convention where I append the unit and then shape of arrays to the end of variable names.
            # So arr_s_X_Y means it's size (X, Y) and unit seconds
            # Get the mean electrons read per subaperture per frame for each subaperture across all 4 WFSs as well as the standard
            # deviation of the electrons read per subaperture per frame for each subaperture across all 4 WFSs. Both should have shape (4, 304)
            # and units of electrons per read although the units of the stds is really sqrt(e-^2 variance)
            sensor_mean_electrons_4_304, sensor_stds_electrons_4_304 = compute_aperture_wise_electron_stats(ocam2k, row_vals)

            # Break the subapertures into a bright and dim population based on the 20th percentile signal level
            # We do this because empirically on all of the telemetry we've looked at, there are some 'less illuminated'
            # subapertures around the edges of the subap mask that don't have nearly as many electrons as the more fully
            # illuminated subapertures and we want to make sure we treat the ones which get electrons as having a signal
            # and the ones which don't as not having a signal. The 20th percentile has worked well before
            thresh = 20 # Use 20th percentile
            
            # Get the indices of the lit and unlit subapertures 
            row_vals['unlit_subap_indices'] = get_unlit_data_indices_for_all_wfs(sensor_mean_electrons_4_304, thresh=thresh)
            row_vals['lit_subap_indices'] = get_lit_data_indices_for_all_wfs(sensor_mean_electrons_4_304, thresh=thresh)

            # Calculate the median signal level from the dim and bright subapertures respectively for each WFS
            _, median_bright_mean_e_4 = compute_median_values_for_all_wfs(sensor_mean_electrons_4_304, thresh=thresh, unlit_data_indices_for_all_wfs=row_vals['unlit_subap_indices'], lit_data_indices_for_all_wfs=row_vals['lit_subap_indices'])
            median_bright_mean_e_4 = turn_list_of_quantities_into_quantity_array(median_bright_mean_e_4)

            # Calculate the standard deviation of signal levels from the dim and bright subapertures respectively for each WFS
            # This tells us about the noise of each pixel
            median_dim_std_e_4, _ = compute_median_values_for_all_wfs(sensor_stds_electrons_4_304, thresh=20, unlit_data_indices_for_all_wfs=row_vals['unlit_subap_indices'], lit_data_indices_for_all_wfs=row_vals['lit_subap_indices'])
            median_dim_std_e_4 = turn_list_of_quantities_into_quantity_array(median_dim_std_e_4)

            # Signal level is in electrons (electrons = photons for us) per subaperture per frame
            row_vals['lgs_wfs_signal_levels'] = median_bright_mean_e_4
            row_vals['avg_lgs_wfs_signal_level'] = np.mean(median_bright_mean_e_4)

            # Compute background levels based on the dim subapertures. We want these values in number of electrons per pixel per
            # frame coming from the background. Each subaperture has 4x4 = 16 pixels. For each pixel of a dim subap, we expect the total
            # noise sigma^2_pixel = N_electrons + RN^2 where N_electrons is the number of background sky and dark electrons and RN is the read noise in electrons.
            # We use N_electrons instead of sigma^2_electrons because we expect the these electrons to have Poisson statistics where sigma^2_electrons = N_electrons
            # We assume RN = 0.36 electrons based on KAON 1337. Thus for the whole subap, sigma^2_subap = 16 * (N_electrons + RN^2)
            # sigma_subap, the total noise in each subap, is in median_dim_std_e_4 for each WFS. Now we can rearrange the equation to solve for N_electrons:
            # N_electrons = ((sigma^2_subap / 16) - RN^2). We could do this / F^2 where F is the extra noise factor from the OCAM2k EMCCD gain, but since we're directly
            # measuring N_electrons we don't have to care about this - it's baked into our measurement of sigma^2_subap 
            RN = 0.36 * u.electron # Read noise
            pix_per_subap_side = 4
            pix_per_subap = pix_per_subap_side**2
            row_vals['lgs_wfs_background_levels'] = ((median_dim_std_e_4**2 / pix_per_subap) - RN**2).to_value(u.electron * u.electron) * u.electron # Because of Poisson, units of N = units of N^2 :|
            row_vals['avg_lgs_wfs_background_level'] = np.mean(row_vals['lgs_wfs_background_levels'])
        elif has_only_single_lgs_data(ocam2k):
            # NOTE: Using convention where I append the unit and then shape of arrays to the end of variable names.
            # So arr_s_X_Y means it's size (X, Y) and unit seconds
            # Get the mean electrons read per subaperture per frame for each subaperture for our one WFS as well as the standard
            # deviation of the electrons read per subaperture per frame for each subaperture for our one WFS. Both should have shape (304)
            # and units of electrons per read although the units of the stds is really sqrt(e-^2 variance)
            sensor_mean_electrons_1_304, sensor_stds_electrons_1_304 = compute_aperture_wise_electron_stats(ocam2k, row_vals, num_WFS=1)

            # Break the subapertures into a bright and dim population based on the 20th percentile signal level
            # We do this because empirically on all of the telemetry we've looked at, there are some 'less illuminated'
            # subapertures around the edges of the subap mask that don't have nearly as many electrons as the more fully
            # illuminated subapertures and we want to make sure we treat the ones which get electrons as having a signal
            # and the ones which don't as not having a signal. The 20th percentile has worked well before
            thresh = 20 # Use 20th percentile
            
            # Get the indices of the lit and unlit subapertures 
            row_vals['unlit_subap_indices'] = get_unlit_data_indices_for_all_wfs(sensor_mean_electrons_1_304, thresh=thresh, num_WFS=1)
            row_vals['lit_subap_indices'] = get_lit_data_indices_for_all_wfs(sensor_mean_electrons_1_304, thresh=thresh, num_WFS=1)

            # Calculate the median signal level from the dim and bright subapertures respectively for each WFS
            _, median_bright_mean_e_1 = compute_median_values_for_all_wfs(sensor_mean_electrons_1_304, thresh=thresh, unlit_data_indices_for_all_wfs=row_vals['unlit_subap_indices'], lit_data_indices_for_all_wfs=row_vals['lit_subap_indices'], num_WFS=1)
            median_bright_mean_e_1 = turn_list_of_quantities_into_quantity_array(median_bright_mean_e_1)

            # Calculate the standard deviation of signal levels from the dim and bright subapertures respectively for each WFS
            # This tells us about the noise of each pixel
            median_dim_std_e_1, _ = compute_median_values_for_all_wfs(sensor_stds_electrons_1_304, thresh=20, unlit_data_indices_for_all_wfs=row_vals['unlit_subap_indices'], lit_data_indices_for_all_wfs=row_vals['lit_subap_indices'], num_WFS=1)
            median_dim_std_e_1 = turn_list_of_quantities_into_quantity_array(median_dim_std_e_1)

            # Signal level is in electrons (electrons = photons for us) per subaperture per frame
            row_vals['lgs_wfs_signal_levels'] = median_bright_mean_e_1
            row_vals['avg_lgs_wfs_signal_level'] = np.mean(median_bright_mean_e_1)

            # Compute background levels based on the dim subapertures. We want these values in number of electrons per pixel per
            # frame coming from the background. Each subaperture has 4x4 = 16 pixels. For each pixel of a dim subap, we expect the total
            # noise sigma^2_pixel = N_electrons + RN^2 where N_electrons is the number of background sky and dark electrons and RN is the read noise in electrons.
            # We use N_electrons instead of sigma^2_electrons because we expect these electrons to have Poisson statistics where sigma^2_electrons = N_electrons
            # We assume RN = 0.36 electrons based on KAON 1337. Thus for the whole subap, sigma^2_subap = 16 * (N_electrons + RN^2)
            # sigma_subap, the total noise in each subap, is in median_dim_std_e_1 for each WFS. Now we can rearrange the equation to solve for N_electrons:
            # N_electrons = ((sigma^2_subap / 16) - RN^2). We could do this / F^2 where F is the extra noise factor from the OCAM2k EMCCD gain, but since we're directly
            # measuring N_electrons we don't have to care about this - it's baked into our measurement of sigma^2_subap 
            RN = 0.36 * u.electron # Read noise
            pix_per_subap_side = 4
            pix_per_subap = pix_per_subap_side**2
            row_vals['lgs_wfs_background_levels'] = ((median_dim_std_e_1**2 / pix_per_subap) - RN**2).to_value(u.electron * u.electron) * u.electron # Because of Poisson, units of N = units of N^2 :|
            row_vals['avg_lgs_wfs_background_level'] = np.mean(row_vals['lgs_wfs_background_levels'])
        
        # Weather handling
        on_sky_fits_file_path = str(row_vals['image_path'])
        folder_with_on_sky_data = "/u/nstieg/work/ao/keck/massdimm/"
        if row_vals['instrument_name'] in keck_ao_instruments:
            on_sky_conditions = mu.estimate_on_sky_conditions(on_sky_fits_file_path, folder_with_on_sky_data)
            r0, turbulence_profile, wind_speed_profile, wind_direction_profile, _, _, _, _, tau0, theta0, _ = on_sky_conditions

            # Put weather data in row
            row_vals['r0'] = r0 # already in meters
            row_vals['turbulence_profile'] = turbulence_profile
            row_vals['wind_speed_profile'] = wind_speed_profile * (u.m / u.s)
            row_vals['wind_direction_profile'] = wind_direction_profile * u.deg
            row_vals['tau0'] = tau0 * u.s
            row_vals['theta0'] = theta0 * u.arcsec

            # Fix any rows which may be uneven since if a column of the table has a list of arrays, all the arrays need to be the same length (so I'll turn them into 2d arrays)
            row_vals['unlit_subap_indices'] = fill_in_uneven_list_of_arrays(row_vals['unlit_subap_indices']) if row_vals['unlit_subap_indices'] is not np.ma.masked else np.ma.masked
            row_vals['lit_subap_indices'] = fill_in_uneven_list_of_arrays(row_vals['lit_subap_indices']) if row_vals['lit_subap_indices'] is not np.ma.masked else np.ma.masked

        # Put in the table
        new_rows.append(row_vals)

    ### Make the table now from the list of dictionaries
    if not new_rows:
        return table if table is not None else new_empty_observation_table(n_rows=0)

    # Pad variable-length columns to a uniform length across all rows before
    # constructing the QTable. fill_in_uneven_list_of_arrays already pads within
    # each row, but rows can still disagree on the padded width if different
    # exposures have different subaperture counts or atmospheric layer counts.
    for _col in ['turbulence_profile', 'wind_speed_profile', 'wind_direction_profile',
                 'unlit_subap_indices', 'lit_subap_indices']:
        pad_uneven_column_across_rows(new_rows, _col)

    # Build an intermediate table from the new rows. By forcing masked=True,
    # Astropy properly creates MaskedColumns for any np.ma.masked values.
    new_table = QTable(rows=new_rows, masked=True)
    
    # If a table was passed in, stack them vertically. 
    if table is not None:
        out = vstack([table, new_table], join_type='exact')
    else:
        out = new_table

    return out


def turn_list_of_quantities_into_quantity_array(list_of_quantities):
    """
    Helper function to turn a list of astropy Quantities into a single Quantity with an array value.

    :param list_of_quantities: List of astropy Quantities, all with the same unit.
    :type list_of_quantities: astropy.units.quantity.Quantity array-like
    :return: An astropy Quantity with an array value containing the values from the input list, all converted to the specified unit.
    :rtype: astropy.units.Quantity
    """
    assert list_of_quantities, "Input list of quantities is empty"
    assert( hasattr(list_of_quantities[0], 'unit') ), "Elements of input list must be astropy Quantities with units"
    unit = list_of_quantities[0].unit
    values = [q.to(unit).value for q in list_of_quantities]
    return np.array(values) * unit

def fill_in_uneven_list_of_arrays(list_of_arrays):
    """
    Helper function to turn a list of 1D numpy arrays of different lengths into a 2D numpy array where the shorter arrays are filled with np.nan values at the end to match the length of the longest array.

    :param list_of_arrays: List of 1D numpy arrays of different lengths.
    :type list_of_arrays: list of numpy.ndarray
    :return: A 2D numpy array where the shorter arrays are filled with np.nan values at the end to match the length of the longest array.
    :rtype: numpy.ndarray
    """
    max_length = max(len(arr) for arr in list_of_arrays)
    filled_array = np.full((len(list_of_arrays), max_length), np.nan)
    for i, arr in enumerate(list_of_arrays):
        filled_array[i, :len(arr)] = arr
    return filled_array


def pad_uneven_column_across_rows(rows, col):
    """
    Pad a variable-length array column to a uniform last-dimension length across
    all rows so that QTable can store it as a proper N-D column. Works for both
    1D arrays (e.g. profile columns) and 2D arrays (e.g. subap index columns with
    shape (num_WFS, n_subaps) where n_subaps varies). Rows with masked values are
    left untouched. Astropy Quantity units are preserved.

    :param rows: List of row dictionaries, as produced before QTable construction.
    :type rows: list of dict
    :param col: Column name whose values are arrays with a variable last dimension.
    :type col: str
    """
    live = [(i, r[col]) for i, r in enumerate(rows) if r[col] is not np.ma.masked]
    if not live:
        return
    # For 2D arrays the varying dimension is axis=-1 (columns); for 1D it's also axis=-1
    max_len = max(np.shape(arr)[-1] for _, arr in live)
    for i, arr in live:
        arr_np = arr.value if hasattr(arr, 'unit') else np.asarray(arr)
        if arr_np.shape[-1] < max_len:
            pad_width = [(0, 0)] * (arr_np.ndim - 1) + [(0, max_len - arr_np.shape[-1])]
            padded_np = np.pad(arr_np.astype(float), pad_width, constant_values=np.nan)
            rows[i][col] = padded_np * arr.unit if hasattr(arr, 'unit') else padded_np