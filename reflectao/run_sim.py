import numpy as np
import astropy.units as u
import os
from pathlib import Path
import paarti.utils.maos_utils as mu
import reflectao.telemetry_utils as tu
import reflectao.kapa_utils as ku
import reflectao.schema as schema

def run_maos_sim(hdr_tbl_row, seeds=[1]):
    ### Ensure that hdr_tbl_row is a row of the build_observation_table table
    schema.validate_table_has_schema(hdr_tbl_row, allow_extra_columns=True)

    ### Setup simulation parameters
    # Figure out the shortest time interval we need to simulate
    lgs_wfs_int_time = (1.0 / hdr_tbl_row['lgs_wfs_rate']).to(u.s)
    tt_int_time = (1.0 / hdr_tbl_row['tt_wfs_rate']).to(u.s) # STRAP or TRICK are tip tilt (tt) on Keck I. Assuming if TRICK is being used that no light goes to STRAP
    lbwfs_int_time = 15 * u.s # We set this equal to 15s because we simulate < 15s and so this doesn't matter. SHould be at least 10s (KAON 1303), with longer integrations for fainter guide stars
    max_frame_rate = max(hdr_tbl_row['lgs_wfs_rate'], hdr_tbl_row['tt_wfs_rate'])
    sim_dt = (1.0 / max_frame_rate).to(u.s) # Set sim_dt to fastest WFS readout rate in Hz
    sim_dtref = (1.0 / 472) * u.s # Fixed reference time step for KAPA; independent of sim_dt
    sim_end_steps = int((2.0 * u.s / sim_dt).to_value(u.dimensionless_unscaled))
    print(f"sim.dt = {sim_dt:.6f}  ({max_frame_rate:.1f} Hz max WFS frame rate)")
    print(f"sim.dtref = {sim_dtref:.6f}  (hardcoded to 472 Hz)")
    print(f"sim.end = {sim_end_steps} steps  (2.0 s at sim.dt)")

    # Set dtrats (ratio of sample period over dt, must be an integer)
    howfs_dtrat = (lgs_wfs_int_time / sim_dt).to_value(u.dimensionless_unscaled)
    tt_dtrat = (tt_int_time / sim_dt).to_value(u.dimensionless_unscaled)
    lbwfs_dtrat = (lbwfs_int_time / sim_dt).to_value(u.dimensionless_unscaled)
    powfs_dtrat_array = [int(howfs_dtrat), int(tt_dtrat), int(lbwfs_dtrat)] # To put in MAOS params

    # Check to make sure the dtrats are valid
    def check_dtrats(dtrat, name):
        assert dtrat.is_integer(), f"{name} must be an integer. Currently {dtrat}. sim-dt: {sim_dt} and {name}_int_time: {lgs_wfs_int_time if 'howfs' in name else tt_int_time if 'tt' in name else lbwfs_int_time} must be set such that their ratio is an integer."
        assert dtrat > 0, f"{name} must be positive. Currently {dtrat}"
    check_dtrats(howfs_dtrat, "howfs_dtrat")
    check_dtrats(tt_dtrat, "tt_dtrat")
    check_dtrats(lbwfs_dtrat, "lbwfs_dtrat")
    print(f"powfs.dtrat = {powfs_dtrat_array}  (HOWFS: {int(howfs_dtrat)}, TT: {int(tt_dtrat)}, LBWFS: {int(lbwfs_dtrat)})")

    # Set WFS names
    lgs_wfs_name = 'LGSWFS-OCAM2K' # For KAPA on Keck 1 5/5/26
    tt_wfs_name = ('TRICK' if 'trick' in str(hdr_tbl_row['OSIRIS_tt_sensor']).lower() else 'STRAP')

    # Set wavelengths for the WFSs. Note that if the science camera is in H-band, then TRICK uses K-band and vice versa
    powfs_wvl_array = [589e-9, 720e-9, 720e-9] # Note H band is ~1.6um and K is ~2.1um (1.6383e-6 & 2.1145e-6 meters for Hbb & Kp)
    if tt_wfs_name == 'TRICK':
        if hdr_tbl_row['wavelength'] > 2.1 * u.micron:
            tt_wfs_name = 'TRICK-H' # We're in K band (~2.12 micron) so TRICK is taking the H-band
            powfs_wvl_array[1] = 1.6383e-6  # H band
        else:
            tt_wfs_name = 'TRICK-K' # vice versa
            powfs_wvl_array[1] = 2.1145e-6  # K band
    print(f"LGS WFS: {lgs_wfs_name},  TT WFS: {tt_wfs_name}")
    print(f"powfs.wvl = [LGS: 589 nm, TT ({tt_wfs_name}): {powfs_wvl_array[1]*1e9:.1f} nm, LBWFS: {powfs_wvl_array[2]*1e9:.1f} nm]")

    # Get data from lbwfs header
    lbwfs_fwhm = hdr_tbl_row['lbwfs_fwhm'].to_value(u.arcsec)
    ngs_guide_star_magnitude = hdr_tbl_row['tt_gs_r_mag']
    print(f"LBWFS FWHM: {lbwfs_fwhm:.3f} arcsec,  NGS R mag: {ngs_guide_star_magnitude:.2f}")

    # Define atmospheric properties
    brooke_atm_data_heights = np.array([100.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 16000.0]) # Heights that Brooke's code gets the atmospheric profiles at
    r0z = hdr_tbl_row['r0'].to_value(u.m) # <-- Estimated from atmospheric telemetry
    print(f"r0 = {r0z:.4f} m")

    # Calculate signal level (siglev), background (bkgrnd) and noise equivalent angle
    # parameters
    # NOTE: 8.1 and 10.4 kinda suck as estimates which is why we try to use telemetry if available
    # Luckily this only effects the HOWFS (LGS WFS) since we measure the tip tilt star magnitude
    # with the LBWFS
    num_LGS = hdr_tbl_row['num_lgs_wfs']
    if num_LGS == 1:
        lgs_guide_star_magnitude = 8.1 # Brooke has 8.1 for SLGS. Maybe this would be brighter now?
    else:
        lgs_guide_star_magnitude = 10.4 # For 4 LGS, the guide star is dimmer since the light is split.
                                        # Based on some very rough estimates from averaged telemetry.
                                        # See /Users/nstieg/work/ao/keck/kapa/telemetry/ocam2k/keck_nea_photon_investigation.ipynb
    print(f"num_LGS = {num_LGS},  LGS estimated magnitude: {lgs_guide_star_magnitude}")

    # Get nearecon, siglev, and background (scaled by dtref)
    _, howfs_nearecon, howfs_siglev, howfs_bkgrnd = mu.keck_nea_photons(lgs_guide_star_magnitude, lgs_wfs_name,
                                                                        r0z, sim_dtref.to_value(u.s))
    _, tt_nearecon, tt_siglev, tt_bkgrnd = mu.keck_nea_photons(ngs_guide_star_magnitude, tt_wfs_name,
                                                                        r0z, sim_dtref.to_value(u.s))
    _, lbwfs_nearecon, lbwfs_siglev, lbwfs_bkgrnd = mu.keck_nea_photons(ngs_guide_star_magnitude, 'LBWFS',
                                                                        r0z, sim_dtref.to_value(u.s),
                                                                        lbwfs_fwhm=lbwfs_fwhm)

    # Noise equivalent angle. Not used for physical optics, only geometric. In milliarcseconds
    # Since the tt-wfs is the only one that uses geometric optics this really only matters
    # for it
    nearecon = [howfs_nearecon, tt_nearecon, lbwfs_nearecon]
    print(f"nearecon (mas):          HOWFS={howfs_nearecon:.3f},  TT={tt_nearecon:.3f},  LBWFS={lbwfs_nearecon:.3f}")

    # Signal level (siglev) is in number of photons from the star per subaperture per dtref step (dtref ~ dt)
    # powfs.siglev is photons/subaperture/dtref step
    # can do powfs.siglevs (one for each WFS)
    # NOTE: Gets changed if we're on KAPA based on telemetry
    siglev = [howfs_siglev, tt_siglev, lbwfs_siglev]
    individual_siglevs = None # Set this to a list of the individual signal levels for each WFS if we want to use those instead of the same siglev for all WFSs of the same type
    print(f"siglev (ph/subap/dtref): HOWFS={howfs_siglev:.1f},  TT={tt_siglev:.1f},  LBWFS={lbwfs_siglev:.1f}")

    # powfs.bckgrnd is background in e-/pix/dtref step
    # Maybe snag from telemetry eventually if we can
    # NOTE: Gets changed if we're on KAPA based on telemetry
    bkgrnd = [howfs_bkgrnd, tt_bkgrnd, lbwfs_bkgrnd] # bkgrnd is in number of background photons per pixel within subaperture
    individual_bkgrnd = None
    print(f"bkgrnd (e-/px/dtref):    HOWFS={howfs_bkgrnd:.4f},  TT={tt_bkgrnd:.4f},  LBWFS={lbwfs_bkgrnd:.4f}")

    # Make adjustments
    if lgs_wfs_name == 'LGSWFS-OCAM2K':
        OCAM2K_QE = 0.88 # Noted here, probably not used because we only care to simulate the photoelectrons we see

        # Put in signal level measured from telemetry
        if num_LGS == 4: # LTAO
            # See if we could get the signal & background level for each WFS from telemetry, otherwise use estimates
            try:
                # If using individual signal levels per WFS
                individual_siglevs =  list(hdr_tbl_row['lgs_wfs_signal_levels'].to_value(u.electron)) + siglev[1:] # Get [siglev_wfs1, wfs2, wfs3, wfs4, tt_siglev, lbwfs_siglev]

                # If using individual background levels per WFS NOTE: Not supported in MAOS as of 5/15/26
                # individual_bkgrnd = list(hdr_tbl_row['lgs_wfs_background_levels'].to_value(u.electron)) + bkgrnd[1:] # Get [bkgrnd_wfs1, wfs2, wfs3, wfs4, tt_bkgrnd, lbwfs_bkgrnd]
                bkgrnd[0] = hdr_tbl_row['avg_lgs_wfs_background_level'].to_value(u.electron)
                print(f"Using per-WFS siglevs from telemetry: {[f'{s:.1f}' for s in individual_siglevs]}")
                print(f"Using avg LGS background from telemetry: bkgrnd[0] = {bkgrnd[0]:.4f} e-/px/dtref")
            except Exception as e:
                print(f"Error getting individual signal and background levels from telemetry: {e}. Using estimates instead.")

                # Use estimated average
                signal_epers = 75e3 * (u.electron / u.second) # electrons per second per subaperture on average measured from fig. 14: https://docs.google.com/document/d/1YQNdJjWz2LFtnWaDdaJxFfnd-Dfjs9zw7nqB67q1gHM/edit?tab=t.0
                signal_photons = signal_epers * sim_dtref
                siglev[0] = signal_photons.to_value(u.photon)

                # 0.20 e-/px/frame estimated from the average sigma_subaperture values in fig. 14: https://docs.google.com/document/d/1YQNdJjWz2LFtnWaDdaJxFfnd-Dfjs9zw7nqB67q1gHM/edit?tab=t.0#heading=h.scqamqs3mj7m
                bkgrnd[0] = 0.2 * (sim_dtref / (2 * u.ms)) # 0.20 is at 500Hz (2ms) so scale to sim_dtref
                print(f"Using estimated LGS siglev = {siglev[0]:.1f} ph/subap/dtref")
                print(f"Using estimated LGS bkgrnd[0] = {bkgrnd[0]:.4f} e-/px/dtref")

            # Locate tip tilt star and laser guide stars
            tt_x = lbwfs_x = hdr_tbl_row['tt_star_offset_x'].to_value(u.arcsec) # in arcsec -> note the LBWFS looks at the same NGS the tip tilt is looking at
            tt_y = lbwfs_y = hdr_tbl_row['tt_star_offset_y'].to_value(u.arcsec) # in arcsec

            # The LGS are on a 10 arcsec radius asterism, with a 7'' offset in the x direction to match on-sky KAPA
            # Note the orderings of the WFSs matches the convention in https://docs.google.com/document/d/1YQNdJjWz2LFtnWaDdaJxFfnd-Dfjs9zw7nqB67q1gHM/edit?tab=t.0#heading=h.scqamqs3mj7m
            wfs_thetax = f"[12.5 12.5 1.5 1.5 {tt_x} {lbwfs_x}]" # wfs.thetax coordinate in arcsec units
            wfs_thetay = f"[5.5 -5.5 -5.5 5.5 {tt_y} {lbwfs_y}]" # wfs.thetay
            print(f"TT/LBWFS star offset: ({tt_x:.3f}\", {tt_y:.3f}\") arcsec from on-axis")
            print(f"wfs.thetax = {wfs_thetax},  wfs.thetay = {wfs_thetay}")
        else:
            # TODO: Get estimates of the single LGS (non-KAPA) background & signal level

            # For 1 LGS, the LGS is on top of the tip tilt star, so the tt & LGS constellation are at the same place: [0 tt_x lbwfs_x] and same for y
            tt_x = lbwfs_x = hdr_tbl_row['tt_star_offset_x'] # in arcsec -> note the LBWFS looks at the same NGS the tip tilt is looking at
            tt_y = lbwfs_y = hdr_tbl_row['tt_star_offset_y'] # in arcsec
            wfs_thetax = f"[0 {tt_x} {lbwfs_x}]"
            wfs_thetay = f"[0 {tt_y} {lbwfs_y}]"
            print(f"TT/LBWFS star offset: ({tt_x}\", {tt_y}\") arcsec from on-axis")
            print(f"wfs.thetax = {wfs_thetax},  wfs.thetay = {wfs_thetay}")

    # Setup final sim parameters
    if num_LGS == 4:
        conf_name = "A_keck_kapa_compare_sky_4lgs_template.conf"
    elif num_LGS == 1:
        conf_name = "A_keck_kapa_compare_sky_1lgs_template.conf"
    print(f"Config: {conf_name}")

    # Final maos command
    # TODO:
    # - recon.glao = 1 if GLAO mode (average gradient from WFS of same type)
    # - define which subapertures are lit or not based on telemetry for each wfs (which feature in MAOS allows this?)
    # - Get siglev & bkgrnd for single lgs mode from telemetry

    # cd to the directory with the configurations
    maos_config_env = os.environ.get("MAOS_CONFIG_PATH")
    if maos_config_env is None:
        raise EnvironmentError(
            "MAOS_CONFIG_PATH environment variable is not set. "
            "Set it to the directory containing the MAOS .conf files, e.g.:\n"
            "  export MAOS_CONFIG_PATH=/path/to/keck_maos/kapa/"
        )
    maos_configuration_path = Path(maos_config_env)
    os.chdir(maos_configuration_path)

    pp = print_array_maos_style # So we have array sprinting as [1 2 3] instead of [1, 2, 3] which is what MAOS expects

    image_path = Path(hdr_tbl_row['image_path'])
    fits_filename = image_path.name
    night = tu.get_night_from_fits_file_path(image_path)
    out_dir = str(tu.get_data_path() / night / "maos" / fits_filename / conf_name[:-5]) + "/"
    print(f"Seeds {pp(seeds)}: output -> {out_dir}")

    # Create the output directory if it doesn't exist
    out_dir_path = Path(out_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)

    # Get some data from the header table row
    atm_wt = hdr_tbl_row['turbulence_profile']
    atm_ws = hdr_tbl_row['wind_speed_profile'].to_value(u.m / u.s)
    atm_wddeg = hdr_tbl_row['wind_direction_profile'].to_value(u.deg)
    zenith_angle = hdr_tbl_row['zenith_angle']
    na_layer_height = hdr_tbl_row['lgs_layer_alt'].to_value(u.m)
    print(f"zenith angle = {zenith_angle:.2f},  NA layer height = {na_layer_height:.1f} m")
    print(f"r0 = {r0z:.4f} m,  atm.wt = {atm_wt},  atm.ws = {atm_ws} m/s,  atm.wddeg = {atm_wddeg} deg")

    # Check if we're using individual signal levels and backgrounds based on telemetry or the fallback estimates
    # If we have individual levels, we can set them with wfs.siglev & bkgrnd, otherwise we use powfs.siglev & background
    # to set them for the whole group together (e.g. all 4 LGS WFSs get the same siglev if using powfs.siglev, but can have different siglevs if using wfs.siglev)
    siglev_command_part = f"wfs.siglev={pp(individual_siglevs)}" if individual_siglevs is not None else f"powfs.siglev={pp(siglev)}"
    bkgrnd_command_part = f"wfs.bkgrnd={pp(individual_bkgrnd)}" if individual_bkgrnd is not None else f"powfs.bkgrnd={pp(bkgrnd)}"
    print(f"{siglev_command_part}")
    print(f"{bkgrnd_command_part}")

    # Final maos command — pass the full seeds list; MAOS loops over seeds internally
    maos_cmd = str(f"maos -c {conf_name} -o {out_dir}"
                f" sim.seeds={pp(seeds)} sim.zadeg={zenith_angle.to_value(u.deg)}" # include spaces as they're not added automatically
                f" sim.dt={sim_dt.to_value(u.s)} sim.dtref={sim_dtref.to_value(u.s)} sim.end={sim_end_steps}"
                f" {siglev_command_part} {bkgrnd_command_part} powfs.nearecon={nearecon}"
                f" powfs.dtrat={pp(powfs_dtrat_array)} powfs.wvl={pp(powfs_wvl_array)}"
                f" wfs.thetax={wfs_thetax} wfs.thetay={wfs_thetay}"
                f" atm.r0z={r0z} atm.wt={atm_wt} atm.ws={atm_ws}"
                f" atm.wddeg={atm_wddeg} atm.ht={pp(brooke_atm_data_heights)}"
                f" powfs.hs=[{na_layer_height} inf inf]"
                f" -O")
    print(f"Running MAOS:\n  {maos_cmd}")
    os.system(maos_cmd)

def print_array_maos_style(list_or_array):
    return "[" + " ".join([str(x) for x in list_or_array]) + "]"