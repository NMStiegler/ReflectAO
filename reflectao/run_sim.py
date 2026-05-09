import numpy as np
import astropy.units as u
import os
import paarti.utils.maos_utils as mu
import reflectao.telemetry_utils as tu
import reflectao.kapa_utils as ku
import reflectao.schema as schema

def run_maos_sim(hdr_tbl_row):
    ### Ensure that hdr_tbl_row is a row of the build_observation_table table
    schema.validate_table_has_schema(hdr_tbl_row, allow_extra_columns=True)

    ### Setup simulation parameters
    # Figure out the shortest time interval we need to simulate
    shwfs_int_time = (1.0 / hdr_tbl_row['lgs_wfs_rate']).to(u.s)
    tt_int_time = (1.0 / hdr_tbl_row['tt_wfs_rate']).to(u.s) # STRAP or TRICK are tip tilt (tt) on Keck I. Assuming if TRICK is being used that no light goes to STRAP
    lbwfs_int_time = 15 * u.s # We set this equal to 15s because we simulate < 15s and so this doesn't matter. SHould be at least 10s (KAON 1303), with longer integrations for fainter guide stars
    max_frame_rate = max(hdr_tbl_row['lgs_wfs_rate'], hdr_tbl_row['tt_wfs_rate']) 
    sim_dt = (1.0 / max_frame_rate).to(u.s) # Set sim_dt to fastest WFS readout rate in Hz 
    
    # Set dtrats (ratio of sample period over dt, must be an integer)
    howfs_dtrat = sim_dt / shwfs_int_time
    tt_dtrat = sim_dt / tt_int_time
    lbwfs_dtrat = sim_dt / lbwfs_int_time

    # Check to make sure the dtrats are valid
    def check_dtrats(dtrat, name):
        assert dtrat.is_integer(), f"{name} must be an integer"
        assert dtrat > 0, f"{name} must be positive"
    check_dtrats(howfs_dtrat, "howfs_dtrat")
    check_dtrats(tt_dtrat, "tt_dtrat")
    check_dtrats(lbwfs_dtrat, "lbwfs_dtrat")

    # To put in the simulation parameters
    powfs_dtrat_array = [int(howfs_dtrat), int(tt_dtrat), int(lbwfs_dtrat)]

    # Calculate signal level (siglev), background (bkgrnd) and noise equivalent angle
    # parameters
    # TODO: Does this change if we have 4 vs. 1 LGS?
    num_LGS = hdr_tbl_row['num_lgs_wfs']
    
    lgs_wfs_name = 'LGSWFS-OCAM2K' # For KAPA on Keck 1 5/5/26
    tt_wfs_name = ('TRICK' if 'trick' in str(hdr_tbl_row['OSIRIS_tt_sensor']).lower() else 'STRAP')
    # TODO: powfs.wvl [lgs, tt, truth] should change based on TRICK-H vs. TRICK-K
    if tt_wfs_name == 'TRICK':
        if hdr_tbl_row['wavelength'] > 2.1 * u.micron:
            tt_wfs_name = 'TRICK-H' # We're in K band (~2.12 micron) so TRICK is taking the H-band
        else:
            tt_wfs_name = 'TRICK-K' # vice versa
    
    # Get data from lbwfs
    lbwfs_fwhm = hdr_tbl_row['lbwfs_fwhm']
    guide_star_magnitude = hdr_tbl_row['tt_gs_r_mag']
    
    r0 = hdr_tbl_row['r0']
    
    # Note, I think nearecon only matters for tip-tilt since it uses geometric 
    # not physical optics, but the signal levels & backgrounds are important still
    """
    SNR
        Not used
    
    sigma_theta : float
        Noise-equivalent angle in milliarcsec
    
    Np          : float
        Number of photons from the star per pixel within subaperture # TODO Note MAOS expects e-/subap/frame

    Nb          : float
        Number of background photons per pixel within subaperture111
    
    """
    _, howfs_nearecon, howfs_siglev, howfs_bkgrnd = mu.keck_nea_photons(guide_star_magnitude, lgs_wfs_name,
                                                                        r0, shwfs_int_time.to_value(u.s))
    _, tt_nearecon, tt_siglev, tt_bkgrnd = mu.keck_nea_photons(guide_star_magnitude, tt_wfs_name, 
                                                                        r0, tt_int_time.to_value(u.s))
    _, lbwfs_nearecon, lbwfs_siglev, lbwfs_bkgrnd = mu.keck_nea_photons(guide_star_magnitude, 'LBWFS',
                                                                        r0, lbwfs_int_time.to_value(u.s),
                                                                        lbwfs_fwhm=lbwfs_fwhm)
    # TODO: make sure nearecon & siglev are correct
    # Potentially just snag siglev from telemetry?
    nearecon = [howfs_nearecon, tt_nearecon, lbwfs_nearecon] # nearecon is in milliarcseconds
    siglev = [howfs_siglev, tt_siglev, lbwfs_siglev] # siglev is in number of photons from the star per pixel within subaperture
    
    # KAPA is ~18.6e- / pixel (bias) @ 500Hz
    # sky background in e/px/LGS frame
    # TODO: Make sure background is OK
    # Potentially just snag this from telemetry too?
    bkgrnd = [howfs_bkgrnd, tt_bkgrnd, lbwfs_bkgrnd] # bkgrnd is in number of background photons per pixel within subaperture

    # Setup final sim parameters

    # TODO: Use AOTSX/Y (mm) & KAI to convert stage offset between center of FOV
    # and tip tilt star - see kai/reduce/kai_util.py getaotsxy(hdr) or aotsxy2pix
    # which is relative between two positions but has math we want (could use x y
    # from engineering star as a base to offset from)
    # ALso nted AOLB XY OFF for offset of LBWFS on NGS (unsure of units)
    # RAOFF & DECOFF (some other kind of offset for dithering?)

    # wfs.thetax/y is where the WFS is located on the sky I guess -> offset of TT NGS in arcsec
    # So we'd change this for tip-tilt & lbwfs? or just tip tilt? does the lbwfs look at the laser spot?
    # Note wfs.thetax = [-5.5, 5.5, -5.5, 5.5] + 7'' for LGS, wfs.thetay = [-5.5, -5.5, 5.5, 5.5]
    # Yeah so we'd have wfs.thetax = [lgs1_x, lgs2_x, lgs3_x, lgs4_x, tt_x, lbwfs_x] and same for y I think
    # TODO: Put LGS ones in MAOS but leave tt & lbwfs to adjust based on AOTSX/Y

    # size = # psfsize should be ~1 arcsec in diameter (half arcsec radius). '
    # Note 25mas pixel scale set automatically by maos to get psfsize (actually it scales I think?)
    # size = ?

    # Set reconstructor parameters for LTAO
    if num_LGS == 4:
        # powfs.nwfs = [4, 1, 1] else [1, 1, 1]
        
        # What hte reconstructor thinks the atmosphere is
        # atmr.r0 = 0.15 # Avinash note
        # atmr.L0 = 30 # Avinash note
        layeralt = [0.,500,1000,2000,4000,8000,16000] # atmr.ht
        layerfrac = [0.4557,0.1295,0.0442,0.0506,0.1167,0.0926,0.1107] # atmr.wt
        # atmr.ht = layeralt # Height of layers to reconstruct
        # atmr.wt = layerfrac # Weight of layers to reconstruct
        


    # atmr.hs = ? # Height of high-order gs (LGS?) set to sodium layer height
    # powfs.hs = ? # Height of guide star -> set to sodium layer height for LGS ones & inf for tt ones

    # Final maos command
    #TODO 
    # - put sim.dt and sim.dtref in maos command. dt is fastest integration time, dtref is the dt of each wfs to use for unit scaling
    # - make sure sim.zadeg (zenith angle) is set
    # - powfs.siglev is e-/subaperture/dtref step (set based on telem?)
    # - can do powfs.siglevs (one for each WFS). Note dtref is dt so it might be different than the LGS frame rate
    # - powfs.rne is read noise in e-/pix/frame (set based on telem?)
    # - powfs.bckgrnd is background in e-/pix/lgs frame (set based on telem?)
    # - atm.wddeg from wind direction (degrees)
    # - set powfs_dtrat_array s 
    # - set num_LGS
    # - sim.end to get 2s of on-sky time
    # - set wavelength of observation evl.wvl
    # - evl.psf & psfsize
    # - make sure atm.r0z, etc are set from telemetry

    # atm.r0z from telemetry
    # atm.ht
    # atm.wt
    # Check atm.ws & wddeg set in telemetry
    # atm.ws
    # atm.wddeg = [array] # set wind directions
    # atm.wdrand = 0 if so (don't randomize wind directions)

    # recon.glao = 1 if GLAO mode (average gradient from WFS of same type)

    

def notes_from_run_maos_comp_to_sim_sky(hdr_tbl_row):

    size = 256
    mode = ''
    if simtype == 'piston':
        mode = 'piston'
        surf_cmd = ["Keck_ncpa_rmswfe130nm.fits"]
        # Fetch name of current input PSD FITS file in MAOS config file keck_sim.conf
        psd_file = ''
    elif simtype == 'psd+ncpa-seen':
        mode = 'surf_wfs1'
        # Revert NCPA surf back to Lianqi's default since tuned version is now extremely out of date/based on
        # parameters that have been significantly changed

        # surf_cmd = ["Keck_ncpa_rmswfe130nm.fits", "'r0=0.36;l0=3.39;ht=40000;slope=-2; SURFWFS=1; SURFEVL=1; seed=10;'"]
        surf_cmd = ["Keck_ncpa_rmswfe130nm.fits", "'r0=0.1;l0=10;ht=40000;slope=-2; SURFWFS=1; SURFEVL=1; seed=10;'"]
        psd_file = "PSD_Keck_ws26.47mas_vib26mas_rad2.fits"
    elif simtype == 'psd+ncpa-unseen':
        mode = 'surf_wfs0'
        # Revert NCPA surf back to Lianqi's default since tuned version is now extremely out of date/based on
        # parameters that have been significantly changed
        # surf_cmd = ["Keck_ncpa_rmswfe130nm.fits", "'r0=0.36;l0=3.39;ht=40000;slope=-2; SURFWFS=0; SURFEVL=1; seed=10;'"]
        surf_cmd = ["Keck_ncpa_rmswfe130nm.fits", "'r0=0.1;l0=10;ht=40000;slope=-2; SURFWFS=1; SURFEVL=1; seed=10;'"]
        psd_file = "PSD_Keck_ws26.47mas_vib26mas_rad2.fits"
    else:
        raise ValueError(f"Invalid MAOS simulation type '{type}'. Valid types are currently: 'piston', 'psd+ncpa-seen', 'psd+ncpa-unseen'. See help() for further info")

    for seed in seeds:
        # Must be in MAOS simulation directory to run successfully
        if os.getcwd() != baseroot.as_posix():
            os.chdir(baseroot)

        maos_cmd = f"maos -o A_keck_scao_lgs_gc_{mode}_comp_{sky[0]}_seed{seed}_epoch{sky[1]} -c A_keck_scao_lgs_gc.conf evl.psfsize={size} sim.seeds={seed} evl.wvl={wvl} powfs.dtrat={dtrat} sim.zadeg={angle} powfs.siglev={siglev} powfs.bkgrnd={bkgrnd} powfs.nearecon={nearecon} sim.wspsd={psd_file} atm.r0z={fried.value} atm.wt={turbpro} atm.ws={windspds} atm.wddeg={winddrcts} surf={surf_cmd} -O"
        os.system(maos_cmd)