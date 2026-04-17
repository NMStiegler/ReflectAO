

def notes_from_run_maos_comp_to_sim_sky:
    # Set sim_dt to fastest WFS readout rate in Hz
    if hdr_shwfs_frame_rate > hdr_strap_frame_rate:
        max_frame_rate = hdr_shwfs_frame_rate
    else:
        max_frame_rate = hdr_strap_frame_rate

    sim_dt = 1.0/max_frame_rate

    howfs_dtrat = (sim_dt / hdr_shwfs_int_time)*1000.0
    strap_dtrat = (sim_dt / hdr_strap_int_time)*1000.0

    # Compose powfs.dtrat array for input into MAOS config command override
    if howfs_dtrat < 1:
        howfs_dtrat = 1
    elif strap_dtrat < 1:
        strap_dtrat = 1
    dtrat = [int(howfs_dtrat), int(strap_dtrat), 7080]

    # Calculate siglev/bkgrnd/nearecon config parameters using variable
    # integration times from headers
    _, howfs_nearecon, howfs_siglev, howfs_bkgrnd = keck_nea_photons(8.1, 'LGSWFS', 
                                                                        hdr_shwfs_int_time/1000.0)
    _, strap_nearecon, strap_siglev, strap_bkgrnd = keck_nea_photons(14.0, 'STRAP', 
                                                                        hdr_strap_int_time/1000.0)
    _, lbwfs_nearecon, lbwfs_siglev, lbwfs_bkgrnd = keck_nea_photons(14.0, 'LBWFS', 15.0) # LBWFS integration time assumed to be 15 seconds

    # If AO LBWFS average spot size (average FWHM) keyword present in header, use this for spot size
    # and redo keck_nea_photons associated calculations for LBWFS
    if 'AOLBFWHM' in hdr:
        band = "R"
        side = 0.563 
        ps = 0.148
        sigma_e = 7.96
        theta_beta = float(hdr['AOLBFWHM'])
        # Convert header spot size to radians
        theta_beta *= ( math.pi/180.0 ) / ( 60.0*60.0 )
        throughput = 0.03
        pix_per_ap = 16.7*16.7
        m = 14
        time = 15
        # Calculate number of photons and background photons
        LBWFS_Np, LBWFS_Nb = n_photons(side, time, m, band, ps, throughput)
        # Calculate SNR
        SNR = LBWFS_Np / np.sqrt(LBWFS_Np + pix_per_ap*LBWFS_Nb + pix_per_ap*sigma_e**2)
        # Noise equivalent angle in milliarcseconds (eq 65)
        LBWFS_sigma_theta = theta_beta/SNR  * ( 180.0/math.pi ) * 60.0 * 60.0 * 1000.0

        nearecon = [howfs_nearecon, strap_nearecon, LBWFS_sigma_theta]
        siglev = [howfs_siglev, strap_siglev, LBWFS_Np]
        bkgrnd = [howfs_bkgrnd, strap_bkgrnd, LBWFS_Nb]
    else:
        nearecon = [howfs_nearecon, strap_nearecon, lbwfs_nearecon]
        siglev = [howfs_siglev, strap_siglev, lbwfs_siglev]
        bkgrnd = [howfs_bkgrnd, strap_bkgrnd, lbwfs_bkgrnd]


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