# Acronyms: AO — adaptive optics; AOFCSALT — AO facility sodium-layer height (FITS keyword, metres);
# AOLBFWHM — AO LBWFS spot FWHM (FITS keyword); ASIC — application-specific integrated circuit;
# CUNIT — FITS WCS axis unit keyword prefix (e.g. CUNIT1); DEC — declination (FITS keyword in extras);
# DM — deformable mirror; DMMRFN — DM matrix / reconstructor filename (FITS keyword);
# ELAPTIME — elapsed on-sky exposure time (FITS keyword); EXPTIME — exposure time (FITS keyword);
# FITS — Flexible Image Transport System; FITS COMMENT — text after “/” on a header card (human-readable unit or note); FWHM — full width at half maximum;
# GUIDFWHM — guide-star FWHM (FITS keyword in extras); HA — hour angle (FITS keyword in extras); HDU — Header/Data Unit;
# IMAG — imager (OSIRIS imaging mode); KAPA — Keck All-sky Precision Adaptive-optics; KOA — Keck Observatory Archive;
# LGRMSWF — LGS RMS wavefront residual (FITS keyword); LGS — laser guide star;
# MAOS — Multi-Atmosphere Adaptive Optics Simulation; MJD-OBS — modified Julian date of observation (FITS keyword);
# O1FPS — high-order WFS frame rate in hertz (FITS keyword); OSIRIS — OH-Suppressing Infra-Red Imaging Spectrograph;
# POXPOS — pointing/position telemetry keyword family prefix on some engineering cards; RA — right ascension (FITS keyword in extras);
# RMS — root mean square; STRAP — Keck low-order AO sensor; STINTTIM — STRAP integration time, on-disk in ms (FITS keyword);
# TARGWAVE — target wavelength reference (FITS keyword, often carried in extras on OSIRIS); TIME-OBS — time of observation keyword (FITS, unused on raw OSIRIS here);
# UTC — Coordinated Universal Time (here: sexagesimal UT keyword on OSIRIS headers); WAVECNTR — band centre wavelength in nm (FITS keyword);
# WCS — world coordinate system; WFE — wavefront error; WFS — wavefront sensor.

"""Keck OSIRIS (raw IMAG) FITS headers → :class:`~reflectao.observation.FrameMetadata`.

Keyword mapping was checked against a KAPA commissioning raw frame
(``DATE-OBS`` + ``UTC`` time, KOA duplicate keywords, no ``TIME-OBS``).

**Full header audit:** All **517** primary-HDU cards from
``i260226_a010002.fits`` were read (keyword, value, COMMENT). Typed fields
below cover **geometry, time, band, high-order WFS rate, legacy MAOS-oriented
AO metrics** (``AOLBFWHM``, ``LGRMSWF``), **LGS sodium height** (``AOFCSALT``),
**DM reconstructor name**, and **tube temperature**. ``extra_keys`` pulls a
second tier (weather, guide star, detector timing, gains, program metadata)
without copying hundreds of engineering-only cards (ASIC voltages, POXPOS*,
power outlets, etc.).

**Units:** FITS does not attach machine-readable units to arbitrary keyword values.
When present, the text after ``/`` on the header card (the COMMENT field) often states
the unit (e.g. ``O1FPS`` comment says hertz; ``STINTTIM`` says milli-sec). There is no
guarantee every card is documented; WCS uses ``CUNIT*`` for axes, but AO telemetry
keywords are convention + comments.

**STRAP vs high-order:** ``STINTTIM`` (ms on disk) → ``ao_strap_integration_time`` (seconds).
``O1FPS`` → ``ao_wfs_frame_rate`` (hertz). Quantities use :mod:`astropy.units`.

**Exposure time:** KOA total on-sky time is ``ELAPTIME``; ``EXPTIME`` is a per-coadd
fallback — order is preserved in :attr:`exptime_keys` below.
"""

from __future__ import annotations

from reflectao.observation.fits_header import DefaultHeaderAdapter, HeaderKeywordMap

# Raw OSIRIS IMAG headers use ``UTC`` (sexagesimal time-of-day), not ``TIME-OBS``.
_OSIRIS_EXTRA_KEYS: tuple[str, ...] = (
    "TARGWAVE",
    "HA",
    "RA",
    "DEC",
    "COADDS",
    "TRUITIME",
    "READTIME",
    "ITIME",
    "SAMPMODE",
    "GUIDFWHM",
    "GUIDTIME",
    "GUIDWAVE",
    "DMGAIN",
    "O1SMGN",
    "AOOPSMOD",
    "AOAOAMED",
    "SYSGAIN",
    "PSCALE",
    "SCALE",
    "WAVEBLUE",
    "WAVERED",
    "SFILTER",
    "IFILTER",
    "IF2NAME",
    "PROGID",
    "OBJECT",
    "DATASET",
    "WXDOMHUM",
    "WXDOMTMP",
    "WXOUTHUM",
    "WXOUTTMP",
    "WXPRESS",
    "WXWNDIR",
    "WXWNDSP",
    "WXTIME",
    "AODTSTAT",
    "AODMSTAT",
    "WFDM",
    "WFTT",
    "AOFOMODE",
    "DTSENSOR",
    "UTGAIN",
    "DTGAIN",
    "AODRENA",
    "FRAMENUM",
)

OSIRIS_HEADER_KEYWORDS = HeaderKeywordMap(
    date_obs="DATE-OBS",
    time_obs="UTC",
    exptime_keys=("ELAPTIME", "EXPTIME"),
    airmass="AIRMASS",
    filter_name="FILTER",
    wavelength_nm="WAVECNTR",
    wavelength_input_unit="nm",
    ao_wfs_frame_rate_hz="O1FPS",
    ao_strap_integration_time_ms="STINTTIM",
    mjd_obs="MJD-OBS",
    telescope_name="TELESCOP",
    telescope_el_deg="EL",
    telescope_az_deg="AZ",
    parallactic_deg="PARANG",
    target_ra_deg="TARGRA",
    target_dec_deg="TARGDEC",
    tube_temp_c="TUBETEMP",
    ao_lbwfs_fwhm_arcsec="AOLBFWHM",
    lgs_rms_wf_residual="LGRMSWF",
    sodium_layer_altitude_m="AOFCSALT",
    dm_reconstructor_file="DMMRFN",
    extra_keys=_OSIRIS_EXTRA_KEYS,
)


class OsirisHeaderAdapter(DefaultHeaderAdapter):
    """Parse Keck OSIRIS raw science primary HDU headers into :class:`~reflectao.observation.FrameMetadata`.

    Uses :data:`OSIRIS_HEADER_KEYWORDS` for FITS keyword names and on-disk units.

    Examples
    --------
    In-memory header::

        from astropy.io import fits
        from reflectao.observation import OsirisHeaderAdapter

        adapter = OsirisHeaderAdapter()
        meta = adapter.extract(fits.Header(), source_path=None)

    From disk, prefer :func:`reflectao.observation.read_fits_header` with ``instrument="osiris"``.
    """

    def __init__(self) -> None:
        super().__init__(OSIRIS_HEADER_KEYWORDS)
