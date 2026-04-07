"""Keck NIRC2 FITS headers for MAOS-style workflows.

Keyword set matches the primary-header cards used in the reference
``run_maos_comp_to_sky_sim`` logic in the archived PAARTI
``paarti/utils/maos_utils.py`` (``AIRMASS``, ``TARGWAVE`` in microns, ``STINTTIM``,
``WSFRRT``, optional ``AOLBFWHM``). Wavelength on disk is **microns**; it is read as
such and stored canonically as nanometres in :class:`~reflectao.observation.FrameMetadata`.
"""

from __future__ import annotations

from reflectao.observation.fits_header import DefaultHeaderAdapter, HeaderKeywordMap

_NIRC2_EXTRA_KEYS: tuple[str, ...] = (
    "EXPSTART",
    "EXPSTOP",
    "INSTRUME",
)

NIRC2_HEADER_KEYWORDS = HeaderKeywordMap(
    date_obs="DATE-OBS",
    time_obs="TIME-OBS",
    exptime_keys=("EXPTIME", "EXPOSURE"),
    airmass="AIRMASS",
    filter_name="FILTER",
    wavelength_nm="TARGWAVE",
    wavelength_input_unit="micron",
    ao_wfs_frame_rate_hz="WSFRRT",
    ao_strap_integration_time_ms="STINTTIM",
    strap_integration_input_unit="ms",
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
    extra_keys=_NIRC2_EXTRA_KEYS,
)


class Nirc2HeaderAdapter(DefaultHeaderAdapter):
    """Parse Keck NIRC2 science or PSF-style primary HDU headers into :class:`~reflectao.observation.FrameMetadata`.

    Uses :data:`NIRC2_HEADER_KEYWORDS` for FITS keyword names and on-disk units.

    Examples
    --------
    From disk::

        from reflectao.observation import read_fits_header

        meta = read_fits_header("path.fits", "nirc2")
    """

    def __init__(self) -> None:
        super().__init__(NIRC2_HEADER_KEYWORDS)
