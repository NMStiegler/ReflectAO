# FrameMetadata field reference

This document describes the **logical** (in-memory) meaning of each field on
[`FrameMetadata`](../reflectao/observation/frame_metadata.py). It is the
authoritative field-by-field reference; the dataclass docstring duplicates the
same information in a form suitable for Sphinx.

Values are produced by [`DefaultHeaderAdapter.extract`](../reflectao/observation/fits_header.py)
(or instrument subclasses) from a FITS header. Which FITS keywords populate each
field is defined per instrument in [`HeaderKeywordMap`](../reflectao/observation/fits_header.py)
instances (e.g. [`OSIRIS_HEADER_KEYWORDS`](../reflectao/observation/instruments/osiris.py),
[`NIRC2_HEADER_KEYWORDS`](../reflectao/observation/instruments/nirc2.py)).

| Field | Type | Canonical unit / notes | Meaning |
| ----- | ---- | ---------------------- | ------- |
| `source_path` | `pathlib.Path \| None` | N/A | Filesystem path of the FITS file when metadata was read from disk (`read_fits_header` / `extract(..., source_path=...)`). |
| `exposure_start` | `datetime.datetime \| None` | UTC (timezone-aware) | Best-effort start time of the exposure from date + time header keywords. Parsing uses [`astropy.time.Time`](https://docs.astropy.org/en/stable/time/index.html); invalid combinations yield `None`. **Not sure** whether all instruments use UTC vs local civil time on the mountain; treat as UTC when ISO-like strings parse as such. |
| `exposure_time` | `astropy.units.Quantity \| None` | second | On-sky integration time. Multiple FITS keywords may be candidates; resolution uses [`KeywordResolutionPolicy.FIRST_NON_NULL`](../reflectao/observation/fits_header.py) on `HeaderKeywordMap.exptime_keys` (e.g. OSIRIS: `ELAPTIME` then `EXPTIME`). |
| `wavelength` | `astropy.units.Quantity \| None` | nanometre | Effective or band-centre wavelength after converting from the instrument’s on-disk unit (`wavelength_input_unit` in the keyword map, e.g. microns for NIRC2 `TARGWAVE`). |
| `filter_name` | `str \| None` | N/A | Filter name string from the header (e.g. `FILTER`). |
| `airmass` | `astropy.units.Quantity \| None` | dimensionless | Airmass as a dimensionless quantity. |
| `ao_wfs_frame_rate` | `astropy.units.Quantity \| None` | hertz | High-order wavefront sensor sampling rate (e.g. OSIRIS `O1FPS`, NIRC2 `WSFRRT`). |
| `ao_strap_integration_time` | `astropy.units.Quantity \| None` | second | STRAP (low-order) integration time; on-disk values are often milliseconds (`STINTTIM`-like keywords). |
| `mjd_obs` | `float \| None` | dimensionless (MJD) | Modified Julian Date when the header provides it (`MJD-OBS` or instrument equivalent). |
| `telescope_name` | `str \| None` | N/A | Telescope identifier (e.g. `TELESCOP`). |
| `telescope_elevation` | `astropy.units.Quantity \| None` | degree | Telescope elevation angle. |
| `telescope_azimuth` | `astropy.units.Quantity \| None` | degree | Telescope azimuth angle. |
| `parallactic_angle` | `astropy.units.Quantity \| None` | degree | Parallactic angle. |
| `target_ra` | `astropy.units.Quantity \| None` | degree | Target right ascension as stored in the header (instrument-specific convention; **not** guaranteed to be ICRS or precessed—verify against your WCS if sub-arcsecond astrometry matters). |
| `target_dec` | `astropy.units.Quantity \| None` | degree | Target declination; same caveat as `target_ra`. |
| `tube_temperature` | `astropy.units.Quantity \| None` | degree Celsius | Tube or relevant telescope temperature; on-disk values (e.g. `TUBETEMP`) are °C and stored without offset. |
| `ao_lbwfs_fwhm` | `astropy.units.Quantity \| None` | arcsecond | LBWFS spot FWHM (`AOLBFWHM` where present). |
| `lgs_rms_wf_residual` | `astropy.units.Quantity \| None` | nanometre | Telemetry-style RMS metric (e.g. `LGRMSWF`). **Semantic uncertainty:** the exact optical definition (WFE vs slope, single-mode vs combined) depends on Keck/AO pipeline documentation; this package only parses the numeric card into nanometres. Asserted in [`test_osiris_header.py`](../tests/test_osiris_header.py) and [`test_nirc2_header.py`](../tests/test_nirc2_header.py) against synthetic headers. |
| `sodium_layer_altitude` | `astropy.units.Quantity \| None` | metre | Sodium / LGS layer altitude context (e.g. OSIRIS `AOFCSALT`). |
| `dm_reconstructor_file` | `str \| None` | N/A | Filename or path fragment for the DM reconstructor / control matrix when logged (e.g. `DMMRFN`). |
| `extra` | `Mapping[str, Any]` | N/A | Raw header values for keys listed in `HeaderKeywordMap.extra_keys` for that instrument. Used for secondary cards (weather, guide star FWHM, etc.) without promoting every engineering keyword to a typed field. |

## Related API

- **Read from disk:** [`read_fits_header(path, instrument, *, hdu=...)`](../reflectao/observation/fits_header.py) with `instrument` in `{"osiris", "nirc2"}`.
- **In-memory header:** instantiate `OsirisHeaderAdapter` or `Nirc2HeaderAdapter` and call `extract(header, source_path=...)`.
- **Custom instruments:** subclass `DefaultHeaderAdapter` (import from `reflectao.observation.fits_header`) and supply a `HeaderKeywordMap`.
