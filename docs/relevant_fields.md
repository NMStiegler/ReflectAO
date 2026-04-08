# Per-frame observation fields (logical schema)

This note lists the **universal internal quantities** we want to attach to each exposure (one frame = one row in a future table). Names here are descriptive only; they are not tied to any particular file format or implementation.

Every parameter has its **own column** in the internal table. If a source does not provide a value, that cell is empty (None or null); there is no separate catch-all bucket for “everything else.”

We try to reduce redundancy as much as possible. For example, we shouldn't store more than one representation of the same parameter (no alternates/backups).

Units in the table are the **canonical internal** choices; whatever format the data came from must be converted when building the table.

| Field | Units | Meaning |
| ----- | ----- | ------- |
| Source location | (path string, if applicable) | Filesystem path or URI of the data product this row was derived from, for provenance. |
| Program identifier | (string) | Observing program or proposal ID, when logged. |
| Object name | (string) | Target or field name string from the header, when present. |
| Dataset identifier | (string) | Dataset or collection label for the exposure, when present. |
| Instrument name | (string) | Name of the instrument (distinct from the telescope or facility identifier). |
| Telescope name | (string) | Identifier for the telescope (e.g. facility name). |
| Exposure start time | UTC timestamp | Best-effort start of the on-sky integration from the primary date and time keywords (calendar date and time-of-day combined). Interpret as UTC when the source uses a universal-time convention; some sites may use local civil time—verify for your data. |
| Exposure time | s | Total on-sky integration time for the frame (not readout overhead only). |
| Coadd count | dimensionless (integer) | Number of coadds or equivalent stacking count, when logged. |
| True elapsed time | s | On-sky or clock-related elapsed time when a distinct keyword supplies it (e.g. true shutter-open duration vs nominal). |
| Readout time | s | Readout duration associated with the exposure, when logged. |
| Detector integration time | s | Integration time from a detector-timing keyword when it is not redundant with the primary exposure time; verify against your instrument manual. |
| Sampling mode | (string or integer) | Detector sampling or read mode code or name. |
| Wavelength | nm | Effective or band-centre wavelength of the observation (primary band-centre keyword for that pipeline). |
| Target wavelength (supplementary reference) | nm | Additional target or band wavelength from a second header value when both a band centre and a target wavelength are recorded. |
| Band short wavelength | nm | Blue or short-wavelength edge of the bandpass, when logged. |
| Band long wavelength | nm | Red or long-wavelength edge of the bandpass, when logged. |
| Filter name | (string) | Primary filter name string from the header. |
| Airmass | dimensionless | Atmospheric airmass along the line of sight. |
| Target right ascension | deg | Target RA as stored with the observation (instrument/pipeline convention; not guaranteed to match a specific astrometric system without checking the data documentation). |
| Target declination | deg | Target Dec; same convention caveat as right ascension. |
| Hour angle | deg | Hour angle at observation (convert from hours to degrees at ingest if the source gives hours: multiply by 15). |
| Telescope elevation | deg | Pointing elevation angle. |
| Telescope azimuth | deg | Pointing azimuth angle. |
| Parallactic angle | deg | Parallactic angle at the time of observation. |
| High-order WFS frame rate | Hz | Sampling rate of the high-order wavefront sensor loop. |
| Tip-Tilt frame rate | Hz | Frame rate for the low-order Tip Tilt (TT) sensor (STRAP or TRICK on Keck); source data can store this in milliseconds per exposure and should be converted here. |
| LBWFS spot FWHM | arcsec | Full width at half maximum of the low-bandwidth WFS spot, when logged. |
| LGS RMS wavefront residual | nm | Telemetry-style RMS metric associated with laser guide star performance (exact optical definition depends on the pipeline; treat the number as the logged residual in nanometres). |
| Sodium / LGS layer altitude | m | Altitude context for the sodium layer or LGS path (e.g. facility model height), when present. |
| Guide star FWHM | arcsec | FWHM of the guide star as logged for the exposure. |
| Guide star integration time | s | Integration time used on the guide sensor, when logged. |
| Guide wavelength | nm | Effective wavelength of the tip-tilt / acquisition guider band (NSG channel on LGS-AO nights). |
| DM reconstructor reference | (string) | Filename or path fragment identifying the deformable-mirror reconstructor or control matrix in use, when logged. |
| Deformable mirror gain | dimensionless | DM loop gain or equivalent logged gain parameter. |
| High-order WFS gain | dimensionless | Gain on the high-order WFS path (e.g. secondary gain relative to DM gain), when logged. |
| System gain | dimensionless | Overall or system-level gain factor, when logged. |
| Up tip tilt gain | dimensionless | UT gain or equivalent keyword, when present. |
| Down tip tilt gain | dimensionless | DT gain or equivalent gain keyword, when present. |
| AO operations mode | (string or integer) | AO operations or control mode code/name. |
| AO amplitude median | (as documented) | Logged median amplitude statistic from the AO system; units follow the data provider (often dimensionless counts or normalized). |
| AO focus mode | (string or integer) | Focus or AO focus-state mode, when logged. |
| AO deformable mirror / telescope status | (string) | Status string for DM/telescope subsystem when logged. |
| AO DM status | (string) | Deformable-mirror-specific status string, when logged. |
| AO reconstructor enable / name fragment | (string) | Keyword indicating reconstructor enable state or related name fragment, when logged. |
| Wavefront DM component metric | (pipeline-specific) | Logged wavefront or error metric associated with the DM path; define units from header comments or documentation. |
| Wavefront tip-tilt component metric | (pipeline-specific) | Logged wavefront or error metric associated with tip-tilt; define units from header comments or documentation. |
| Detector sensor identifier | (string) | Sensor or detector identifier string, when logged. |
| Plate scale | arcsec / pixel | Plate scale when logged as an angular size per pixel (verify linear vs spectral axis usage from documentation). |
| Dome humidity | % relative humidity | Humidity inside the dome or coude area, when logged. |
| Dome air temperature | °C | Air temperature inside the dome (or equivalent indoor environment), when logged. |
| Outside humidity | % relative humidity | Ambient humidity outside, when logged. |
| Outside air temperature | °C | Ambient air temperature outside, when logged. |
| Barometric pressure | hPa | Atmospheric pressure at the site (or as defined by the weather keyword set), when logged. |
| Wind direction | deg | Wind direction (TODO figure out and document whether azimuth convention is meteorological or astronomical). |
| Wind speed | m/s | Wind speed (convert from other units at ingest if necessary). |
| Weather sample time | (timestamp string) | Time stamp associated with the weather block, when logged. |
| Tube temperature | °C | Temperature inside the telescope tube or equivalent relevant bulk temperature for the observation. |
| Frame number | dimensionless (integer) | Running frame or sequence number within a dataset, when logged. |

## Notes

- Pixel data are intentionally omitted; this schema is for metadata and telemetry needed to reproduce or simulate the observation context.
- The column set is the union of everything we currently preserve from the supported Keck-style loaders: there is no separate “extra” map—each listed quantity is either filled from the source or left empty.
