# ReflectAO

Comparison of on-sky adaptive optics data with digital-twin / simulation results (e.g. MAOS), with a telescope- and instrument-agnostic core. 

This package is heavily work in progress as of April 2026

This repository is separate from [PAARTI](https://github.com/MovingUniverseLab/paarti/tree/dev_brooke_digia_mar_2026), but draws significant inspiration from it, especially from that of [@Brooke Digia](https://github.com/bdigia). Scientific credit, licensing, and citation for algorithms that originate in PAARTI follow that upstream project where applicable.

This package was developed in part as an exercise in the use of generative AI coding tools for scientific software development using Cursor 3.0.12.  **Note**: code and documentation has been written by AI coding assistants and *may not be accurate* at the moment.

## Physical quantities

[FrameMetadata](reflectao/observation/frame_metadata.py) stores **Astropy** [Quantity](https://docs.astropy.org/en/stable/units/quantity.html) values with package-specific **canonical units** (e.g. wavelength in nanometres, exposure time in seconds, angles in degrees, frame rates in hertz). Instrument adapters declare the **on-disk** FITS unit (microns for NIRC2 `TARGWAVE`, milliseconds for `STINTTIM`, etc.); conversion happens when headers are read. Use `read_fits_header(path, "osiris" | "nirc2", *, hdu=0)` ([fits_header.py](reflectao/observation/fits_header.py)) to load a header from disk; see [docs/frame_metadata.md](docs/frame_metadata.md) for a field-by-field reference.

## File Structure Reference

**Deduced layout** for NIRC2 files

```text
{skyroot}/{YYYYMMDD}nirc2_kp/{cccc}_psf.fits
```

Here `{YYYYMMDD}` and `nirc2_kp` are **one directory name** (no slash between date and `nirc2_kp`). `{cccc}` is the frame id (e.g. `c2061`)

**Integration test:** by default, `tests/test_nirc2_header.py` uses the lab example  
`/u/bdigia/work/ao/airopa_input/20170823nirc2_kp/c2061_psf.fits` when that path exists (skipped on machines without it). Override with `REFLECTAO_NIRC2_PSF` pointing at any other `*_psf.fits` with the usual NIRC2 header cards.

## Development

Keep personal or team planning notes under `planning/` at the repo root; that directory is listed in `.gitignore` so it is never committed.

Use your existing **conda/mamba** environment. Install **astropy**, **pytest**, and **numpy** there if needed. From this directory:

```bash
pytest
```

Use `**python -m pytest**` (your conda `python`) so Astropy and pytest match the environment.

The `reflectao` package is plain Python on the repo root (no editable `pip install` required for local work). Packaging for PyPI can be added later.

## Status

Observation metadata layer: `reflectao.observation` reads FITS headers into simulation-oriented `FrameMetadata` (`Quantity`-based types; see `tests/` and `docs/frame_metadata.md`).