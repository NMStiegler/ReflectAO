# ReflectAO

Comparison of on-sky adaptive optics data with digital-twin / simulation results (e.g. MAOS), with a telescope- and instrument-agnostic core.

This repository is separate from the archived **PAARTI** reference tree (`brooke_paarti/`); it will grow incrementally. Scientific credit, licensing, and citation for algorithms that originate in PAARTI follow that upstream project where applicable.

Keep personal or team planning notes under `planning/` at the repo root; that directory is listed in `.gitignore` so it is never committed.

## Physical quantities

[`FrameMetadata`](reflectao/observation/frame_metadata.py) stores **Astropy** [`Quantity`](https://docs.astropy.org/en/stable/units/quantity.html) values with **canonical units** (e.g. wavelength in nanometres, exposure time in seconds, angles in degrees, frame rates in hertz). Instrument adapters declare the **on-disk** FITS unit (microns for NIRC2 `TARGWAVE`, milliseconds for `STINTTIM`, etc.); conversion happens when headers are read. Use `read_fits_header(path, "osiris" | "nirc2", *, hdu=0)` ([`fits_header.py`](reflectao/observation/fits_header.py)) to load a header from disk; see [`docs/frame_metadata.md`](docs/frame_metadata.md) for a field-by-field reference.

## Reference NIRC2 layout (PAARTI material)

There is **no public URL** for these FITS files; PAARTI hard-codes **machine-local** roots and builds paths the same way in code and notebooks.

**Deduced layout** (see `fetch_sky_frames` / `run_maos_comp_to_sky_sim` in `brooke_paarti/paarti/paarti/utils/maos_utils.py`): a root directory `skyroot` (notebooks often call it `onskyroot` or `skyroot`, e.g. `.../airopa_input/`) plus:

```text
{skyroot}/{YYYYMMDD}nirc2_kp/{cccc}_psf.fits
```

Here `{YYYYMMDD}` and `nirc2_kp` are **one directory name** (no slash between date and `nirc2_kp`). `{cccc}` is the frame id (e.g. `c2061`). Example full paths from the reference notebooks (original author’s workstation):

- `.../airopa_input/20170823nirc2_kp/c2061_psf.fits` — `keck_nea_photons_validation.ipynb`
- `.../airopa_input/20170823nirc2_kp/c2058_psf.fits` (and `c2059`…`c2062`) — `run_maos_telemetry.ipynb`
- `.../airopa_input/20100815nirc2_kp/c0120_psf.fits` — `run_maos.ipynb`

**Integration test:** by default, `tests/test_nirc2_header.py` uses the lab example  
`/u/bdigia/work/ao/airopa_input/20170823nirc2_kp/c2061_psf.fits` when that path exists (skipped on machines without it). Override with `REFLECTAO_NIRC2_PSF` pointing at any other `*_psf.fits` with the usual MAOS header cards.

## Development

Use your existing **conda/mamba** environment. Install **astropy**, **pytest**, and **numpy** there if needed. From this directory:

```bash
pytest
```

Use **`python -m pytest`** (your conda `python`) so Astropy and pytest match the environment.

The `reflectao` package is plain Python on the repo root (no editable `pip install` required for local work). Packaging for PyPI can be added later.

## Status

Observation metadata layer: `reflectao.observation` reads FITS headers into simulation-oriented `FrameMetadata` (`Quantity`-based types; see `tests/` and `docs/frame_metadata.md`).
