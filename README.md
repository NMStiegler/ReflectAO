# ReflectAO

Comparison of on-sky adaptive optics data with digital-twin / simulation results (e.g. MAOS), with a telescope- and instrument-agnostic core. 

This package is heavily work in progress as of April 2026

This repository is separate from [PAARTI](https://github.com/MovingUniverseLab/paarti/tree/dev_brooke_digia_mar_2026), but draws significant inspiration from it, especially from that of [@Brooke Digia](https://github.com/bdigia). Scientific credit, licensing, and citation for algorithms that originate in PAARTI follow that upstream project where applicable.

This package was developed in part as an exercise in the use of generative AI coding tools for scientific software development using Cursor 3.0.12.  **Note**: code and documentation has been written by AI coding assistants and *may not be accurate* at the moment.

## File Structure Reference

**OSIRIS Files**

An example file can be found at `/g3/data/kapa/2026feb26/raw/i260226_a010002.fits` on certain machines. Use it for testing & learning what is in a KAPA OSIRIS .fits file if it exists on the machine you're on (MULab cluster).

## Development

Keep personal or team planning notes under `planning/` at the repo root; that directory is listed in `.gitignore` so it is never committed.

Use your existing **conda/mamba** environment. Install **astropy**, **pytest**, and **numpy** there if needed. From this directory:

```bash
pytest
```

Use `**python -m pytest**` (your conda `python`) so Astropy and pytest match the environment.

The `reflectao` package is plain Python on the repo root (no editable `pip install` required for local work). Packaging for PyPI can be added later.