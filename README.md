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

Use your existing **conda/mamba** environment. From this directory, install the package in editable mode so `import reflectao` works from any working directory (not only when the repo root is on `PYTHONPATH`):

```bash
pip install -e ".[dev]"
```

That installs runtime dependencies (`numpy`, `astropy`) and the optional `dev` extra (`pytest`). If you only need to run tests from the repo root without a global install, you can instead rely on `pytest.ini` (`pythonpath = .`) and install **astropy**, **pytest**, and **numpy** in the environment yourself.

Then:

```bash
python -m pytest
```

Use `python -m pytest` (your conda `python`) so Astropy and pytest match the environment.

**Imports:** `reflectao/__init__.py` re-exports `build_observation_table` and schema helpers, so `from reflectao import build_observation_table` and `reflectao.build_observation_table` refer to the same function. Submodules such as `reflectao.kapa_utils` are imported normally, for example `import reflectao.kapa_utils as ku`.