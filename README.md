# ReflectAO

Comparison of on-sky adaptive optics data with digital-twin / simulation results (e.g. MAOS), with a telescope- and instrument-agnostic core. 

This package is in a minimum viable product state as of May 17th, 2026.

This repository is separate from [PAARTI](https://github.com/MovingUniverseLab/paarti/tree/dev_brooke_digia_mar_2026), but draws significant inspiration from it, especially from the work of [@Brooke Digia](https://github.com/bdigia). Scientific credit, licensing, and citation for algorithms that originate in PAARTI follow that upstream project where applicable.

This earlier versions of this package was developed in part as an exercise in the use of generative AI coding tools for scientific software development using Cursor 3.0.12.  **Note**: code and documentation has been written by AI coding assistants from early commits *may not be accurate*. As of May 2026, generative AI coding tools have been phased out of the project and most of the code has been extensively reviewed by eye. As such, errors are likely now more human in origin than hallucinatory.

## Installation and Dependencies
**Installation**

As this package is WIP, it's not on any package repositories. I'll describe here the method we've been using to install the package locally. There may be a better way to do this. 

To install the package, clone this repository to your machine using `git clone git@github.com:NMStiegler/ReflectAO.git` and then add the location of the `reflectao` subdirectory to your python path in your `bashrc`, `zshenv`, or similar shell environment/setup script like: `export PYTHONPATH="$PYTHONPATH:/location/of/ReflectAO/"`. For the latest features, make sure to `git switch dev` (or whichever branch makes sense). Make sure to `rebash` or `source ~/.zshenv` to make sure the changes update before trying to use the package. 

**Dependencies**

This package relies on several other packages developed by the [Moving Universe Lab](https://github.com/MovingUniverseLab) at Berkeley
- `KAI` (potentially a dev branch, currently `dev-for-ao-sims`) the Keck Adaptive optics Imaging reduction pipeline
- `paarti` (potentially a dev branch, currently `dev`)
- maybe others / other packages not from MULab

See `environment.yml` for the Python packages we've been using to develop this package. You can recreate the environment using `conda env create -f environment.yml`. 

## File Structure Reference

**Package Files**

`reflectao`: The python package source files
- `build_observation_table.py`: Contains the main function `build_observation_table` which collects information about an on-sky exposure taken with adaptive optics. It pulls information from the FITS file header, OCAM2K LGSWFS telemetry, LBWFS telemetry, MASS/DIMM measurements, and weather data from various sources around the Keck observatory site.
- `kapa_utils.py`: Contains utility functions specific to the KAPA LGSAO system on Keck I, such as conversions for the OCAM2K LGS WFS detector, as well as collated metadata about our KAPA commissioning, like a list of dates from which we have data
- `maos_utils.py`: Contains utility functions for working with [MAOS](https://lianqiw.github.io/maos/index.html), the Multi-threaded Adaptive Optics Simulator
- `run_sim.py`: Contains functions for running simulations using various simulation suites. Currently only configured for MAOS with the run_maos_sim command. Gets the information for how to run the simulation from the observation table
- `schema.py`: Contains the schema for how the `observation_table` is built, including units and comments. The table may contain more columns, but this is a starting place. It's meant to define the basic properties that a row should have. It should match all the columns assigned in `build_observation_table` but may not since it hasn't closely been checked
- `telemetry_utils.py`: Contains utility functions for working with data from the Keck I telemetry system as well as utility functions for locating files on the MULab group cluster (at some point these should probably be separated). Also contains collated metadata / notes about the nights, for example which ones don't have telemetry.

`scripts`: Contains some scripts used in development of the package
- `dump_fits_header.py`: Saves a FITS header to a file for inspection
- `find_frames_with_darks_LTAO`: Tries to find OCAM2K telemetry data from while darks were being taken in order to characterize read noise and dark current. Method wasn't super successful
- `find_frames_with_good_LTAO_telemetry`: Prints a dictionary of (night->set->image) trees to sort images which were taken with closed-loops and on a target, etc from darks or images with open loop operation (which aren't used as science images)
- `get_all_ocam2k_read_stats.py`: Goes through all the images from `find_frames_with_good_LTAO_telemetry` (stored in `telemetry_utils.py`), analyzes data from their ocam2k telemetry, and saves statistics about the signal and background levels in subapertures and pixels to a file
- `run_first_sim.py`: Script which runs a MAOS sim on an example image file

`tests`: Unit tests of the package. Not very complete at the moment

**Documentation**

All files should have documentation in them and in each function in the sphynx/doxygen format. The goal was to be able to generate documentation from these comments easily, although this hasn't yet been done. Reach out to the authors with any specific questions.

**OSIRIS Files**

An example file can be found at `/g3/data/kapa/2026feb26/raw/i260226_a010002.fits` on certain machines. Use it for testing & learning what is in a KAPA OSIRIS .fits file if it exists on the machine you're on (MULab cluster).

This package also takes input from OSIRIS / Keck I telemetry, which on the MULab cluster is located at, `/g3/data/kapa/NIGHT/telemetry/IMAG/IMAGE_FILENAME.fits/...`, `/g3/data/kapa/2026feb26/telemetry/IMAG/i260226_a010002.fits/` for example.

## Development

Use the conda environment from `environment.yml`.

**Testing**
Testing is done with pytest which defines unit tests for functions in the package. Testing has lagged behind production and so many of the functions are not well tested if they are tested at all.

Run `python -m pytest` to run the unit testing suite.

**Notes**
Keep personal or team planning notes under `planning/` at the repo root; that directory is listed in `.gitignore` so it is never committed.