"""
Read a finished MAOS simulation directory and print the Strehl, empirical FWHM,
Gaussian FWHM, and encircled-energy radii at each field position and wavelength.

Mirrors the analysis in MovingUniverseLab/keck_maos kapa/results.ipynb.
Uses paarti.psf_metrics.metrics.calc_psf_metrics_single for all measurements —
the same empirical approach used for on-sky PSF characterisation.

Usage:
    python print_maos_strehl.py <sim_directory> [--seed SEED] [--oversamp N]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from paarti.psf_metrics import metrics


def _position_from_header(hdr) -> tuple[float, float]:
    """Return (x_arcsec, y_arcsec) from a MAOS PSF HDU header.

    MAOS stores the evaluation direction as FITS complex keyword THETA where
    real part = y offset and imaginary part = x offset, both in arcseconds
    (matching paarti.psfs.MAOS_PSF_stack convention).
    """
    theta = hdr.get("THETA")
    if theta is None:
        raise KeyError("THETA")
    if isinstance(theta, complex):
        return float(theta.imag), float(theta.real)
    # Astropy sometimes returns a '(re,im)' string for FITS complex keywords.
    m = re.match(r"\(\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*\)", str(theta))
    if m:
        return float(m.group(2)), float(m.group(1))  # imag=x, real=y
    raise ValueError(f"Cannot parse THETA keyword: {theta!r}")


def _position_from_filename(path: Path) -> tuple[float, float]:
    """Parse (x, y) arcsec offsets from filename 'evlpsfcl_N_xX_yY.fits'."""
    m = re.search(r"_x(-?[\d]+)_y(-?[\d]+)\.fits$", path.name)
    if not m:
        return 0.0, 0.0
    return float(m.group(1)), float(m.group(2))


def get_psf_metrics_over_field(
    directory: Path, seed: int = 1, oversamp: int = 3
) -> tuple[np.ndarray, ...]:
    """Compute PSF metrics at every field position and wavelength.

    Equivalent to maos_utils.get_psf_metrics_over_field from the keck_maos
    results notebook, implemented via paarti.psf_metrics.metrics.

    Returns
    -------
    x, y   : (n_pos, n_wvl) field positions in arcsec
    w      : (n_pos, n_wvl) wavelengths in microns
    strehl : (n_pos, n_wvl) Strehl ratio (peak pixel — MAOS normalisation)
    fwhm_g : (n_pos, n_wvl) Gaussian-fit FWHM in arcsec
    fwhm_e : (n_pos, n_wvl) empirical FWHM in arcsec (area of pixels > 0.5 * peak)
    r_ee50 : (n_pos, n_wvl) 50% encircled-energy radius in arcsec
    r_ee80 : (n_pos, n_wvl) 80% encircled-energy radius in arcsec
    """
    fits_files = sorted(directory.glob(f"evlpsfcl_{seed}_x*_y*.fits"))
    if not fits_files:
        print(
            f"Error: no files matching 'evlpsfcl_{seed}_x*_y*.fits' in {directory}",
            file=sys.stderr,
        )
        sys.exit(1)

    with fits.open(fits_files[0]) as h:
        n_wvl = len(h)
    n_pos = len(fits_files)

    x      = np.zeros((n_pos, n_wvl))
    y      = np.zeros((n_pos, n_wvl))
    w      = np.zeros((n_pos, n_wvl))
    strehl = np.zeros((n_pos, n_wvl))
    fwhm_g = np.zeros((n_pos, n_wvl))
    fwhm_e = np.zeros((n_pos, n_wvl))
    r_ee50 = np.zeros((n_pos, n_wvl))
    r_ee80 = np.zeros((n_pos, n_wvl))

    for i, fpath in enumerate(fits_files):
        with fits.open(fpath) as hdul:
            try:
                xi, yi = _position_from_header(hdul[0].header)
            except (KeyError, ValueError):
                xi, yi = _position_from_filename(fpath)

            for j, hdu in enumerate(hdul):
                hdr = hdu.header
                psf = hdu.data

                mets = metrics.calc_psf_metrics_single(
                    psf, hdr["DP"], oversamp=oversamp
                )

                x[i, j] = xi
                y[i, j] = yi
                w[i, j] = hdr["WVL"] * 1e6        # metres -> microns
                strehl[i, j] = mets["strehl"]
                fwhm_g[i, j] = mets["fwhm"]        # Gaussian-fit FWHM (arcsec)
                fwhm_e[i, j] = mets["emp_fwhm"]    # empirical FWHM (arcsec)
                r_ee50[i, j] = mets["ee50"]         # 50% EE radius (arcsec)
                r_ee80[i, j] = mets["ee80"]         # 80% EE radius (arcsec)

    return x, y, w, strehl, fwhm_g, fwhm_e, r_ee50, r_ee80


def print_results(
    x: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    strehl: np.ndarray,
    fwhm_g: np.ndarray,
    fwhm_e: np.ndarray,
    r_ee50: np.ndarray,
    r_ee80: np.ndarray,
) -> None:
    n_pos, n_wvl = strehl.shape
    wl_nm = w[0, :] * 1e3  # microns -> nm

    r = np.sqrt(x[:, 0] ** 2 + y[:, 0] ** 2)
    order = np.argsort(r)

    header = (
        f"{'X (arcsec)':>12}  {'Y (arcsec)':>12}  {'R (arcsec)':>12}"
        f"  {'Strehl':>8}"
        f"  {'FWHM_g (mas)':>13}"
        f"  {'FWHM_e (mas)':>13}"
        f"  {'EE50 (mas)':>11}"
        f"  {'EE80 (mas)':>11}"
    )
    sep = "-" * len(header)

    for j in range(n_wvl):
        print(f"\nWavelength: {wl_nm[j]:.0f} nm")
        print(header)
        print(sep)
        for i in order:
            print(
                f"{x[i, j]:>12.4f}  {y[i, j]:>12.4f}  {r[i]:>12.4f}"
                f"  {strehl[i, j]:>8.4f}"
                f"  {fwhm_g[i, j] * 1e3:>13.2f}"
                f"  {fwhm_e[i, j] * 1e3:>13.2f}"
                f"  {r_ee50[i, j] * 1e3:>11.2f}"
                f"  {r_ee80[i, j] * 1e3:>11.2f}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Print MAOS simulation Strehl and PSF metrics at each field position "
            "and wavelength using paarti's empirical measurement (same as on-sky). "
            "Reads closed-loop PSF files (evlpsfcl_*) from a finished simulation directory."
        )
    )
    parser.add_argument(
        "directory",
        help="Path to the finished MAOS simulation results directory",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Simulation seed used when MAOS was run (default: 1)",
    )
    parser.add_argument(
        "--oversamp",
        type=int,
        default=3,
        help="PSF oversampling factor passed to paarti (default: 3)",
    )
    args = parser.parse_args()

    sim_dir = Path(args.directory).expanduser().resolve()
    if not sim_dir.is_dir():
        print(f"Error: directory not found: {sim_dir}", file=sys.stderr)
        return 1

    x, y, w, strehl, fwhm_g, fwhm_e, r_ee50, r_ee80 = get_psf_metrics_over_field(
        sim_dir, seed=args.seed, oversamp=args.oversamp
    )
    print_results(x, y, w, strehl, fwhm_g, fwhm_e, r_ee50, r_ee80)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
