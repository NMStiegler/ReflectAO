"""
Read a finished MAOS simulation directory and print Strehl, empirical FWHM,
Gaussian FWHM, and encircled-energy radii at each field position and wavelength.

Mirrors the analysis in MovingUniverseLab/keck_maos kapa/results.ipynb by calling
paarti.utils.maos_utils.get_psf_metrics_over_field directly (requires paarti dev branch).

Usage:
    python print_maos_strehl.py <sim_directory> [--seed SEED] [--oversamp N]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from paarti.utils import maos_utils


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
    # w is in microns; fwhm_g, fwhm_e, r_ee50, r_ee80 are already in mas
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

    for j in range(strehl.shape[1]):
        print(f"\nWavelength: {wl_nm[j]:.0f} nm")
        print(header)
        print(sep)
        for i in order:
            print(
                f"{x[i, j]:>12.4f}  {y[i, j]:>12.4f}  {r[i]:>12.4f}"
                f"  {strehl[i, j]:>8.4f}"
                f"  {fwhm_g[i, j]:>13.2f}"
                f"  {fwhm_e[i, j]:>13.2f}"
                f"  {r_ee50[i, j]:>11.2f}"
                f"  {r_ee80[i, j]:>11.2f}"
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

    x, y, w, strehl, fwhm_g, fwhm_e, r_ee50, r_ee80 = (
        maos_utils.get_psf_metrics_over_field(
            str(sim_dir) + "/", seed=args.seed, oversamp=args.oversamp
        )
    )
    print_results(x, y, w, strehl, fwhm_g, fwhm_e, r_ee50, r_ee80)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
