"""
Read a finished MAOS simulation directory and print the Strehl ratio at each
field position and wavelength.

Usage:
    python print_maos_strehl.py <sim_directory> [--seed SEED]

MAOS normalises closed-loop PSF cubes so that the peak pixel value equals the
Strehl ratio, so no further calculation is needed beyond reading the FITS data.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


def _position_from_header(hdr) -> tuple[float, float]:
    """Return (x_arcsec, y_arcsec) from a MAOS PSF HDU header.

    MAOS stores the evaluation direction as a FITS complex keyword 'THETA'
    where the real part is the y offset and the imaginary part is the x offset,
    both in arcseconds.
    """
    theta = hdr.get("THETA")
    if theta is None:
        raise KeyError("THETA")
    if isinstance(theta, complex):
        return float(theta.imag), float(theta.real)
    # Astropy may return a string '(re,im)' for FITS complex keywords.
    m = re.match(r"\(\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*\)", str(theta))
    if m:
        return float(m.group(2)), float(m.group(1))  # imag=x, real=y
    raise ValueError(f"Cannot parse THETA keyword: {theta!r}")


def _position_from_filename(path: Path) -> tuple[float, float]:
    """Parse (x, y) field position in arcsec from filename 'evlpsfcl_N_xX_yY.fits'."""
    m = re.search(r"_x(-?[\d]+)_y(-?[\d]+)\.fits$", path.name)
    if not m:
        return 0.0, 0.0
    return float(m.group(1)), float(m.group(2))


def read_maos_results(
    directory: Path, seed: int
) -> tuple[list[tuple[float, float]], np.ndarray, np.ndarray]:
    """Return positions, wavelengths (micron), and Strehl arrays.

    Returns
    -------
    positions : list of (x_arcsec, y_arcsec) tuples, length N_pos
    wavelengths : ndarray shape (N_wvl,) in microns
    strehls : ndarray shape (N_pos, N_wvl)
    """
    fits_files = sorted(directory.glob(f"evlpsfcl_{seed}_x*_y*.fits"))
    if not fits_files:
        print(
            f"Error: no files matching 'evlpsfcl_{seed}_x*_y*.fits' in {directory}",
            file=sys.stderr,
        )
        sys.exit(1)

    positions: list[tuple[float, float]] = []
    wavelengths: np.ndarray | None = None
    strehls: list[np.ndarray] = []

    for fpath in fits_files:
        with fits.open(fpath) as hdul:
            n_wvl = len(hdul)
            wvls = np.zeros(n_wvl)
            strehl_row = np.zeros(n_wvl)

            try:
                pos = _position_from_header(hdul[0].header)
            except (KeyError, ValueError):
                pos = _position_from_filename(fpath)
            positions.append(pos)

            for j, hdu in enumerate(hdul):
                hdr = hdu.header
                # WVL is stored in metres; convert to microns for display
                wvls[j] = hdr.get("WVL", 0.0) * 1e6
                strehl_row[j] = float(hdu.data.max())

        if wavelengths is None:
            wavelengths = wvls
        strehls.append(strehl_row)

    if wavelengths is None:
        wavelengths = np.array([])

    return positions, wavelengths, np.array(strehls)


def print_results(
    positions: list[tuple[float, float]],
    wavelengths: np.ndarray,
    strehls: np.ndarray,
) -> None:
    n_wvl = len(wavelengths)
    wl_nm = wavelengths * 1e3  # micron -> nm

    col_w = 10
    wl_col_w = 14

    header = (
        f"{'X (arcsec)':>{col_w}}  {'Y (arcsec)':>{col_w}}  {'R (arcsec)':>{col_w}}"
    )
    for wl in wl_nm:
        header += f"  {f'Strehl@{wl:.0f}nm':>{wl_col_w}}"
    print(header)
    print("-" * len(header))

    radii = [np.hypot(x, y) for x, y in positions]
    order = np.argsort(radii)

    for i in order:
        x, y = positions[i]
        r = radii[i]
        line = f"{x:>{col_w}.4f}  {y:>{col_w}.4f}  {r:>{col_w}.4f}"
        for j in range(n_wvl):
            line += f"  {strehls[i, j]:>{wl_col_w}.4f}"
        print(line)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Print MAOS simulation Strehl at each field position and wavelength. "
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
    args = parser.parse_args()

    sim_dir = Path(args.directory).expanduser().resolve()
    if not sim_dir.is_dir():
        print(f"Error: directory not found: {sim_dir}", file=sys.stderr)
        return 1

    positions, wavelengths, strehls = read_maos_results(sim_dir, seed=args.seed)
    print_results(positions, wavelengths, strehls)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
