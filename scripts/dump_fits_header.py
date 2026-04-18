"""
Dump all FITS HDU headers to a text file under ReflectAO/temp.

Default input path matches the README example. You can override it with --fits.
"""

from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_FITS = "/g3/data/kapa/2025dec04/raw/i251204_a001002.fits"


def dump_headers_to_txt(*, fits_path: Path, out_path: Path) -> None:
    try:
        from astropy.io import fits  # type: ignore
    except Exception as e:  # pragma: no cover
        raise SystemExit(
            "Missing dependency: astropy. Install it in your environment (see README)."
        ) from e

    if not fits_path.exists():
        raise SystemExit(f"FITS file not found: {fits_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with fits.open(fits_path, memmap=True) as hdul:
        lines: list[str] = []
        for idx, hdu in enumerate(hdul):
            hdu_name = getattr(hdu, "name", None) or ""
            header = hdu.header

            title = f"HDU {idx}"
            if hdu_name:
                title += f" ({hdu_name})"

            lines.append("=" * 80)
            lines.append(title)
            lines.append("=" * 80)
            lines.append(header.tostring(sep="\n", endcard=True, padding=False))
            lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(
        description="Write full FITS header(s) to ReflectAO/temp/*.txt"
    )
    p.add_argument(
        "--fits",
        default=DEFAULT_FITS,
        help="Path to input .fits (default: README example path)",
    )
    p.add_argument(
        "--out",
        default=None,
        help="Output .txt path (default: ReflectAO/temp/<fits_stem>_header.txt)",
    )
    args = p.parse_args()

    fits_path = Path(args.fits).expanduser()

    if args.out is None:
        out_path = Path(__file__).resolve().parents[1] / "temp" / f"{fits_path.stem}_header.txt"
    else:
        out_path = Path(args.out).expanduser()

    dump_headers_to_txt(fits_path=fits_path, out_path=out_path)
    print(f"Wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
