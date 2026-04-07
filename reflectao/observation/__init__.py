# Acronyms: FITS — Flexible Image Transport System; HDU — Header/Data Unit (FITS extension).

"""On-sky observation metadata for simulation setup (no pixel science products).

Public objects:

- :class:`FrameMetadata` — logical per-frame metadata.
- :func:`read_fits_header` — read a FITS header from disk for a supported instrument.
- :class:`OsirisHeaderAdapter`, :class:`Nirc2HeaderAdapter` — parse in-memory
  :class:`astropy.io.fits.Header` instances.
- :class:`HeaderKeywordMap`, :class:`KeywordResolutionPolicy` — configuration
  for custom adapters (see :mod:`reflectao.observation.fits_header` for
  :class:`DefaultHeaderAdapter`, which is not re-exported here).

"""

from __future__ import annotations

from reflectao.observation.fits_header import (
    HeaderAdapter,
    HeaderKeywordMap,
    KeywordResolutionPolicy,
    read_fits_header,
)
from reflectao.observation.frame_metadata import FrameMetadata
from reflectao.observation.instruments import Nirc2HeaderAdapter, OsirisHeaderAdapter

__all__ = [
    "FrameMetadata",
    "HeaderAdapter",
    "HeaderKeywordMap",
    "KeywordResolutionPolicy",
    "Nirc2HeaderAdapter",
    "OsirisHeaderAdapter",
    "read_fits_header",
]
