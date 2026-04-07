# Acronyms: AO — adaptive optics; FITS — Flexible Image Transport System; ReflectAO — on-sky vs simulation comparison tooling (name plays on “reflect” + AO).

"""ReflectAO: observation metadata, atmosphere inputs, and simulation glue.

The public surface for FITS-based observation metadata is
:class:`reflectao.observation.FrameMetadata` and :func:`reflectao.observation.read_fits_header`.
Instrument-specific adapters live under :mod:`reflectao.observation.instruments`.
"""

from __future__ import annotations

from reflectao.observation import (
    FrameMetadata,
    HeaderKeywordMap,
    KeywordResolutionPolicy,
    Nirc2HeaderAdapter,
    OsirisHeaderAdapter,
    read_fits_header,
)

__all__ = [
    "FrameMetadata",
    "HeaderKeywordMap",
    "KeywordResolutionPolicy",
    "Nirc2HeaderAdapter",
    "OsirisHeaderAdapter",
    "read_fits_header",
]
