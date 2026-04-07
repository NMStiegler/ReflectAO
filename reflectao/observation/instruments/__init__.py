# Acronyms: NIRC2 — Near-Infrared Camera 2 (Keck II); OSIRIS — OH-Suppressing Infra-Red Imaging Spectrograph.

"""Instrument-specific :class:`~reflectao.observation.fits_header.HeaderAdapter` implementations."""

from __future__ import annotations

from reflectao.observation.instruments.nirc2 import (
    NIRC2_HEADER_KEYWORDS,
    Nirc2HeaderAdapter,
)
from reflectao.observation.instruments.osiris import (
    OSIRIS_HEADER_KEYWORDS,
    OsirisHeaderAdapter,
)

__all__ = [
    "NIRC2_HEADER_KEYWORDS",
    "Nirc2HeaderAdapter",
    "OSIRIS_HEADER_KEYWORDS",
    "OsirisHeaderAdapter",
]
