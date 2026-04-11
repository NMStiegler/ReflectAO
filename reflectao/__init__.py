"""
ReflectAO

This package contains the core utilities for building a universal, per-exposure
metadata representation (one row per FITS file) to support adaptive optics
digital-twin simulations and comparisons.

The public surface area is intentionally small at this stage.
"""

from .build_observation_table import build_observation_table
from .schema import SCHEMA, ColumnDef, new_empty_observation_table, validate_table_has_schema

__all__ = [
    "SCHEMA",
    "ColumnDef",
    "build_observation_table",
    "new_empty_observation_table",
    "validate_table_has_schema",
]

