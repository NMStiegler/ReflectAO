"""
Logical schema for ReflectAO's per-frame observation table.

This module defines the column set ("schema") for the canonical internal table:
one exposure corresponds to one row, and every schema field always has a column.

Missing values are represented using **masked** table entries (not Python
``None``). This lets us keep columns typed consistently (floats, strings,
quantities, etc.) while still marking missing data in a standard way.
"""

from astropy.table import MaskedColumn, Table


class ColumnDef(object):
    """
    Definition of one logical schema column.

    :param name: Column name in the canonical internal table.
    :type name: str
    :param unit: Canonical unit for the column, stored as metadata.
                 This is a string for now because many entries are not strict
                 `astropy.units.Unit` expressions.
    :type unit: str or None
    :param meaning: Human-readable description of the field.
    :type meaning: str
    """

    def __init__(self, name, unit, meaning):
        self.name = name
        self.unit = unit
        self.meaning = meaning


# NOTE: This list is copied from `docs/relevant_fields.md` (Field, Units, Meaning).
# It is intentionally verbose and explicit: schema changes should be obvious diffs.
SCHEMA = (
    # Information from the system where the script was run
    ColumnDef(
        name="image_path",
        unit="string",
        meaning="Filesystem path or URI of the data product this row was derived from, for provenance.",
    ),
    
    # Telescope / instrument / set information
    ColumnDef(
        name="telescope_name",
        unit="string",
        meaning="Identifier for the telescope (e.g. facility name).",
    ),
    ColumnDef(
        name="instrument_name",
        unit="string",
        meaning="Name of the instrument (distinct from the telescope or facility identifier).",
    ),
    ColumnDef(
        name="telescope_elevation",
        unit="deg",
        meaning="Pointing elevation angle of the telescope when the observation was taken",
    ),
    ColumnDef(
        name="telescope_azimuth",
        unit="deg",
        meaning="Pointing azimuth angle of the telescope when the observation was taken",
    ),
    ColumnDef(
        name="plate_scale",
        unit="arcsec / pixel",
        meaning="Plate scale as an angular size per pixel",
    ),
    ColumnDef(
        name="frame_number",
        unit="int",
        meaning="Running frame or sequence number within a dataset",
    ),
    ColumnDef(
        name="set_number",
        unit="int",
        meaning="Running set or sequence number within a dataset",
    ),
    ColumnDef(
        name="set_name",
        unit="string",
        meaning="Name of the set or sequence",
    ),
    
    # Target information
    ColumnDef(
        name="target_name",
        unit="string",
        meaning="Target name (from tracking system)",
    ),
    ColumnDef(
        name="object_name",
        unit="string",
        meaning="Object name (from observer's notes)",
    ),
    ColumnDef(
        name="target_ra",
        unit="deg",
        meaning="Right ascension of the target",
    ),
    ColumnDef(
        name="target_dec",
        unit="deg",
        meaning="Declination of the target",
    ),
    ColumnDef(
        name="epoch",
        unit="string",
        meaning="Epoch of the recorded RA/Dec",
    ),
    
    
    # Exposure information
    ColumnDef(
        name="t_exposure_start",
        unit="UTC timestamp",
        meaning="Best-effort start of the on-sky integration from the primary date and time keywords (calendar date and time-of-day combined).",
    ),
    ColumnDef(
        name="t_exposure_duration",
        unit="s",
        meaning="True, total on-sky integration time for the frame (should be true integration time * num_coadds)",
    ),
    ColumnDef(
        name="t_int",
        unit="s",
        meaning="(true) integration time per coadd",
    ),
    ColumnDef(
        name="num_coadds",
        unit="int",
        meaning="Number of coadds or equivalent stacking count",
    ),
    ColumnDef(
        name="wavelength",
        unit="nm",
        meaning="Effective or band-center wavelength of the observation",
    ),
    ColumnDef(
        name="filter_name",
        unit="string",
        meaning="Name of the filter",
    ),
    ColumnDef(
        name="airmass",
        unit="dimensionless",
        meaning="Atmospheric airmass along the line of sight.",
    ),
    ColumnDef(
        name="aborted",
        unit="boolean",
        meaning="Whether the exposure was stopped before completion",
    ),
    ColumnDef(
        name="dither_name",
        unit="string",
        meaning="Name of the dither sequence, if used",
    ),
    
    # AO system information
    
    ColumnDef(
        name="lgs_wfs_rate",
        unit="Hz",
        meaning="Sampling/frame rate of the high order laser guide star wavefront sensor loop",
    ),
    # ColumnDef(
    #     name="tt_wfs_rate",
    #     unit="Hz",
    #     meaning="Sampling/frame rate of the low order tip tilt wavefront sensor loop",
    # ),
    ColumnDef(
        name="lgs_rms_wfe",
        unit="nm",
        meaning="Average RMS residual wavefront error of the LGS WFS for the observation",
    ),
    ColumnDef(
        name="lgs_layer_alt",
        unit="m",
        meaning="Altitude of the sodium layer",
    ),
    ColumnDef(
        name="ngs_fwhm",
        unit="arcsec",
        meaning="FWHM size of the TT-NGS on the sky for the observation", # Is this true/logged?
    ),
    ColumnDef(
        name="ngs_wavelength",
        unit="nm",
        meaning="Effective wavelength of the TT-NGS star", # Is this in the headers/logged? Do we need it?
    ),
    # ColumnDef(
    #     name="t_int_ngs",
    #     unit="s",
    #     meaning="Integration time of the sensor looking at the TT-NGS",
    # ),
    ColumnDef(
        name="reconstructor_name",
        unit="string",
        meaning="Which reconstructor was used for the observation (full LTAO, semi-GLAO, etc)",
    ),
    ColumnDef(
        name="dm_gain",
        unit="dimensionless",
        meaning="Gain of the DM loop",
    ),
    ColumnDef(
        name="lgs_wfs_gain",
        unit="dimensionless",
        meaning="Gain of the LGS WFS loop",
    ),
    ColumnDef(
        name="system_gain",
        unit="dimensionless",
        meaning="System gain factor",
    ),
    ColumnDef(
        name="utt_gain",
        unit="dimensionless",
        meaning="Gain of the up tip tilt loop",
    ),
    ColumnDef(
        name="dtt_gain",
        unit="dimensionless",
        meaning="Gain of the down tip tilt loop",
    ),
    ColumnDef(
        name="ao_mode",
        unit="string",
        meaning="AO operations or control mode code/name.", # What is this / where did it come from? Same as reconstructor_name?
    ),
    
    # Weather information
    ColumnDef(
        name="humidity_dome",
        unit="percent relative humidity (float)",
        meaning="Relative humidity inside the dome",
    ),
    ColumnDef(
        name="T_dome_air",
        unit="°C",
        meaning="Air temperature inside the dome",
    ),
    ColumnDef(
        name="humidity_outside",
        unit="percent relative humidity (float)",
        meaning="Relative humidity outside the dome",
    ),
    ColumnDef(
        name="T_outside_air",
        unit="°C",
        meaning="Air temperature outside the dome",
    ),
    ColumnDef(
        name="P_barometric",
        unit="Pa",
        meaning="Barometric pressure at the site",
    ),
    ColumnDef(
        name="wind_direction",
        unit="deg",
        meaning="Wind direction",
    ),
    ColumnDef(
        name="wind_speed",
        unit="m/s",
        meaning="Wind speed",
    ),
    ColumnDef(
        name="t_weather_sample",
        unit="string",
        meaning="Time stamp associated with the weather block",
    ),
    ColumnDef(
        name="T_tube",
        unit="°C",
        meaning="Temperature inside the telescope tube",
    )
)


def schema_column_names(schema=SCHEMA):
    """
    Return the required schema column names in schema order.

    :param schema: Iterable of `ColumnDef` objects.
    :type schema: iterable
    :return: List of column names in order.
    :rtype: list[str]
    """
    return [c.name for c in list(schema)]


def new_empty_observation_table(n_rows=0, schema=SCHEMA):
    """
    Create a new observation table with all schema columns present.

    The table is created with `MaskedColumn`s so that missing values are stored
    as **masked entries**, not as ``None``.

    :param n_rows: Number of rows to allocate. All values start masked.
    :type n_rows: int
    :param schema: Column definitions to use.
    :type schema: iterable
    :return: Table with all schema columns present, in schema order.
    :rtype: astropy.table.Table
    """
    tbl = Table()

    # Table-level metadata is useful for provenance and future schema migration.
    tbl.meta["reflectao_schema_version"] = 0
    tbl.meta["reflectao_schema_source"] = "docs/relevant_fields.md"

    for col in schema:
        # For now, we use dtype=object because we are not yet enforcing per-field
        # numeric/string/quantity typing. The important behavior for this slice
        # is: (1) missing values are masked, and (2) columns always exist.
        #
        # Later, we can tighten dtypes field-by-field once we are confident in
        # the instrument keyword conventions and unit conversions.
        tbl[col.name] = MaskedColumn([None] * n_rows, mask=[True] * n_rows, dtype=object)
        tbl[col.name].meta["meaning"] = col.meaning
        if col.unit is not None:
            tbl[col.name].meta["unit"] = col.unit

    return tbl


def validate_table_has_schema(table, schema=SCHEMA, allow_extra_columns=True):
    """
    Validate that a table contains (at least) the schema-defined columns.

    This is intentionally a strict, readable check: it does not attempt to infer
    or rename columns. If a required column is missing, we raise a ValueError.

    :param table: Table to validate.
    :type table: astropy.table.Table
    :param schema: Schema definitions.
    :type schema: iterable
    :param allow_extra_columns: If True, allow columns beyond the schema; only
        require that all schema columns exist. If False, require an exact match
        of column sets.
    :type allow_extra_columns: bool, optional
    :raises ValueError: If required schema columns are missing, or if extra
        columns are present when `allow_extra_columns` is False.
    :return: None
    :rtype: None
    """
    required = schema_column_names(schema)
    missing = [c for c in required if c not in table.colnames]
    if missing:
        raise ValueError(
            "Provided table does not match ReflectAO schema. Missing required columns: "
            + ", ".join(missing)
        )

    if not allow_extra_columns:
        extra = [c for c in table.colnames if c not in required]
        if extra:
            raise ValueError(
                "Provided table has columns not in the ReflectAO schema: " + ", ".join(extra)
            )
