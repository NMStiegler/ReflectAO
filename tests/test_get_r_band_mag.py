from pathlib import Path

import pytest

from reflectao.build_observation_table import build_observation_table

from reflectao.telemetry_utils import get_tt_guide_star_r_mag

def test_get_tt_guide_star_r_mag(tmp_path: Path):
    TEST_PATH = Path("/g3/data/kapa/2026feb26/raw/i260226_a010002.fits")
    tbl = build_observation_table(TEST_PATH, instrument="OSIRIS")
    r_mag = get_tt_guide_star_r_mag(TEST_PATH, tbl['t_exposure_start'][0], tbl['t_exposure_duration'][0])
    assert r_mag == 15.61
