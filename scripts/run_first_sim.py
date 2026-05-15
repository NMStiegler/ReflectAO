from pathlib import Path
import reflectao as rao
from reflectao import telemetry_utils as tu
from reflectao import run_sim as rs

# Pick the night to analyze
night = "2026mar04"
set = 21
image = 2

# Get the path to the image file
path_to_image_file = Path(tu.get_path_to_image(night, set, image))

# Make an observation table collecting all the data about this image
hdr_tbl = rao.build_observation_table(path_to_image_file, "OSIRIS", verbose=True)

# Use maos to simulate the AO performance
rs.run_maos_sim(hdr_tbl[0])