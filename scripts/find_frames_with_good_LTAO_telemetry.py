# Import useful packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import reflectao as rao
from reflectao import telemetry_utils as tu
from astropy import units as u
from astropy.time import Time

# Pick the night to analyze
night = tu.good_kapa_night_list[5]
print(night)

# Get telemetry filenames for all images taken on this night
path_to_telemetry_file_dirs = tu.get_path_to_telemetry_dir(night) / "IMAG"
all_telemetry_file_dirs = np.array(list(path_to_telemetry_file_dirs.glob("*.fits")))
print(f"Found {len(all_telemetry_file_dirs)} telemetry files for night {night}.")

# Find the paths of the telemetry .npz files for each image
telemetry_files = np.array([tu.read_image_telemetry(dir_name) for dir_name in all_telemetry_file_dirs], dtype='object')

# Check that there's telemetry data in the image telemetry folders
num_telemetry_files_per_image = np.array([len(telem_files) for telem_files in telemetry_files])
print(f"For {night}, number of telemetry files per image:\n{num_telemetry_files_per_image}")

# Check if there is an ocam2k telemetry file for each image
has_ocam2k_telemetry = np.array([tu.has_ocam2k_data(telem_file) for telem_file in telemetry_files])
print(sum(has_ocam2k_telemetry), "have ocam2k telemetry")

# Load in the WFS telemetry file for each image
ocam2k_telemetry_files = [tu.load_ocam2k_data(telem_files) for telem_files in telemetry_files[has_ocam2k_telemetry]]
print(f"Number of Ocam2k telemetry files loaded: {len(ocam2k_telemetry_files)}")
print(f"Shape of first Ocam2k telemetry file: {ocam2k_telemetry_files[0][tu.get_ocam2k_intensity_key(1)].shape}")

# Look for ocam2k telemetry with data for four WFSs (KAPA)
four_lgs_images = []
four_lgs_image_indices = []
for i, file in enumerate(ocam2k_telemetry_files):
    if tu.has_four_lgs_data(file):
        four_lgs_images.append(all_telemetry_file_dirs[has_ocam2k_telemetry][i])
        four_lgs_image_indices.append(i)

four_lgs_image_mask = np.zeros(len(all_telemetry_file_dirs[has_ocam2k_telemetry]), dtype=bool)
four_lgs_image_mask[four_lgs_image_indices] = True

print(f"Found {len(four_lgs_images)} images with SHWFS4 telemetry (KAPA).")

# Load in all the headers for each of the images associated with this telemetry
image_paths = [tu.get_image_path_from_telemetry_path(telem_file) for telem_file in all_telemetry_file_dirs[has_ocam2k_telemetry]]
hdr_tbl = rao.build_observation_table(image_paths, "OSIRIS", verbose=False)
print("Length of the table is:", len(hdr_tbl), "entries")

# Get the indices for images with WFE < 1000nm (closed-loop ish)
wfe_values = np.array(hdr_tbl["lgs_rms_wfe"].value.tolist()) * u.nm
closed_loop_mask = wfe_values < 1000 * u.nm
print(f"Found {sum(closed_loop_mask)} images with RMS LGS WFE < 1000 nm.")

# Define a function which tells us which sets are not good for a given date
def date_dependent_setup_mask(set_numbers, date):
    if date.to_datetime().date() == Time("2026-02-28").to_datetime().date():
        return np.array(set_numbers) < 6
    elif date.to_datetime().date() == Time("2026-03-04").to_datetime().date():
        return np.array(set_numbers) < 5
    elif date < Time("2026-03-04"):
        return np.array(set_numbers) == 0 | np.array(set_numbers) == 2
    else:
        return np.array(set_numbers) == 0

# Get mask for images with set numbers known to be darks or other non-standard images
set_numbers = [int(num) for num in hdr_tbl['set_number']]
dates = [Time(date) for date in hdr_tbl['t_exposure_start']]
setup_mask = np.array([date_dependent_setup_mask(set_numbers[i], dates[i]) for i in range(len(set_numbers))])
not_setup_mask = ~setup_mask
print(f"Found {sum(not_setup_mask)} images with set number != 0 and != 2")

# Mask everything @ ra & dec of 0 (defaults means not real images of targets)
ra_values = np.array(hdr_tbl['target_ra'].value.tolist())
dec_values = np.array(hdr_tbl['target_dec'].value.tolist())
ra_mask = ra_values == 0
dec_mask = dec_values == 0
target_mask = ra_mask & dec_mask
not_target_mask = ~target_mask
print(f"Found {sum(not_target_mask)} images with target ra and dec != 0.")

# And which is not at a target name of 'stow' or 'wfsblank1'
target_names = [str(name) for name in hdr_tbl['target_name']]
stow_mask = np.array(['stow' in str(name) for name in target_names])
wfsblank1_mask = np.array(['wfsblank1' in str(name) for name in target_names])
not_stow_mask = ~(stow_mask | wfsblank1_mask)
print(f"Found {sum(not_stow_mask)} images with target name not containing 'stow' or 'wfsblank1'.")

# Mask anything which has 'darks' in the set name
set_names = [str(name) for name in hdr_tbl['set_name']]
darks_mask = np.array(['darks' in str(name) for name in set_names])
not_darks_mask = ~darks_mask
print(f"Found {sum(not_darks_mask)} images with set name not containing 'darks'.")

# And which don't have pinhole mask as an object name
object_names = [str(name) for name in hdr_tbl['object_name']]
pinhole_mask = np.array(['pinhole' in str(name) for name in object_names])
not_pinhole_mask = ~pinhole_mask
print(f"Found {sum(not_pinhole_mask)} images with object name not containing 'pinhole'.")

# And which haven't aborted
aborted_mask = np.array(hdr_tbl['aborted'].value.tolist())
not_aborted_mask = ~aborted_mask
print(f"Found {sum(not_aborted_mask)} images with aborted flag not set.")

from collections import defaultdict

def get_set_dictionary(table):
    # 1. Convert columns to int directly on the table object
    # No need for .data.data - Astropy columns handle .astype() directly
    assert(sum(table['set_number'].mask) == 0), "There are masked values in the 'set_number' column. Please handle them before converting to int."
    assert(sum(table['frame_number'].mask) == 0), "There are masked values in the 'frame_number' column. Please handle them before converting to int."

    table['set_number'] = table['set_number'].filled().astype(int)
    table['frame_number'] = table['frame_number'].filled().astype(int)

    # 2. Sort numerically
    table.sort(['set_number', 'frame_number'])
    
    set_dict = defaultdict(list)
    for row in table:
        set_dict[row['set_number']].append(row['frame_number'])
        
    return dict(set_dict)

# Mask everything that doesn't meet the criteria
good_mask = four_lgs_image_mask & closed_loop_mask & not_setup_mask & not_target_mask & not_darks_mask & not_stow_mask & not_aborted_mask & not_pinhole_mask

print(f"There are {sum(good_mask)} images with all required criteria")

# 1. Extract the subset into its own table variable
sub_table = hdr_tbl[good_mask]

# Generate the dictionary
set_frames_map = get_set_dictionary(sub_table)

# Print with a custom format
print("{")
for key, value in set_frames_map.items():
    print(f" {key}: {value},")
print("}")