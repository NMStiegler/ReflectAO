# Note: Fails on /g3/data/kapa/2025dec06/telemetry/IMAG/i251206_a001002.fits

# Import useful packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, LogNorm, Normalize
from pathlib import Path
import reflectao as rao
from reflectao import kapa_utils as ku
from reflectao import telemetry_utils as tu
import astropy.units as u
from astropy.time import Time
import paarti.utils.maos_utils as mu

# Loop through all good images in all sets for all nights to build a big header
# table we can use to collect all information

# Get all the FITS file and telemetry file paths
count_of_all_images = 0 # Count how many good images we have
fits_paths = []
telemetry_paths = []
for night in tu.good_kapa_night_list:
    set_numbers = tu.good_kapa_img_wfs_telemetry[night].keys()
    for set_number in set_numbers:
        frame_numbers = tu.good_kapa_img_wfs_telemetry[night][set_number]
        for frame_number in frame_numbers:
            print(f"Configuring night {night}, set {set_number}, frame {frame_number}...")
            # if count_of_all_images > 10:
            #     print(f"Reached 10 images, stopping early for testing purposes. Remove this condition to process all images.")
            #     break
            count_of_all_images += 1
            fits_path = tu.get_path_to_image(night, set_number, frame_number)
            telemetry_path = tu.get_path_to_image_telemetry(night, set_number, frame_number)
            fits_paths.append(fits_path)
            telemetry_paths.append(telemetry_path)

# Build the table
hdr_tbl = rao.build_observation_table(fits_paths, instrument="OSIRIS", telemetry_paths=telemetry_paths, verbose=True)

print(f"Total number of good images to analyze: {count_of_all_images}")
print("-" * 10, f"Done creating table at {Time.now().to_datetime()}, now analyzing files", "-" * 10)

# For each good image, get the median of the mean and standard deviation of the
# electrons recorded per read for illuminated and unilluminated subapertures for
# all four WFSs. These will be lists of lists of 4 values, one for each WFS, for each image
dim_mean_electrons = [] # Collect the median of the per subaperture mean electrons of the dim population.
dim_std_electrons = [] # Collect the median of the per subaperture standard deviations of the dim population.
bright_mean_electrons = [] # Collect the median of the per subaperture mean electrons of the bright population.
bright_std_electrons = [] # Collect the median of the per subaperture standard deviations of the bright population.

# Helper function to safely extract values whether they are Quantities or standard numbers
def strip_units(data_list):
    return [x.to_value(u.electron) if hasattr(x, 'to_value') else x for x in data_list]

for i, row in enumerate(hdr_tbl):
    # Print info
    print(f"Analyzing {i} / {len(hdr_tbl)} night {str(row['t_exposure_start'].to_datetime().date())}, set {row['set_number']}, frame {row['frame_number']} at {Time.now().to_datetime()}")

    # Load the telemetry and fits header for that image
    path_to_telemetry = Path(row["telemetry_file_path"])
    telem_files = tu.read_image_telemetry(path_to_telemetry, verbose=True)
    ocam2k = tu.load_ocam2k_data(telem_files)

    # Compute the mean and standard deviation of the number of electrons read per frame
    # for each subaperture across all 4 WFSs
    sensor_mean_electrons, sensor_stds_electrons = tu.compute_aperture_wise_electron_stats(ocam2k, row)
    median_dim_mean_e, median_bright_mean_e = tu.compute_median_values_for_all_wfs(sensor_mean_electrons, thresh=20)
    median_dim_std_e, median_bright_std_e = tu.compute_median_values_for_all_wfs(sensor_stds_electrons, thresh=20)

    # Add to lists as unitless numpy arrays by iterating through the returned lists
    dim_mean_electrons.append(np.array(strip_units(median_dim_mean_e)))
    dim_std_electrons.append(np.array(strip_units(median_dim_std_e)))
    bright_mean_electrons.append(np.array(strip_units(median_bright_mean_e)))
    bright_std_electrons.append(np.array(strip_units(median_bright_std_e)))

print("-" * 10, f"Done analyzing files at {Time.now().to_datetime()}, now saving to table", "-" * 10)

# Convert lists of 1D arrays into 2D numpy arrays
# Shape will be (N_images, 4) where 4 is the number of WFSs
dim_mean_arr = np.array(dim_mean_electrons)
dim_std_arr = np.array(dim_std_electrons)
bright_mean_arr = np.array(bright_mean_electrons)
bright_std_arr = np.array(bright_std_electrons)

# Put arrays into hdr_tbl as separate columns for each WFS (0 through 3)
num_wfs = 4
for i in range(num_wfs):
    # Assign the columns
    hdr_tbl[f'wfs{i}_dim_mean_electrons'] = dim_mean_arr[:, i]
    hdr_tbl[f'wfs{i}_dim_std_electrons'] = dim_std_arr[:, i]
    hdr_tbl[f'wfs{i}_bright_mean_electrons'] = bright_mean_arr[:, i]
    hdr_tbl[f'wfs{i}_bright_std_electrons'] = bright_std_arr[:, i]
    
    # Reapply the Astropy units to the table columns directly so metadata is preserved
    hdr_tbl[f'wfs{i}_dim_mean_electrons'].unit = u.electron
    hdr_tbl[f'wfs{i}_dim_std_electrons'].unit = u.electron
    hdr_tbl[f'wfs{i}_bright_mean_electrons'].unit = u.electron
    hdr_tbl[f'wfs{i}_bright_std_electrons'].unit = u.electron

# Save hdr_tbl
place_to_save = "/u/nstieg/work/ao/keck/kapa/telemetry/kapa_ocam2k_read_stats.ecsv"
hdr_tbl.write(place_to_save, format="ascii.ecsv", overwrite=True)