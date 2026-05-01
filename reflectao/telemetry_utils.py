"""
Utilities for working with data from the Keck I telemetry system, specifically
focused on the post-KAPA upgrade

Noah Stiegler
4/15/26

"""

# Import useful packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, LogNorm, Normalize
from pathlib import Path
from reflectao import kapa_utils as ku
import reflectao as rao
import astropy.units as u
from astropy.time import Time
from astropy.table import Table, Row
import paarti.utils.maos_utils as mu

### Useful References ###

# List of nights with telemetry available
good_kapa_night_list = ["2025nov06",
                 #   "2025nov08", <-- Telemetry doesn't match data (which is is in a weird place too)
                   "2025dec04",
                   "2025dec06",
                   "2026feb26",
                   "2026feb28",
                 #   "2026jan12", <-- No telemetry (TRS server down)
                 #   "2026jan31", <-- Telemetry files are empty (lingering TRS issue)
                   "2026mar04"]

# Some nights didn't record t_exposure_duration right. Keeping track
nights_with_exposure_duration_issues = [
    "2025nov06", # All images have t_exposure_duration that is much longer than t_int * num_coadds <(5000.0 s) does not match t_int * num_coadds (4.425 s)>
    "2025dec04", # Likewise here, although I think it was a unit conversion error because everything is 5x longer than it should be
    "2025dec06", # again
]

# Some nights have data stored in a slightly different format than usual
nights_with_weird_image_paths = [
    "2025nov08", # Images are in raw/ instead of raw_images/
]

# Darks
# 2/28/26 sets 2 & 36
# 3/4/26 set 3

# For each night, we want to know which sets and which images have good ocam2k files to
# analyze. Here we'll organize that information by dictionaries where keys are
# night dates (as YYYYmonDD strings) and values are dictionaries with set_num
# keys and lists of their images with good telemetry files as values. Each image
# listed here is one that has telemetry from 4 laser guide stars and wfe < 1000nm
# (ensuring closed loop). I tried to be generous on the closed loop criterion
# so this might not be the *best* telemetry
# Data generated from ReflectAO/scripts/find_frames_with_good_LTAO_telemetry.py
# This is data and should not be autocorrected by GitHub Copilot or other AI tools
good_kapa_img_wfs_telemetry = {
    "2025nov06": {
        1: [117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159],
    },
    "2025nov08": {},
    "2025dec04": {
        3: [2, 3, 4, 5, 6, 7, 8, 9],
        6: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        7: [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13],
        8: [2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123],
    },
    "2025dec06": {
        4: [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
        5: [2, 3, 5, 6, 7, 8, 9, 10, 11],
        6: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        7: [2, 3, 4, 5],
        8: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        9: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],
    },
    "2026feb26": {
        9: [2, 3, 4, 5, 6, 7, 8, 9],
        10: [2, 3, 4, 5, 6, 7, 8, 9],
        11: [2, 3, 4, 5, 6, 7, 8, 9],
        12: [2, 3, 4, 5, 6, 7, 8, 9],
        13: [2, 3, 4, 5, 6, 7, 8, 9],
        14: [2, 3, 4, 5, 6, 7, 8, 9],
        15: [2, 3, 4, 5, 6, 7, 8, 9],
        16: [2, 3, 4, 5, 6, 7, 8, 9],
        17: [2, 3, 4, 5, 6, 7, 8, 9],
        18: [2, 3, 4, 5, 6, 7, 8, 9],
        19: [2, 3, 4, 5, 6, 7, 8, 9],
        20: [2, 3, 4, 5, 6, 7, 8, 9],
        21: [2, 3, 4, 5, 6, 7, 8, 9],
        22: [2, 3, 4, 5, 6, 7, 8, 9],
        23: [2, 3, 4, 5, 6, 7, 8, 9],
        24: [2],
        25: [2],
        26: [2],
        27: [2],
        28: [2, 3, 4, 6],
        29: [2],
        30: [2],
        31: [2],
        32: [2, 3, 4, 5, 6, 7, 8, 9, 10],
        33: [2],
        34: [2, 3, 4, 5],
        35: [2],
        36: [2],
        37: [2],
        38: [2],
        39: [2, 3, 4],
        40: [2],
        41: [3, 4, 5, 6, 7, 8, 9, 10],
        42: [2, 3, 4, 5, 6, 7, 8, 9, 10],
        43: [2, 3, 4, 5, 6, 7],
    },
    "2026feb28": {
        13: [2],
        14: [2],
        15: [2, 3, 4, 5, 6, 7, 8, 9, 10],
        16: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        17: [2],
        18: [2, 3, 4, 5, 6],
        20: [2],
        21: [2, 3, 4, 5, 6, 7, 8, 9, 10],
        23: [2, 3, 4, 5, 6, 7, 8, 9, 10],
        24: [2, 3, 4, 5, 6, 7, 9, 10],
    },
    "2026mar04": {
        5: [2, 3, 4, 5, 6, 7, 8, 9],
        6: [2, 6, 7, 8, 9],
        7: [2, 3, 4, 5, 6, 7, 8, 9],
        8: [2, 3, 4, 5, 6, 7],
        9: [2, 3, 4, 5],
        10: [2, 3, 4, 5, 6, 7, 8, 9],
        11: [2, 3],
        12: [2, 3, 4, 5, 6, 7, 8, 9],
        13: [2, 3, 4, 5, 6, 7, 8, 9],
        14: [2, 3, 4, 5, 6, 7, 8, 9],
        15: [2, 3, 4, 5, 6, 7, 8, 9],
        16: [2, 3, 4, 5, 6, 7, 8, 9],
        17: [2, 3, 4, 5, 6, 7, 8, 9],
        18: [2, 3, 4, 5],
        19: [2, 3, 4, 5],
        20: [2, 3, 4, 5],
        21: [2, 3, 4, 5],
        22: [2, 3, 4, 5],
        23: [2, 3, 4, 5],
        24: [2, 3, 4, 5],
        25: [2, 3, 4, 5],
        26: [2, 4, 5],
        27: [2, 3, 4, 5],
        28: [2, 3, 4, 5],
        29: [2, 3, 4, 5],
        30: [2, 3, 4, 5],
        31: [2, 3, 4, 5],
        32: [2],
        33: [2, 4, 5],
        34: [2, 3, 4, 5],
        35: [2, 3, 4, 5],
        36: [2, 3, 4, 5],
        37: [2],
        38: [2, 3, 4, 5],
        39: [2, 3, 4, 5],
        40: [2, 3, 4, 5],
        41: [2, 3, 4, 5],
        42: [2],
        43: [2],
        44: [2],
        45: [2, 3, 4, 5, 6, 7, 8],
        46: [2, 3, 4, 5],
    }
}

### Internal Helper Functions ###

def _check_night_param(night):
    # Check night is in YYYYmonDD format
    assert(isinstance(night, str)), "Night must be a string in the format 'YYYYmonDD', e.g. '2026mar04'"
    assert(len(night) == 9), "Night must be in the format 'YYYYmonDD', e.g. '2026mar04'"
    assert(night[4:7].isalpha()), "Night must be in the format 'YYYYmonDD', e.g. '2026mar04'"
    assert(night[0:4].isdigit()), "Night must be in the format 'YYYYmonDD', e.g. '2026mar04'"
    assert(night[7:9].isdigit()), "Night must be in the format 'YYYYmonDD', e.g. '2026mar04'"
    assert(night[4:7].lower() in month_dict.keys()), f"Month must be one of {list(month_dict.keys())}, got {night[4:7]}"
    assert(int(night[7:9]) >= 1 and int(night[7:9]) <= 31), "Day must be between 01 and 31"
    assert(int(night[0:4]) >= 2025), "Year must be after 2025" # First commissioning data in 2025

def _check_set_param(set_num):
    assert(isinstance(set_num, int)), "Set must be an integer"

def _check_image_param(image_num):
    assert(isinstance(image_num, int)), "Image must be an integer"

from datetime import datetime

def convert_to_ut_seconds(iso_ts=None):
    """
    Converts OSIRIS style timestamps to UT seconds after midnight.
    """
    # if raw_ts:
    #     # Raw timestamp is in nanoseconds (19 digits)
    #     # Convert to seconds, then get the remainder of a day (86400 seconds)
    #     seconds_since_epoch = raw_ts / 1e9
    #     ut_seconds = seconds_since_epoch % 86400
    #     return round(ut_seconds, 6)

    if iso_ts:
        # Parses 2025-11-08T08:50:27.000425
        dt = datetime.fromisoformat(iso_ts)
        # Calculate seconds since the start of that specific day
        ut_seconds = (dt.hour * 3600) + (dt.minute * 60) + dt.second + (dt.microsecond / 1e6)
        return ut_seconds

# Conversions from 3-letter month to 2-digit month for parsing telemetry file names
month_dict = {
    "jan": "01",
    "feb": "02",
    "mar": "03",
    "apr": "04",
    "may": "05",
    "jun": "06",
    "jul": "07",
    "aug": "08",
    "sep": "09",
    "oct": "10",
    "nov": "11",
    "dec": "12"
}

def _convert_month(month):
    """
    Convert a 3-letter month (ie mar) to a 2-digit month (ie 03)
    
    :param month: The 3-letter month
    :type month: str
    
    :return: The 2-digit month
    :rtype: str
    """
    return month_dict.get(month.lower(), f"Invalid month: {month}")

### Path and file handling functions ###

def get_subap_map_path():
    return Path("/u/bdigia/code/python/paarti/paarti/utils/sub_ap_map.txt")

def load_subap_map():
    subap_map = mu.load_sub_ap_map(subap_map=get_subap_map_path())
    subap_map = np.where(subap_map == 0, np.nan, subap_map) # Set non-subaperture pixels to NaN so they don't mess with the colorbar
    return subap_map

def get_actuator_map_path():
    return Path("/u/bdigia/code/python/paarti/paarti/utils/actuator_map.txt")

def load_actuator_map():
    actuator_map = mu.load_act_map(actuator_map=get_actuator_map_path())
    actuator_map = np.where(actuator_map == 0, np.nan, actuator_map) # Set non-actuator pixels to NaN so they don't mess with the colorbar
    return actuator_map

def get_data_path():
    """
    Get the path to the data directory on the MULab filesystem
    
    :return: The path to the data directory
    :rtype: pathlib.Path
    """
    return Path("/g3/data/kapa/")


def get_path_to_night(night):
    """
    Get the path to the night directory on the MULab filesystem
    
    :return: The path to the night directory
    :rtype: pathlib.Path
    """

    _check_night_param(night)

    return get_data_path() / night

def get_fits_filename(night, set_num, image_num, iso_ts=None):
    """
    Get the fits filename for a given night, set, and image number
    
    :param night: The night of the observation, in the format 'YYYYmonDD', e.g. '2026mar04
    :type night: str'
    :param set: The set number of the observation
    :type set: int
    :param image: The image number of the observation
    :type image: int
    :param iso_ts: The ISO format timestamp from the telemetry file. Only needed for nights with weird image paths (currently just 2025nov08) where the filename includes the timestamp instead of the set and image number.
    :type iso_ts: str
    :return: The fits filename for the given night, set, and image number
    :rtype: str
    """

    _check_night_param(night)
    _check_set_param(set_num)
    _check_image_param(image_num)

    two_digit_year = night[2:4]
    two_digit_month = _convert_month(night[4:7])
    two_digit_day = night[7:9]

    if night not in nights_with_weird_image_paths:
        filename = f"i{two_digit_year}{two_digit_month}{two_digit_day}_a{set_num:03d}{image_num:03d}.fits"
    elif night == "2025nov08":
        four_digit_year = night[0:4]
        UT_seconds = convert_to_ut_seconds(iso_ts=iso_ts)
        filename = f"OI.{four_digit_year}{two_digit_month}{two_digit_day}.{round(UT_seconds, 2)}.fits"

    return filename

def get_path_to_telemetry_dir(night):
    """
    Get the path to the telemetry directory for a given night
    
    :param night: The night of the observation, in the format 'YYYYmonDD', e.g. '2026mar04
    :type night: str'
    
    :return: The path to the telemetry directory for the given night
    :rtype: pathlib.Path
    """

    _check_night_param(night)

    path_to_data = get_data_path()
    path_to_night = path_to_data / night
    return path_to_night / "telemetry"

def get_path_to_image_telemetry(night, set_num, image_num):
    """
    Get the path to the telemetry folder for a given night, set, and image number
    
    :param night: The night of the observation, in the format 'YYYYmonDD', e.g. '2026mar04
    :type night: str'
    :param set: The set number of the observation
    :type set: int
    :param image: The image number of the observation
    :type image: int
    
    :return: The path to the telemetry folder for the given night, set, and image number
    :rtype: pathlib.Path
    """

    _check_night_param(night)
    _check_set_param(set_num)
    _check_image_param(image_num)

    path_to_telemetry = get_path_to_telemetry_dir(night)
    path_to_imaging_telemetry = path_to_telemetry / "IMAG"
    fits_filename = get_fits_filename(night, set_num, image_num)
    return path_to_imaging_telemetry / fits_filename

def get_path_to_image(night, set_num, image_num):
    """
    Get the path to the observed image fits file for a given night, set, and image number

    :param night: The night of the observation, in the format 'YYYYmonDD', e.g. '2026mar04
    :type night: str'
    :param set: The set number of the observation
    :type set: int
    :param image_num: The image number of the observation
    :type image_num: int
    """

    # Check params
    _check_night_param(night)
    _check_set_param(set_num)
    _check_image_param(image_num)

    path_to_night = get_path_to_night(night)
    fits_filename = get_fits_filename(night, set_num, image_num)

    if night not in nights_with_weird_image_paths:
        fits_path = path_to_night / "raw" / fits_filename
    elif night == "2025nov08":
        fits_path = path_to_night / "raw" / "cal" / fits_filename
    return fits_path


def get_image_path_from_telemetry_path(telemetry_path):
    """
    Get the path to the observed image fits file from the path to the telemetry folder for a given night, set, and image number

    :param telemetry_path: The path to the telemetry folder for a given night, set, and image number
    :type telemetry_path: pathlib.Path
    :return: The path to the observed image fits file for the given night, set, and image number
    :rtype: pathlib.Path
    """

    # Check params
    assert(isinstance(telemetry_path, Path)), "telemetry_path must be a pathlib.Path"
    assert(telemetry_path.parent.name == "IMAG"), "telemetry_path must be in the IMAG folder of the telemetry directory"
    assert(telemetry_path.name.startswith("i") and telemetry_path.name.endswith(".fits")), "telemetry_path must be in the format 'iYYMMDD_aSSSIII.fits', e.g. 'i260304_a001078.fits'"

    # Parse night, set, and image number from telemetry path
    fits_filename = telemetry_path.name
    night = telemetry_path.parent.parent.parent.name
    set_num = int(fits_filename[9:12])
    image_num = int(fits_filename[12:15])

    return get_path_to_image(night, set_num, image_num)

def read_image_telemetry(path_to_telemetry_dir, verbose=False):
    """
    Read in the telemetry associated with a given image from the place
    it's stored on the MULab filesystem
    
    :param path_to_telemetry_dir: The path to the telemetry directory
    :type path_to_telemetry_dir: pathlib.Path
    :param verbose: Whether to print verbose output
    :type verbose: bool
    :return: A list of the telemetry file paths associated with the given image
    :rtype: list

    """
    assert(isinstance(path_to_telemetry_dir, Path)), "path_to_telemetry_dir must be a pathlib.Path"

    if verbose: print("Looking at telemetry from", path_to_telemetry_dir)

    # Get the files in the telemetry folder for this image
    telem_files = [file for file in path_to_telemetry_dir.glob("*")]
    if verbose: 
        for file in telem_files:
            print(file)
    
    return telem_files

def get_header(night, set_num, image_num, verbose=False):
    """
    Get the useful parts of the fits file header from observed image

    :param night: The night of the observation, in the format 'YYYYmonDD', e.g. '2026mar04
    :type night: str
    :param set_num: The set number of the observation
    :type set_num: int
    :param image_num: The image number of the observation
    :type image_num: int
    :param verbose: Whether to print verbose output
    :type verbose: bool
    """

    _check_night_param(night)
    _check_set_param(set_num)
    _check_image_param(image_num)

    hdr_tbl = rao.build_observation_table(get_path_to_image(night, set_num, image_num), "OSIRIS", verbose=True)
    return hdr_tbl

def get_telemetry(night, set_num, image_num, verbose=False):
    """
    Get the telemetry file paths for a given night, set, and image number

    :param night: The night of the observation, in the format 'YYYYmonDD', e.g. '2026mar04
    :type night: str
    :param set_num: The set number of the observation
    :type set_num: int
    :param image_num: The image number of the observation
    :type image_num: int
    :param verbose: Whether to print verbose output
    :type verbose: bool
    :return: A list of the telemetry file paths associated with the given night, set, and image number
    :rtype: list
    """

    _check_night_param(night)
    _check_set_param(set_num)
    _check_image_param(image_num)

    path_to_telemetry = get_path_to_image_telemetry(night, set_num, image_num)
    return read_image_telemetry(path_to_telemetry, verbose=verbose)

def load_straingauge_data(telem_files):
    """
    Load straingauge data from a list of telemetry files

    :param telem_files: List of telemetry files
    :type telem_files: list
    :return: Strain gauge data
    :rtype: numpy.ndarray
    """

    straingauge_filename_loc = np.where(np.array([file if 'straingauge' in str(file) else None for file in telem_files]) != None)[0][0]
    return np.load(telem_files[straingauge_filename_loc])

def load_dtt_data(telem_files):
    """
    Load downtiptilt data from a list of telemetry files

    :param telem_files: List of telemetry files
    :type telem_files: list
    :return: Downtip tilt data
    :rtype: numpy.ndarray
    """

    dtt_loc = np.where(np.array([file if 'downtiptilt' in str(file) else None for file in telem_files]) != None)[0][0]
    downtiptilt = np.load(telem_files[dtt_loc])
    return downtiptilt

def has_ocam2k_data(telem_files):
    """
    Check if a list of telemetry files contains ocam2k data
    
    :param telem_files: List of telemetry files (strings or Path objects)
    :type telem_files: list
    :return: True if the list contains a file with 'ocam2k' in its name/path, False otherwise
    :rtype: bool
    """
    
    return any('ocam2k' in str(f) for f in telem_files)

def load_ocam2k_data(telem_files):
    """
    Load ocam2k data from a list of telemetry files.
    
    Searches the provided list for the first file containing 'ocam2k' 
    in its name/path and loads it using numpy.

    :param telem_files: List of telemetry files (strings or Path objects)
    :type telem_files: list
    :return: Ocam2k data
    :rtype: numpy.lib.npyio.NpzFile from ocam2k telemetry
    :raises ValueError: If the input list is empty or the file cannot be loaded.
    :raises FileNotFoundError: If no matching file is found in the list, 
                               or if the file does not exist on disk.
    """
    
    # 1. Check if the list is empty or None
    if not telem_files:
        raise ValueError("The provided telem_files list is empty or invalid.")

    # 2. Find the first file matching 'ocam2k'
    # next() grabs the first match. If no match is found, it returns None.
    ocam2k_file = next((f for f in telem_files if 'ocam2k' in str(f)), None)

    # Handle the case where the string isn't found in any filenames
    if ocam2k_file is None:
        raise FileNotFoundError(
            "Could not find any file containing 'ocam2k' in the provided telemetry list."
        )

    # 3. Safely attempt to load the numpy array
    try:
        ocam2k = np.load(ocam2k_file)
        return ocam2k
        
    except FileNotFoundError:
        # The file was in the string list, but doesn't actually exist on the hard drive
        raise FileNotFoundError(
            f"The file '{ocam2k_file}' was in the list, but does not exist on disk."
        )
    except (OSError, ValueError) as e:
        # The file exists, but isn't a valid .npy or .npz file, or is corrupted
        raise ValueError(
            f"Failed to load numpy data from '{ocam2k_file}'. Original error: {e}"
        )

def has_four_lgs_data(telem_file):
    """
    Check if a numpy.lib.npyio.NpzFile file contains data for all four LGS

    :param telem_file: The numpy.lib.npyio.NpzFile to check
    :type telem_file: numpy.lib.npyio.NpzFile
    :return: True if the file contains data for all four LGS, False otherwise
    :rtype: bool
    """

    return all(get_ocam2k_intensity_key(i) in telem_file.files for i in range(1, 5))

def get_ocam2k_intensity_key(i):
    """
    Get the key corresponding to WFS i in the ocam2k telemetry data. For example,
    for i=1, return "SHWFS1SubAper
    
    :param i: The WFS number, from 1 to 4
    :type i: int
    :return: The key corresponding to WFS i in the ocam2k telemetry data
    :rtype: str
    """

    assert(i >= 1 and i <= 4), "WFS number must be between 1 and 4"
    return f"SHWFS{i}SubAper"


def load_straprtc_data(telem_files):
    """ Load straprtc data from a list of telemetry files

    :param telem_files: List of telemetry files
    :type telem_files: list
    :return: Straprtc data
    :rtype: numpy.ndarray
    """

    strap_loc = np.where(np.array([file if 'strap' in str(file) else None for file in telem_files]) != None)[0][0]
    straprtc = np.load(telem_files[strap_loc])
    return straprtc

def load_utt_data(telem_files):
    """ Load utt data from a list of telemetry files

    :param telem_files: List of telemetry files
    :type telem_files: list
    :return: Utt data
    :rtype: numpy.ndarray
    """

    utt_loc = np.where(np.array([file if 'uptiptilt' in str(file) else None for file in telem_files]) != None)[0][0]
    uptiptilt = np.load(telem_files[utt_loc])
    return uptiptilt

def get_xinetics_data(telem_files):
    """ Load xinetics data from a list of telemetry files

    :param telem_files: List of telemetry files
    :type telem_files: list
    :return: Xinetics data
    :rtype: numpy.ndarray
    """

    xinetics_loc = np.where(np.array([file if 'xinetics' in str(file) else None for file in telem_files]) != None)[0][0]
    xinetics = np.load(telem_files[xinetics_loc])
    return xinetics

### Plotting helper functions ###

def get_wfs_index_map():
    """
    Maps from a plotting index in a 2x2 grid of axes to the corresponding KAPA LGS
    WFS identifier in the telemetry data. For example the top left subplot (index 0)
    corresponds to shwfs4 in the telemetry data.

    The 2x2 grid of axes from plt.subplots(2, 2)is arranged like this:
    
    0 | 1
    ------
    2 | 3

    if you use ax = axs[i // 2, i % 2] feeding in i from 0 to 3.

    Note that this mapping was derived by ensuring fratricide spikes point to
    the center of the plot when the intensity from each WFS is plotted according
    to the reference mapping from populate_subap_map_with_data

    :return: A dictionary mapping from plotting index to WFS identifier (int to int)
    :rtype: dict
    """
    
    return {0: 4,
            1: 1,
            2: 3,
            3: 2}

def populate_subap_map_with_data(frame_data, subap_map=None):
    """
    Put 304 data points from a single frame of Keck I OSIRIS / KAPA
    ocam2k telemetry intensity data into a 2D map of the subapertures
    for plotting. Convention derived from "Subaperture conventions.docx", copied to
    https://docs.google.com/document/d/1B75WPM39Olx3Qv3C_gIpOyOSGmmICqj5/edit?rtpof=true&tab=t.0

    :param frame_data: The 304 data points from a single frame of telemetry data
    :type frame_data: numpy.ndarray or list of length 304
    :param subap_map: The 2D map of the subapertures, with shape (20, 20) and values of 1 for subaperture pixels and NaN for non-subaperture pixels. If None, the function will load the subap map from the default path.
    :type subap_map: numpy.ndarray
    """

    # Check params
    assert(len(frame_data) == 304), f"frame_data must be a list or numpy array of length 304, got {len(frame_data)}"
    if subap_map is None:
        subap_map = load_subap_map()
    assert(subap_map.shape == (20, 20)), f"subap_map must have shape (20, 20), got {subap_map.shape}"
    assert((subap_map == 1).sum() == 304), f"subap_map must have 304 pixels with value 1, got {(subap_map == 1).sum()}"

    result = np.zeros_like(subap_map, dtype=float)
    result[subap_map == 1] = frame_data
    return result

def plot_data_on_KAPA_WFSs(data, norm=None, cmap=None, cbar_label=None, axes_label="Subaperture"):
    """
    Plot data from the 4 KAPA LGS WFSs. The data
    should be organized from WFS1 to 4. If passed in, the normalization norm and 
    colorbar cmap are shared across all 4 subplots

    :param data: A list of 4 numpy arrays, each with shape (x, y), or numpy array of size (4, x, y), and values for unused display pixels (real pixels or subapertures) having NaN values
    :type data: list of numpy.ndarray
    :param norm: The normalization to use for the colorbar
    :type norm: matplotlib.colors.Normalize
    :param cmap: The colormap to use for the plots
    :type cmap: matplotlib.colors.Colormap
    :param cbar_label: The label for the colorbar (what the values of data are measuring)
    :type cbar_label: str
    """

    # Verify input
    if isinstance(data, list):
        assert len(data) == 4, "data must be a list of 4 numpy arrays"
    elif isinstance(data, np.ndarray):
        assert data.shape[0] == 4, "data must be a numpy array of size (4, x, y)"
    else:
        raise TypeError("data must be a list of 4 numpy arrays or a numpy array of size (4, x, y)")
    assert isinstance(norm, (Normalize, LogNorm, type(None))), "norm must be an instance of matplotlib.colors.Normalize, matplotlib.colors.LogNorm, or None"
    assert isinstance(cmap, (Colormap, type(None))), "cmap must be an instance of matplotlib.colors.Colormap or None"

    # Set up the plot
    fig, axs = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    
    # Get mapping
    sensor_map = get_wfs_index_map()

    # Normalize size of data. Should be (4, x, y) after this step
    if isinstance(data, list):
        data = np.array(data)

    # Setup default colormap
    if cmap == None:
            cmap = plt.cm.viridis.copy()
            cmap.set_bad(color="black")

    # Set default Norm
    if norm == None:
        norm = Normalize(vmin=np.nanmin(data), vmax=np.nanmax(data))\

    # Put data on the plot
    for i in range(4):
        ax = axs[i // 2, i % 2]
        im = ax.imshow(data[sensor_map[i] - 1], norm=norm, cmap=cmap)
        ax.set_title(f"SHWFS {sensor_map[i]}")
        ax.set_xlabel(axes_label + " x index")
        ax.set_ylabel(axes_label + " y index")
        ax.set_xticks(range(0, data.shape[1], int(data.shape[1] / 10)))
        ax.set_yticks(range(0, data.shape[2], int(data.shape[2] / 10)))

    # 4. Add the colorbar using the 'axs' grid
    # By passing the whole array 'axs', it places the colorbar to the right of the entire grid
    cbar = fig.colorbar(im, ax=axs, shrink=0.7, aspect=30, pad=0.02)
    cbar.set_label(cbar_label, rotation=270, labelpad=15)

    plt.show()

def load_ocam2k_wfs_data(ocam2k, ocam2k_index, units=u.adu):
    """
    Load the ocam2k data for each of the 4 WFSs.

    :param ocam2k: A dictionary containing the ocam2k data for each WFS
    :type ocam2k: dict
    :param ocam2k_index: The index of the data to load
    :type ocam2k_index: int
    :param units: The units to convert the data to of u.adu, u.electron, or u.electron / u.second
    :type units: astropy.units.Unit
    :return: A list of 4 numpy arrays, each with shape (304,) containing the subaperture values in the specified units
    :rtype: list of numpy.ndarray
    """
    
    # get data
    all_data = [] # (4, 304) array of subaperture values as units
    for i in range(1, 5):
        adu_data = ocam2k[get_ocam2k_intensity_key(i)][ocam2k_index] * u.adu
        
        # Get the values to append in the right units
        if units == u.electron:
            to_append = ku.convert_adu_to_photo_electrons(adu_data, Time(ocam2k['timestamp'][0]), "4LGS").value
        elif units == u.electron / u.second:
            frame_differences = np.diff(ocam2k['rawtimestamp'] * u.ns)
            frame_differences = frame_differences[frame_differences > 0 * u.ns] # Remove the negative differences (which are likely due to resets)
            implied_frequencies = (1 / frame_differences).to(u.Hz)
            frame_rate = implied_frequencies.mean()
            to_append = ku.convert_adu_to_flux(adu_data, frame_rate, Time(ocam2k['timestamp'][0]), "4LGS").value
        elif units == u.adu:
            to_append = adu_data.value
        else:
            raise ValueError("Unsupported units")
        all_data.append(to_append)
    
    return all_data

def plot_ocam2k_data(ocam2k, ocam2k_index, Norm=LogNorm, units=u.adu):
    """
    Plot the ocam2k data for each of the 4 WFSs.

    :param ocam2k: A dictionary containing the ocam2k data for each WFS
    :type ocam2k: dict
    :param ocam2k_index: The index of the data to plot
    :type ocam2k_index: int
    :param Norm: The normalization to use for the colorbar
    :type Norm: matplotlib.colors.Normalize
    :param units: The units to convert the data to of u.adu, u.electron, or u.electron / u.second
    :type units: astropy.units.Unit
    """
    
    all_data = load_ocam2k_wfs_data(ocam2k, ocam2k_index, units=units)
    all_values = np.concatenate([d.flatten() for d in all_data])

    global_vmin = np.nanpercentile(all_values, 5)
    global_vmax = np.nanmax(all_values)
    shared_norm = Norm(vmin=global_vmin, vmax=global_vmax)

    # Put data onto subapertures
    data = []
    for i in range(4):
        this_wfs_data = populate_subap_map_with_data(all_data[i])
        this_wfs_data = np.where(this_wfs_data == 0, np.nan, this_wfs_data)
        data.append(this_wfs_data)

    # Plot data
    plot_data_on_KAPA_WFSs(data, norm=shared_norm, cmap=None, cbar_label=f"{units}")


def read_ocam2k_background(path):
    """
    Read in the ocam2k background data from the given path. The background data is stored as a binary file with big-endian unsigned 16-bit integers. The first 6400 pixels (80x80) are the non-zero background values for the subapertures, and the rest of the pixels are zero.

    :param path: The path to the ocam2k background data file
    :type path: str or pathlib.Path
    :return: A numpy array of shape (80, 80) containing the background values for the subapertures in electrons
    :rtype: numpy.ndarray
    """

    # Check inputs
    assert(isinstance(path, (str, Path))), "path must be a string or pathlib.Path"
    assert(Path(path).exists()), f"File not found: {path}"

    # Read in raw data 228x228 binary file of big-endian unsigned 16-bit integers
    # '>' means Big-Endian, 'u2' means unsigned 2-byte (16-bit)
    raw_data = np.fromfile(path, dtype='>u2')
    square_side_length = int(np.sqrt(raw_data.size))
    assert(square_side_length == 228)

    # Check that the first 6400 pixels are non-zero and the rest are zero
    # corresponds to the 80x80 pixels
    assert(np.sum(raw_data != 0) == 6400) # Only 6400 non-zero pixels
    assert(np.where(raw_data.flatten() == 0)[0][0] == 6400) # The first 6400 pixels are non-zero, the rest are zero
    
    # Let's get the non-zero pixels
    pupil_data = raw_data[raw_data != 0].reshape((80, 80))

    # Convert them to the right units
    pupil_data_adu = pupil_data * u.adu
    pupil_data_electrons = ku.convert_adu_to_photo_electrons(pupil_data_adu, Time('2026-03-04T10:46:14'), '4LGS')
    return pupil_data_electrons

### Data Processing Helper Functions ###

def compute_median_values(data, thresh=20):
    """
    Compute median lit subaperture and unlit subaperture values across
    all 304 subapertures on a single WFS. The cutoff between the lit and unlit
    subpopulations are determined by the thresh parameter which defines a 
    percentile above and below which to separate the data. Values here could be
    intensities, electron counts, or standard deviation of electron counts, etc.
    Just some data defined across all subapertures which are defined into a 
    high and low population
    
    :param data: A numpy array or list of subaperture intensities for a single WFS, with length 304
    :param thresh: The percentile threshold for separating lit and unlit subapertures
    :return: A tuple of (median_unlit_intensity, median_lit_intensity)
    """
    
    # Verify inputs
    assert(isinstance(data, (np.ndarray, list))), "data must be a numpy array or list"
    assert(len(data) == 304), f"data must have length 304, got {len(data)}"
    assert(isinstance(thresh, (int, float))), "thresh must be an int or float"
    assert(thresh >= 0 and thresh <= 100), "thresh must be a percentile between 0 and 100"

    # Split the data into lit and unlit subapertures based on the threshold percentile
    thresh = np.percentile(data, thresh)
    lit_data = data[data > thresh]
    unlit_data = data[data <= thresh]

    # Get the medians
    median_lit_intensity = np.median(lit_data)
    median_unlit_intensity = np.median(unlit_data)

    # Return final values
    return median_unlit_intensity, median_lit_intensity

def compute_median_values_for_all_wfs(data, thresh=20):
    """
    Compute median lit subaperture and unlit subaperture values across
    all 304 subapertures on all 4 WFSs. The cutoff between the lit and unlit
    subpopulations are determined by the thresh parameter which defines a 
    percentile above and below which to separate the data. Values here could be
    intensities, electron counts, or standard deviation of electron counts, etc.
    Just some data defined across all subapertures which are defined into a 
    high and low population

    :param data: A list of numpy arrays, each containing the data for a single WFS with length 304, or a numpy array of shape (4, 304) where the first dimension corresponds to the WFS number and the second dimension corresponds to the subapertures. 
    :param thresh: The percentile threshold for separating lit and unlit subapertures
    :type thresh: int or float
    :return: A tuple of (median_unlit_values, median_lit_values) where each is a list of length 4 containing the median unlit and lit values for each WFS respectively
    :rtype: tuple of (list, list)

    """
    
    # Check params
    if isinstance(data, list):
        assert len(data) == 4, "data must be a list of 4 numpy arrays"
        for d in data:
            assert len(d) == 304, "each array in data must have length 304"
    elif isinstance(data, np.ndarray):
        assert data.shape == (4, 304), "data must be a numpy array of shape (4, 304)"
    else:
        raise TypeError("data must be a list of 4 numpy arrays or a numpy array of shape (4, 304)")
    assert isinstance(thresh, (int, float)), "thresh must be an int or float"
    assert thresh >= 0 and thresh <= 100, "thresh must be a percentile between 0 and 100"

    # Loop over all 4 WFSs and compute median values for each
    median_unlit_values = []
    median_lit_values = []    
    for wfs_number in range(1, 5):
        wfs_index = wfs_number - 1
        median_unlit_val, median_lit_val = compute_median_values(data[wfs_index], thresh)
        median_unlit_values.append(median_unlit_val)
        median_lit_values.append(median_lit_val)

    return median_unlit_values, median_lit_values

def compute_aperture_wise_electron_stats(ocam2k, hdr_tbl):
    """
    Compute the mean and standard deviation of the number of electrons read per frame
    for each subaperture across all 4 WFSs. The mean and standard deviation are
    returned as two numpy arrays with shape (4, 304) and units of electrons

    :param ocam2k: The ocam2k telemetry data loaded from the .npy file, as a numpy.lib.npyio.NpzFile
    :type ocam2k: numpy.lib.npyio.NpzFile
    :param hdr_tbl: The header table loaded from the observed image fits file, as an astropy Table
    :type hdr_tbl: astropy.table.Table
    :return: A tuple of (sensor_mean_electrons, sensor_stds_electrons) where each is a numpy array with shape (4, 304) and units of electrons per read of the mean and standard deviation respectively
    :rtype: tuple of (numpy.ndarray, numpy.ndarray)
    """

    # Verify inputs
    assert(isinstance(ocam2k, np.lib.npyio.NpzFile)), "ocam2k must be a dictionary"
    for i in range(1, 5):
        key = get_ocam2k_intensity_key(i)
        assert(key in ocam2k), f"ocam2k must contain key {key}"
        assert(isinstance(ocam2k[key], np.ndarray)), f"ocam2k[{key}] must be a numpy array"
        assert(ocam2k[key].ndim == 2), f"ocam2k[{key}] must be a 2D numpy array with shape (num_frames, 304), got shape {ocam2k[key].shape}"
        assert(ocam2k[key].shape[1] == 304), f"ocam2k[{key}] must have shape (num_frames, 304), got shape {ocam2k[key].shape}"
    assert(isinstance(hdr_tbl, (Table, Row))), "hdr_tbl must be an astropy Table or Row"
    assert('lgs_wfs_rate' in hdr_tbl.columns), "hdr_tbl must contain column 'lgs_wfs_rate', the frame rate of the camera"
    assert('t_exposure_start' in hdr_tbl.columns), "hdr_tbl must contain column 't_exposure_start', the start time of the exposure"

    # Get date of observation
    if isinstance(hdr_tbl, Table):
        exposure_start = hdr_tbl["t_exposure_start"][0]
    else:
        exposure_start = hdr_tbl["t_exposure_start"]

    # Compute mean and standard deviation of electrons for each subaperture across all 4 WFSs
    sensor_mean_electrons = []
    sensor_stds_electrons = []
    for i in range(1, 5):
        # Get adus, fluxes, and electrons for this wfs
        adus = ocam2k[get_ocam2k_intensity_key(i)] * u.adu
        electrons = ku.convert_adu_to_photo_electrons(adus, date=Time(exposure_start), mode="KAPA")
        sensor_mean_electrons.append(electrons.mean(axis=0))
        sensor_stds_electrons.append(electrons.std(axis=0))
    
    # Convert lists to numpy arrays with the right units
    sensor_mean_electrons = np.array(sensor_mean_electrons) * u.electron
    sensor_stds_electrons = np.array(sensor_stds_electrons) * u.electron

    # Return final values
    return sensor_mean_electrons, sensor_stds_electrons