"""
Utilities for working with data from the Keck I telemetry system

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
import paarti.utils.maos_utils as mu

### Useful References ###

# List of nights with telemetry available
night_list = ["2025nov06",
              "2025nov08",
              "2025dec04",
              "2026dec06",
              "2026feb26",
              "2026feb28",
            #   "2026jan12", <-- No telemetry (TRS server down)
            #   "2026jan31", <-- Telemetry files are empty (lingering TRS issue)
              "2026mar04"]

# For each night, we want to know which sets and which images have good files to
# analyze. Here we'll organize that information by dictionaries where keys are
# night dates (as YYYYmonDD strings) and values are dictionaries with set_num
# keys and lists of their images with good telemetry files as values. We'll
# start by just listing all the nights and having empty dictionaries
good_img_telemetry = {
    # "2025nov06": {}, <--- set 1 image 78-173 exist. Not great logging. 
    "2025nov08": {},
    "2025dec04": {},
    "2026dec06": {},
    "2026feb26": {},
    "2026feb28": {},
    "2026mar04": {}
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



def get_fits_filename(night, set_num, image_num):
    """
    Get the fits filename for a given night, set, and image number
    
    :param night: The night of the observation, in the format 'YYYYmonDD', e.g. '2026mar04
    :type night: str'
    :param set: The set number of the observation
    :type set: int
    :param image: The image number of the observation
    :type image: int
    
    :return: The fits filename for the given night, set, and image number
    :rtype: str
    """

    _check_night_param(night)
    _check_set_param(set_num)
    _check_image_param(image_num)

    two_digit_year = night[2:4]
    two_digit_month = _convert_month(night[4:7])
    two_digit_day = night[7:9]
    return f"i{two_digit_year}{two_digit_month}{two_digit_day}_a{set_num:03d}{image_num:03d}.fits"

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
    return path_to_night / "raw" / fits_filename

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

def read_image_telemetry(night, set_num, image_num, verbose=False):
    """
    Read in the telemetry associated with a given image from the place
    it's stored on the MULab filesystem
    
    :param night: The night of the observation, in the format 'YYYYmonDD', e.g. '2026mar04
    :type night: str'
    :param set: The set number of the observation
    :type set: int
    :param image: The image number of the observation
    :type image: int
    
    """
    
    _check_night_param(night)
    _check_set_param(set_num)
    _check_image_param(image_num)

    path_to_image = get_path_to_image_telemetry(night, set_num, image_num)
    if verbose: print("Looking at telemetry from", path_to_image)

    # Get the files in the telemetry folder for this image
    telem_files = [file for file in path_to_image.glob("*")]
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

def load_ocam2k_data(telem_files):
    """
    Load ocam2k data from a list of telemetry files

    :param telem_files: List of telemetry files
    :type telem_files: list
    :return: Ocam2k data
    :rtype: numpy.ndarray
    """

    ocam2k_filename_loc = np.where(np.array([file if 'ocam2k' in str(file) else None for file in telem_files]) != None)[0][0]
    ocam2k = np.load(telem_files[ocam2k_filename_loc])
    return ocam2k

def get_ocam2k_key(i):
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

def populate_subap_map_with_data(frame_data, subap_map):
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
    assert(subap_map.shape == (20, 20)), f"subap_map must have shape (20, 20), got {subap_map.shape}"
    assert(sum(subap_map == 1) == 304), f"subap_map must have 304 pixels with value 1, got {sum(subap_map == 1)}"

    result = np.zeros_like(subap_map, dtype=float)
    result[subap_map == 1] = frame_data
    return result

def plot_KAPA_subapertures_with_data(data, norm=None, cmap=None, cbar_label=None):
    """
    Plot data from the 4 KAPA LGS WFSs in a 2x2 grid of subapertures. The data
    should be organized from WFS1 to 4. If passed in, the normalization norm and 
    colorbar cmap are shared across all 4 subplots

    :param data: A list of 4 numpy arrays, each with shape (20, 20), or numpy array of size (4, 20, 20), and values for the subapertures of the corresponding WFS with unused subapertures having NaN values
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
        for d in data:
            assert d.shape == (20, 20), "each array in data must have shape (20, 20)"
    elif isinstance(data, np.ndarray):
        assert data.shape == (4, 20, 20), "data must be a numpy array of size (4, 20, 20)"
    else:
        raise TypeError("data must be a list of 4 numpy arrays or a numpy array of size (4, 20, 20)")
    assert isinstance(norm, (Normalize, LogNorm, type(None))), "norm must be an instance of matplotlib.colors.Normalize, matplotlib.colors.LogNorm, or None"
    assert isinstance(cmap, (Colormap, type(None))), "cmap must be an instance of matplotlib.colors.Colormap or None"

    # Set up the plot
    fig, axs = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    
    # Get mapping
    sensor_map = get_wfs_index_map()
    subap_map = _load_subap_map()

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
        im = ax.imshow(data, norm=norm, cmap=cmap)
        ax.set_title(f"SHWFS {sensor_map[i]}")
        ax.set_xlabel("Subaperture x index")
        ax.set_ylabel("Subaperture y index")
        ax.set_xticks(range(0, subap_map.shape[1], 2))
        ax.set_yticks(range(0, subap_map.shape[0], 2))

    # 4. Add the colorbar using the 'axs' grid
    # By passing the whole array 'axs', it places the colorbar to the right of the entire grid
    cbar = fig.colorbar(im, ax=axs, shrink=0.7, aspect=30, pad=0.02)
    cbar.set_label(cbar_label, rotation=270, labelpad=15)

    plt.show()

### Data Processing Helper Functions ###

def compute_median_intensities(data, thresh=20):
    """
    Compute the median lit subaperture and unlit subaperture intensity across
    all 304 subapertures on a single WFS. The cutoff between the lit and unlit
    subpopulations are determined by the thresh parameter which defines a 
    percentile above and below which to separate the data.
    
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
    return median_unlit_intensity, 


def compute_aperture_wise_electron_stats(ocam2k, hdr_tbl):
    """
    Compute the mean and standard deviation of the number of electrons read per frame
    for each subaperture across all 4 WFSs. The mean and standard deviation are
    returned as two numpy arrays with shape (4, 304) and units of electrons
    """

    # Verify inputs
    assert(isinstance(ocam2k, dict)), "ocam2k must be a dictionary"
    for i in range(1, 5):
        key = get_ocam2k_key(i)
        assert(key in ocam2k), f"ocam2k must contain key {key}"
        assert(isinstance(ocam2k[key], np.ndarray)), f"ocam2k[{key}] must be a numpy array"
        assert(ocam2k[key].ndim == 3), f"ocam2k[{key}] must be a 3D numpy array with shape (num_frames, 20, 20)"
    assert(isinstance(hdr_tbl, pd.DataFrame)), "hdr_tbl must be a pandas DataFrame"
    assert('lgs_wfs_rate' in hdr_tbl.columns), "hdr_tbl must contain column 'lgs_wfs_rate', the frame rate of the camera"
    assert('t_exposure_start' in hdr_tbl.columns), "hdr_tbl must contain column 't_exposure_start', the start time of the exposure"

    # Compute mean and standard deviation of electrons for each subaperture across all 4 WFSs
    sensor_mean_electrons = []
    sensor_stds_electrons = []
    for i in range(1, 5):
        # Get adus, fluxes, and electrons for this wfs
        adus = ocam2k[get_ocam2k_key(i)] * u.adu
        electrons = ku.convert_adu_to_photo_electrons(adus, date=Time(hdr_tbl["t_exposure_start"][0]), mode="KAPA")
        sensor_mean_electrons.append(electrons.mean(axis=0))
        sensor_stds_electrons.append(electrons.std(axis=0))
    
    # Convert lists to numpy arrays with the right units
    sensor_mean_electrons = np.array(sensor_mean_electrons) * u.electron
    sensor_stds_electrons = np.array(sensor_stds_electrons) * u.electron

    # Return final values
    return sensor_mean_electrons, sensor_stds_electrons