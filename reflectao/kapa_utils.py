"""
Utilities for working with data from the KAPA upgrade on Keck I LGS for OSIRIS.

Noah Stiegler
4/15/26

"""

import numpy as np
from astropy import units as u
from astropy import time

# List of nights with telemetry available
kapa_night_list = ["2025nov06",
                   "2025nov08",
                   "2025dec04",
                   "2025dec06",
                   "2026feb26",
                   "2026feb28",
                   "2026jan12",
                   "2026jan31",
                   "2026mar04"]

def get_wfso1_gain(date, mode):
    """
    Get the WFSO1 gain for a given mode and date. This is the gain of the OCAM2k detector, which is how many electrons cascade out of a single photoelectron. This is not the same as the system gain, which is how many electrons correspond to one ADU (count) on the detector.

    :param date: The date of the observation as an astropy DateTime object
    :type date: astropy.time.Time
    :param mode: The mode of the WFSO1, either 'LGS' or '4LGS' (KAPA)
    :type mode: str
    :return: The WFSO1 gain for the given mode and date (unitless)
    :rtype: float
    """

    # ~Date of KAPA installation
    # Later you can put the dates of any changes here and change the
    # function to return the right gain for the given detector epoch
    KAPA_INSTALL_DATE = time.Time('2024-11-01') 
    assert(date > KAPA_INSTALL_DATE), f"Date must be after {KAPA_INSTALL_DATE.iso} for KAPA gain values to apply"

    # From Avinash, see https://docs.google.com/document/d/1yU_FpeVar_OlH-0ezAmpHIG7s87cANQaT-MNuvUhbBk/edit?tab=t.0
    # and from KAPA commissioning file header /g3/data/KECK/20260112_OSIRIS/KOA_13410/OSIRIS/raw/i260112_a000070.fits
    # Also in the O1SMGN spot on the image FITS header
    if mode == 'LGS':
        return 200 # Gain for single lasers (LGS)
    elif mode == '4LGS' or mode == 'KAPA':
        return 600 # Gain for multiple lasers (KAPA)
    else:
        raise ValueError(f"Invalid mode '{mode}'. Valid modes are currently: 'LGS' and '4LGS' AKA 'KAPA'")

def get_system_gain(date):
    """
    Get the system gain for the Keck I OCAM2k LGS WFS detector.
    This is how many electrons correspond to one ADU (count) on the detector.
    This is a constant value that does not depend on mode

    :param date: The date of the observation as an astropy DateTime object
    :type date: astropy.time.Time
    :return: The system gain for the Keck I OCAM2k detector
    :rtype: astropy.units.Quantity with units of electrons / ADU
    """

    # ~Date of KAPA installation
    # Later you can put the dates of any changes here and change the
    # function to return the right gain for the given detector epoch
    KAPA_INSTALL_DATE = time.Time('2024-11-01') 
    assert(date > KAPA_INSTALL_DATE), f"Date must be after {KAPA_INSTALL_DATE.iso} for KAPA gain values to apply"

    # From Avinash, see https://docs.google.com/document/d/1yU_FpeVar_OlH-0ezAmpHIG7s87cANQaT-MNuvUhbBk/edit?tab=t.0
    # and from KAPA commissioning file header /g3/data/KECK/20260112_OSIRIS/KOA_13410/OSIRIS/raw/i260112_a000070.fits
    return 28 * u.electron / u.adu # e- / ADU, how many electrons correspond to one ADU (count). Constant

def convert_adu_to_photo_electrons(adu_per_subap, date, mode):
    """
    Description
        Convert Keck I OCAM2k subaperture counts from ADU to photoelectrons
        recorded by the detector in that 2x2 subaperture (probably works per pixel too)

        :param adu_per_subap: Intensities recorded in each subaperture, in ADU
        :type adu_per_subap: numpy.ndarray of astropy.units.Quantity with units of ADU
        :param date: The date of the observation as an astropy DateTime object
        :type date: astropy.time.Time
        :param mode: The mode of the WFSO1, either 'LGS' or '4LGS' (KAPA)
        :type mode: str
        :return: Intensities in each subaperture converted to photoelectrons
        :rtype: numpy.ndarray of astropy.units.Quantity with units of electrons
    """

    assert(type(adu_per_subap) == u.Quantity and adu_per_subap.unit == u.adu), "Input must be a numpy array of astropy Quantity with units of ADU"

    wfso1_gain = get_wfso1_gain(date, mode) # Unitless, how many e- cascade out of a single photoelectron
    system_gain = get_system_gain(date) # e- / ADU, how many electrons correspond to one ADU (count)
    return adu_per_subap * (system_gain / wfso1_gain)

def convert_photo_electrons_to_adu(photo_electrons_per_subap, date, mode):
    """
    Description
        Convert from photoelectrons to ADU recordedby the  Keck I OCAM2k

        :param photo_electrons_per_subap: Intensities recorded in each subaperture, in photoelectrons
        :type photo_electrons_per_subap: numpy.ndarray of astropy.units.Quantity with units of electrons
        :param date: The date of the observation as an astropy DateTime object
        :type date: astropy.time.Time
        :param mode: The mode of the WFSO1, either 'LGS' or '4LGS' (KAPA)
        :type mode: str
        :return: Intensities in each subaperture converted to ADU
        :rtype: numpy.ndarray of astropy.units.Quantity with units of ADU
    """

    assert(type(photo_electrons_per_subap) == u.Quantity and photo_electrons_per_subap.unit == u.electron), "Input must be a numpy array of astropy Quantity with units of electrons"

    wfso1_gain = get_wfso1_gain(date, mode) # Unitless, how many e- cascade out of a single photoelectron
    system_gain = get_system_gain(date) # e- / ADU, how many electrons correspond to one ADU (count)
    return photo_electrons_per_subap * (wfso1_gain / system_gain)


def convert_photo_electrons_to_flux(photo_electrons_per_subap, frame_rate):
    """
    Description
        Convert from photoelectrons to flux (photoelectrons per second) recorded by the Keck I OCAM2k

        :param photo_electrons_per_subap: Intensities recorded in each subaperture, in photoelectrons
        :type photo_electrons_per_subap: numpy.ndarray of astropy.units.Quantity with units of electrons
        :param frame_rate: Frame rate of the WFS, in Hz
        :type frame_rate: astropy.units.Quantity with units of Hz
        :return: Intensities in each subaperture converted to flux (photoelectrons per second)
        :rtype: numpy.ndarray of astropy.units.Quantity with units of electrons / second
    """

    assert(type(photo_electrons_per_subap) == u.Quantity and photo_electrons_per_subap.unit == u.electron), "Input must be a numpy array of astropy Quantity with units of electrons"
    assert(type(frame_rate) == u.Quantity and frame_rate.unit == u.Hz), "Frame rate must be an astropy Quantity with units of Hz"

    # Brooke Digia: Flux = Gain * Counts / Exptime (https://mirametrics.com/help/mira_al_8/source/magnitude_calculations.htm)
    return photo_electrons_per_subap * frame_rate

def convert_adu_to_flux(adu_per_subap, frame_rate, date, mode):
    """
    Description
        Convert from ADU to flux (photoelectrons per second) recorded by the Keck I OCAM2k

        :param adu_per_subap: Intensities recorded in each subaperture, in ADU
        :type adu_per_subap: numpy.ndarray of astropy.units.Quantity with units of ADU
        :param frame_rate: Frame rate of the WFS, in Hz
        :type frame_rate: astropy.units.Quantity with units of Hz
        :param date: The date of the observation as an astropy DateTime object
        :type date: astropy.time.Time
        :param mode: The mode of the WFSO1, either 'LGS' or '4LGS' (KAPA)
        :type mode: str
        :return: Intensities in each subaperture converted to flux (photoelectrons per second)
        :rtype: numpy.ndarray of astropy.units.Quantity with units of electrons / second
    """

    assert(type(adu_per_subap) == u.Quantity and adu_per_subap.unit == u.adu), "Input must be a numpy array of astropy Quantity with units of ADU"
    assert(type(frame_rate) == u.Quantity and frame_rate.unit == u.Hz), "Frame rate must be an astropy Quantity with units of Hz"

    # Brooke Digia: Flux = Gain * Counts / Exptime (https://mirametrics.com/help/mira_al_8/source/magnitude_calculations.htm)
    photo_electrons_per_subap = convert_adu_to_photo_electrons(adu_per_subap, date, mode)
    return convert_photo_electrons_to_flux(photo_electrons_per_subap, frame_rate)