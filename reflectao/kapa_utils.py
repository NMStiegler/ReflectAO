"""
Utilities for working with data from the KAPA upgrade on Keck I LGS for OSIRIS.

Noah Stiegler
4/15/26

"""

import numpy as np
from astropy import units as u

def convert_adu_to_photo_electrons(adu_per_subap):
    """
    Description
        Convert Keck I OCAM2k subaperture counts from ADU to photoelectrons
        recorded by the detector in that 2x2 subaperture (probably works per pixel too)

        :param adu_per_subap: Intensities recorded in each subaperture, in ADU
        :type adu_per_subap: numpy.ndarray of astropy.units.Quantity with units of ADU
        :return: Intensities in each subaperture converted to photoelectrons
        :rtype: numpy.ndarray of astropy.units.Quantity with units of electrons
    """

    # From Avinash (see https://docs.google.com/document/d/1yU_FpeVar_OlH-0ezAmpHIG7s87cANQaT-MNuvUhbBk/edit?tab=t.0)
    # and from KAPA commissioning file header /g3/data/KECK/20260112_OSIRIS/KOA_13410/OSIRIS/raw/i260112_a000070.fits
    wfso1_gain = 600 # Unitless, how many e- cascade out of a single photoelectron. 600 for multiple lasers (KAPA), 200 for single lasers (LGS)
    system_gain = 28 * u.electron / u.adu # e- / ADU, how many electrons correspond to one ADU (count). Constant

    return adu_per_subap * (system_gain / wfso1_gain)

def convert_photo_electrons_to_adu(photo_electrons_per_subap):
    """
    Description
        Convert from photoelectrons to ADU recordedby the  Keck I OCAM2k

        :param photo_electrons_per_subap: Intensities recorded in each subaperture, in photoelectrons
        :type photo_electrons_per_subap: numpy.ndarray of astropy.units.Quantity with units of electrons
        :return: Intensities in each subaperture converted to ADU
        :rtype: numpy.ndarray of astropy.units.Quantity with units of ADU
    """

    # From Avinash (see https://docs.google.com/document/d/1yU_FpeVar_OlH-0ezAmpHIG7s87cANQaT-MNuvUhbBk/edit?tab=t.0)
    # and from KAPA commissioning file header /g3/data/KECK/20260112_OSIRIS/KOA_13410/OSIRIS/raw/i260112_a000070.fits
    wfso1_gain = 600 # Unitless, how many e- cascade out of a single photoelectron. 600 for multiple lasers (KAPA), 200 for single lasers (LGS)
    system_gain = 28 * u.electron / u.adu # e- / ADU, how many electrons correspond to one ADU (count). Constant

    return photo_electrons_per_subap * (wfso1_gain / system_gain)