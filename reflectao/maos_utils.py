import astropy.units as u
import numpy as np

def calc_psfgridsize(wavelength, evl_dx, sampling):
    """
    Calculate the MAOS evl.psfgridsize parameter based on the wavelength
    (evl.wvl), evl.dx (the spatial resolution of sampling on the mirror),
    and angular sampling rate we want on the PSF

    Calculated as 
    sampling = wavelength / (psfgridsize * evl_dx)

    So if we want to get the psfgridsize, we can rearrange the formula to:
    psfgridsize = wavelength / (sampling * evl_dx)

    :param wavelength: The wavelength of the light in meters
    :type wavelength: astropy.units.Quantity
    :param evl_dx: The spatial resolution of sampling on the mirror in meters
    :type evl_dx: astropy.units.Quantity in meters per pixel
    :param sampling: The sampling rate in radians or arcseconds
    :type sampling: astropy.units.Quantity
    :return: The MAOS psfgridsize parameter
    :rtype: int
    """

    sampling = sampling.to(u.radian)
    wavelength = wavelength.to(u.meter)
    evl_dx = evl_dx.to(u.meter)

    

    psfgridsize = wavelength / (sampling.value * evl_dx)

    # Make sure psfgridsize is a power of 2 (rounding up is safe)
    return int(2 ** np.ceil(np.log2(psfgridsize)))