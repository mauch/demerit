#! /usr/bin/env python
#
# Generate a 'demerit' map in FITS format.
# 
# The demerit score is computed at each pixel in an output map by summing in
# quadrature the individual demerit contributions of each source in the input
# catalogue within a given radius of the input pixel. Suitable input catalogues
# would be the NVSS or SUMSS as only brighter sources contribute appreciably
# to the demerit score at any given position.
#
# The included file AllSkyVZ.FIT contains an AIPS VZ table comprising the NVSS,
# SUMSS and MGPS-II catalogues which together cover the entire sky for sources
# with peak amplitude greater than 10 mJy/beam (extrapolated from the reference
# frequency of each respective survey to 1284 MHz using a spectral index of -0.7).
# This is read by default by the script and again extrapolated to the desired 
# frequency of the demerit map using the given default spectral index
# (default -0.7).
#
# Details of the demerit calculation can be found in Section 3 of Mauch et al.
# (2020).
# 
# Tom Mauch
# May 2023

import logging

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy import coordinates as c
import numpy as np

from matplotlib import pyplot as plt

import multiprocessing
import concurrent.futures


log = logging.getLogger('demerit')

def _gaussian(rho, FWHM=1. * u.deg):
    """Return the attenuation of a Gaussian with given FWHM for a given offset rho."""
    const = -4. * np.log(2)
    fwhmratio = rho / FWHM
    exponent = const * fwhmratio * fwhmratio
    return np.exp(exponent)


def make_fits_header(pix_scale):
    """ Generate a All-Sky Mollweide HDU object for a given pixel scale with empty data array.

    Parameters
    ----------
    pix_scale : float
        The desired pixel scale in degrees. Only square pixels supported.

    Returns
    -------
    hdu : :class:PrimaryHDU
        A fits header with data array that can be filled
    """

    w = WCS(naxis=2)
    w.wcs.ctype = [s + 'MOL' for s in ['RA---', 'DEC--']]
    w.wcs.cdelt = [-pix_scale, pix_scale]
    w.wcs.crval = [180., 0.]

    num_pixels_long = int(180. / pix_scale)

    num_pixels_lat = 2 * num_pixels_long

    w.wcs.crpix = [num_pixels_lat//2, num_pixels_long//2]

    d = np.empty((num_pixels_long, num_pixels_lat), dtype=np.float32)

    h = fits.PrimaryHDU(data=d, header=w.to_header())

    return h


def calculate_demerit(positions, catalogue, flux, beamfwhm, pointing_rms=0.5*u.deg, gain_error=0.01):
    """ Calculate the demerit score D for positions using sources from given catalogue"""


    pconst = 8. * np.log(2.) / (beamfwhm * beamfwhm)

    radius = 3. * beamfwhm


    # Trap nan positions
    demerit = np.full(positions.shape, np.nan, dtype=np.float32) * u.Jy / u.beam
    true_position_idx = np.where(np.isfinite(positions.ra) | np.isfinite(positions.dec))[0]
    true_positions = positions[true_position_idx]
    true_demerit = np.zeros(true_positions.shape, dtype=np.float32) * u.Jy / u.beam

    idx1, idx2, seps, _ = catalogue.search_around_sky(true_positions, radius)

    if len(idx1) == 0:
        return demerit

    # Get indices of targets with matches (these will have demerit>0)
    demerit_indices = np.unique(idx1)

    allatten = _gaussian(seps, FWHM=beamfwhm)
    allflux = flux[idx2]
    allattenflux = allflux * allatten
    # Where does idx1 move to its next target
    idx1change = np.where(idx1[:-1] != idx1[1:])[0] + 1
    idx1change = np.insert(idx1change, 0, 0)

    # Pointing flux error
    dsp = pconst * allattenflux * seps * pointing_rms
    # Gain error
    dsg = allattenflux * gain_error

    ds_squared = (dsp * dsp) + (dsg * dsg)

    true_demerit[demerit_indices] = np.sqrt(np.add.reduceat(ds_squared, idx1change, dtype=np.float32))

    demerit[true_position_idx] = true_demerit

    return demerit


def load_catalogue(filename, freq):
    """ Load an input catalogue (FITS binary table) into a SkyCoord array. """

    ffcat = fits.open(filename)

    ffdata = ffcat[1].data
    cat_freq = ffcat[1].header['REF_FREQ'] * u.Hz
    if type(freq) != u.Quantity:
        log.warn('Unknown unit for given frequency, assuming Hz')
    freq = freq << u.Hz
    flux_scale = (cat_freq / freq)**-0.7
    flux = flux_scale * ffdata['PEAK'] * (u.Jy / u.beam)
    catalogue = c.SkyCoord(ffdata['RA(2000)'] * u.deg, ffdata['DEC(2000)'] * u.deg)

    return catalogue, flux


def get_demerit(fits_output, catalogue, flux, beamfwhm):

    out_array = fits_output.data
    w = WCS(fits_output)
    y_pix = np.arange(out_array.shape[1])
    for i, _ in enumerate(out_array):
        coords = w.array_index_to_world(i, y_pix)
        thisdec = w.array_index_to_world(i, out_array.shape[1]//2).dec.deg
        print(thisdec)
        out_array[i] = calculate_demerit(coords, catalogue, flux, beamfwhm)
    fits_output.data[:] = out_array
    fits_output.writeto('test.fits', overwrite=True)
    return fits_output
