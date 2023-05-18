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
# frequency of the demerit map using the given spectral index (default -0.7).
#
# Details of the demerit calculation can be found in Section 3 of Mauch et al.
# (2020).
#
# Tom Mauch
# May 2023

import argparse
import logging
import os

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from demerit import demeritlib

logging.basicConfig(level=logging.INFO)
log = logging.getLogger('demerit')


def create_parser():
    parser = argparse.ArgumentParser(description="Make an all-sky FITS demerit map.",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", "--fits-output", type=str,
                        help="Name of output FITS file. (default: './demerit_<band>.fits)")
    parser.add_argument("--pixel-scale", type=float, default=1.,
                        help="Pixel scale of the output FITS map in arcminutes. (default: 1 arcmin)")
    parser = demeritlib.common_options(parser)
    return parser


def make_fits_header(pix_scale, projection='MOL'):
    """ Generate a All-Sky Mollweide HDU object for a given pixel scale with empty data array.

    Parameters
    ----------
    pix_scale : float
        The desired pixel scale in degrees. Only square pixels supported.
    projection : string
        The projection to use for the FITS image. Anything accepted by wcslib
        works here, but really you will want an equal area all-sky projection.

    Returns
    -------
    hdu : :class:PrimaryHDU
        A fits header with data array that can be filled
    """

    w = WCS(naxis=2)
    w.wcs.ctype = [s + projection for s in ['RA---', 'DEC--']]
    w.wcs.cdelt = [-pix_scale, pix_scale]
    w.wcs.crval = [180., 0.]
    num_pixels_long = int(180. / pix_scale)
    num_pixels_lat = 2 * num_pixels_long
    w.wcs.crpix = [num_pixels_lat//2, num_pixels_long//2]
    data = np.empty((num_pixels_long, num_pixels_lat), dtype=np.float32)
    hdu = fits.PrimaryHDU(data=data, header=w.to_header())
    return hdu


def get_demerit_fits(fits_output, catalogue, flux, beamfwhm, radius, pointing_rms, gain_rms):

    out_array = fits_output.data
    w = WCS(fits_output)
    y_pix = np.arange(out_array.shape[1])
    for i, _ in enumerate(out_array):
        coords = w.array_index_to_world(i, y_pix)
        thisdec = w.array_index_to_world(i, out_array.shape[1]//2).dec.deg
        if np.isfinite(thisdec):
            log.info('Computing demerit for Declination: %5.1f', thisdec)
        out_array[i] = demeritlib.calculate_demerit_array(coords, catalogue, flux, beamfwhm,
                                                          radius, pointing_rms, gain_rms)
    fits_output.data[:] = out_array
    return fits_output


def main():
    parser = create_parser()
    args = parser.parse_args()
    if args.fits_output is None:
        args.fits_output = f'demerit_{args.band}.fits'
    fits_output = make_fits_header(args.pixel_scale / 60.)
    freq, beamfwhm = demeritlib.BAND_FREQ[args.band]
    log.info('Reading catalogue from: %s', os.path.basename(args.catalogue))
    catalogue, flux = demeritlib.load_catalogue(args.catalogue, freq)
    fits_output = get_demerit_fits(fits_output, catalogue, flux,
                                   beamfwhm, args.search_radius,
                                   args.pointing_rms,
                                   args.gain_rms)
    log.info('Writing output to %s', args.fits_output)
    fits_output.header['BANDCODE'] = args.band
    fits_output.writeto(args.fits_output, overwrite=True)


if __name__ == '__main__':
    main()
