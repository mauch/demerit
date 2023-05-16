#! /usr/bin/env python
#
# Derive a plot of the cumulative demerit distribution given a
# FITS demerit map.
#
# The FITS demerit map should have been created using the 
# make_demerit_map.py script in this package. The script
# will produce png files containing the cumulative demerit
# distribution as well as a png file of the demerit map
# image itself.
#
# The script will also produce an npz file containing
# a spline fit to the cumulatve demerit distribution as is
# contained in the 'static' directory of this package, and
# read by the get_demerit.py script.
#
# Tom Mauch 
# May 2023
#
import argparse
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pylab as plt
from matplotlib.ticker import MaxNLocator


def create_parser():
    parser = argparse.ArgumentParser(description="Derive cumulative demerit plot from a demerit map.")
    parser.add_argument('demerit_map', type=str, help="Demerit FITS map filename. This file should have been \n"
                                                      "created using the make_demerit_map.py script.")
    return parser


def plot_cumulative_demerit(data, band):

    datamin = np.amin(data)
    datamax = np.amax(data)
    ndata = len(data)

    bins = np.logspace(np.log10(datamin), np.log10(datamax), 5000)
    # Demerit in mJy / beam
    hist, bins = np.histogram(data, bins=bins)

    cumhist = np.cumsum(hist)
    # Normalise cumhist - its maximum bin should be 100% of sky area in the image
    cumhist = cumhist / cumhist[-1] * 100.0

    # Point where cumhist is 95%
    plotmax = np.interp(95., cumhist, bins[1:])

    plt.figure(figsize=(5,2.5), tight_layout=True, dpi=300)
    plt.semilogy(bins[1:], cumhist)
    plt.xlim((datamin, plotmax))
    plt.ylim((0.1, 200))
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.text((datamin + plotmax) / 2, 1, f'{band}-band')
    plt.xlabel('Demerit Score (D / mJy)')
    plt.ylabel('% of sky < D')
    plt.savefig(f'cum_demerit_{band}.png')

    # Export the cumulative histogram
    output = np.stack((bins[1:], cumhist))
    np.save(f'cum_demerit_{band}', output)


parser = create_parser()
args = parser.parse_args()

ffhdu = fits.open(args.demerit_map)
# Assume data is Jy/beam and convert to mJy/beam
data = ffhdu[0].data[~np.isnan(ffhdu[0].data)] * 1000.
# Without BANDCODE the FITS file was likely not made with the right script anyhow.
band = ffhdu[0].header['BANDCODE']

plot_cumulative_demerit(data, band)
