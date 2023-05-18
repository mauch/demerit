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
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.axes import Axes
from mpl_toolkits.axes_grid1 import make_axes_locatable


def create_parser():
    parser = argparse.ArgumentParser(description="Derive cumulative demerit plot from a demerit map.")
    parser.add_argument('demerit_map', type=str, help="Demerit FITS map filename. This file should have been \n"
                                                      "created using the make_demerit_map.py script.")
    return parser


def get_cum(data, bins=5000):
    d = data[~np.isnan(data)]
    datamin = np.amin(d)
    datamax = np.amax(d)
    bins = np.logspace(np.log10(datamin), np.log10(datamax), 5000)
    # Demerit in mJy / beam
    hist, bins = np.histogram(d, bins=bins)
    cumhist = np.cumsum(hist)
    return bins, cumhist


def plot_cumulative_demerit(hdu):

    band = hdu.header['BANDCODE']
    data = hdu.data[~np.isnan(hdu.data)] * 1000.
    datamin = np.amin(data)

    bins, cumhist = get_cum(data)
    # Normalise cumhist - its maximum bin should be 100% of sky area in the image
    cumhist = cumhist / cumhist[-1] * 100.0

    # Point where cumhist is 95%
    plotmax = np.interp(95., cumhist, bins[1:])

    plt.figure(figsize=(5, 2.5), tight_layout=True, dpi=300)
    plt.semilogy(bins[1:], cumhist)
    plt.xlim((datamin, plotmax))
    plt.ylim((0.1, 200))
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.text((datamin + plotmax) / 2, 1, f'{band}-band')
    plt.xlabel('Demerit Score (D / mJy)')
    plt.ylabel('% of sky < D')
    plt.savefig(f'cum_demerit_{band}.png')
    plt.clf()
    plt.close('all')

    # Export the cumulative histogram
    output = np.stack((bins[1:], cumhist))
    np.save(f'cum_demerit_{band}', output)


def plot_image(hdu):
    dpi = 300
    width = 4000
    wcs = WCS(hdu)
    band = hdu.header['BANDCODE']
    data = hdu.data * 1000.
    bins, cumhist = get_cum(data)
    cumhist = cumhist / cumhist[-1] * 100.0
    vmax = np.interp(95., cumhist, bins[1:])
    vmin = np.interp(5., cumhist, bins[1:])
    image_height, image_width = data.shape
    finite_data = np.where(np.isfinite(data))
    if finite_data[0].size > 0:
        ymin = np.min(finite_data[0])
        ymax = np.max(finite_data[0])
        xmin = np.min(finite_data[1])
        xmax = np.max(finite_data[1])
    aspect = image_width / image_height
    height = width / aspect * 1.08
    fig = plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
    ax = fig.add_subplot(projection=wcs, frame_class=EllipticalFrame)
    ax.set_xlim(-0.5 + xmin, xmax + 0.5)
    ax.set_ylim(-0.5 + ymin, ymax + 0.5)
    im = ax.imshow(data, origin='lower', aspect='equal', vmin=vmin, vmax=vmax)
    ax.grid()
    ax.coords['ra'].set_ticklabel(color='white')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('bottom', pad='3%', size='5%', axes_class=Axes)
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label(label='Demerit (mJy / beam)', size=15)
    ax.set_title(f'{band}-band', size=20)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f'demerit_{band}.png')


def main():
    parser = create_parser()
    args = parser.parse_args()

    ffhdu = fits.open(args.demerit_map)

    plot_cumulative_demerit(ffhdu[0])
    plot_image(ffhdu[0])


if __name__ == '__main__':
    main()
