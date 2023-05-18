#! /usr/bin/env python
#
# Derive the 'demerit' score at a position or a list of positions.
#
# The demerit score is computed at each input position by summing in
# quadrature the individual demerit contributions of each source in an input
# catalogue within a given radius. Suitable input catalogues
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
import textwrap
import os
import re

import numpy as np

from astropy import coordinates as c
from astropy import units as u

from demerit import demeritlib, static_dir


def create_parser():
    parser = argparse.ArgumentParser(description="Calculate demerit score at a position or list of positions",
                                     usage="Either: get_demerit.py [-h] [--catalog CATALOG] filename\n"
                                           "Or: get_demerit.py [-h] [--catalog CATALOG] RA DEC",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("position", type=str, nargs="+",
                        help=textwrap.dedent("""
                                             Either one or two arguments.

                                             One argument:
                                             A file containing a list of positions (two columns,
                                             one position per row).

                                             Two arguments:
                                             A single position (Right Ascension, Declination).

                                             The demerit score will be computed at each position.
                                             Positions can be given either as two floats in
                                             decimal degrees, or two floats with units appended
                                             (eg. 10h 20d), or in sexagesimal RA, Dec (hms dms).
                                             Sexagesimal subdivisions can be separated by ':', by
                                             spaces, or by units (eg. XXhXXmXX.Xs XXdXXmXX.Xs).
                                             Comment lines starting with '#' are ignored.'
                                             """))
    parser.add_argument("-n", "--num-sources", type=int, default=3,
                        help="Number of top-scoring demerit sources to print at each position.\n"
                             "(default: %(default)s)")
    parser = demeritlib.common_options(parser)
    return parser


def _ra_to_deg(ras):
    """Convert RA-like coordinates to decimal degrees array"""
    ra_deg = np.empty(len(ras), dtype=np.float)
    for i, ra in enumerate(ras):
        unit = None
        # String is just degrees without units
        try:
            float(ra)
        except ValueError:
            pass
        else:
            unit = u.deg
        # Check if we need to specify hourangle units
        if ' ' in ra or ':' in ra:
            unit = u.hourangle
        ra_deg[i] = c.Angle(ra, unit=unit).deg

    return ra_deg * u.deg


def _dec_to_deg(decs):
    """Convert Dec-like coordinates to decimal degrees array"""
    dec_deg = np.empty(len(decs), dtype=np.float)
    for i, dec in enumerate(decs):
        unit = None
        # String is just degrees without units
        try:
            float(dec)
        except ValueError:
            pass
        else:
            unit = u.deg
        # Check if we need to specify dms units
        if ' ' in dec or ':' in dec:
            unit = u.deg
        dec_deg[i] = c.Angle(dec, unit=unit).deg

    return dec_deg * u.deg


def main():
    parser = create_parser()
    args, leftovers = parser.parse_known_args()
    # Negative declination strings can cause upsets so just look for
    # negative sign and a number and call that a declination.
    for arg in leftovers:
        if re.compile('-[0-9]').match(arg) is not None:
            args.position.append(arg)
    if len(args.position) > 2:
        raise ValueError(f"Only 1 or 2 arguments accepted. See '{parser.prog} -h' for details.")
    if len(args.position) == 1:
        ras, decs = np.loadtxt(args.position[0], dtype=str, comments='#', unpack=True)
    else:
        ras, decs = [args.position[0]], [args.position[1]]

    input_coords = c.SkyCoord(_ra_to_deg(ras), _dec_to_deg(decs))

    frequency, beamfwhm = demeritlib.BAND_FREQ[args.band]

    catalogue, peak_flux = demeritlib.load_catalogue(args.catalogue, frequency)
    limit = beamfwhm * args.search_radius

    # Cumulative demerit to interpolate
    cum = np.load(os.path.join(static_dir, f'cum_demerit_{args.band}.npy'))

    catfile = os.path.basename(args.catalogue)
    print(f"Demerit score results at {args.band}-band:")

    for i, source in enumerate(input_coords):
        print(f"{'=' * 80}")
        print(f"{'Input Position' : ^20} {f'Demerit Score' : ^13} {f'Num. sources from {catfile}' : ^45}")
        print(f"{' ' *20} {'mJy / beam' : ^13} {f'within {limit : .1f} ({args.search_radius} x {(beamfwhm << u.arcmin) : <4.1f} FWHM)' : ^45}")
        print(f"{'=' * 80}")
        seps = catalogue.separation(source)
        keep = np.where(seps < limit)[0]
        seps = seps[keep]
        allatten = demeritlib.gaussian(seps, FWHM=beamfwhm)
        allflux = peak_flux[keep]
        allattenflux = allflux * allatten
        ds_squared = demeritlib.demerit_score_squared(beamfwhm, allattenflux, seps, args.pointing_rms, args.gain_rms)
        sort_args = np.argsort(ds_squared)[:-args.num_sources - 1:-1]
        this_d = np.sqrt(np.sum(ds_squared)) << u.mJy / u.beam
        this_sky_frac = np.interp(this_d.value, cum[0], cum[1])
        print(f"{source.to_string('hmsdms', precision=0)} {this_d.value : 8.1f} {len(keep) : 27d}")
        print(f"\n{this_sky_frac : .1f}% of sky has a lower demerit score.")
        print("\n"
              f"    The {len(sort_args)} sources of greatest demerit:\n"
              f"    {'-' * 70}\n"
              f"    {'Position' : ^20} {'Separation' : ^14} {'Amplitude' : ^14} {'Demerit' : ^14}\n"
              f"    {' ' * 20} {'arcmin' : ^14} {'mJy / beam' : ^14} {'mJy / beam' : ^14}\n"
              f"    {'-' * 70}")

        for arg in sort_args:
            this_d = np.sqrt(ds_squared[arg]) << u.mJy / u.beam
            max_flux = allflux[arg] << u.mJy / u.beam
            max_pos = catalogue[keep][arg].to_string('hmsdms', precision=0)
            max_sep = seps[arg]
            print(f"    {max_pos} {max_sep.arcmin : 10.1f} {max_flux.value : 14.1f} {this_d.value : 14.2f}")

        print(f"{'=' * 80}")


if __name__ == '__main__':
    main()
