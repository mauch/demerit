import logging
import os
import textwrap

from astropy.io import fits
from astropy import units as u
from astropy import coordinates as c
import numpy as np

from demerit import static_dir

log = logging.getLogger('demerit')

# Default MeerKAT rms dish pointing in degrees.
RMS_POINTING = 30. * u.arcsec
# Default MeerKAT rms gain fluctuation.
RMS_GAIN = 0.01
# Band to frequency and beam FWHM conversion
BAND_FREQ = {'UHF': [816. * u.MHz, 107./60. * u.deg],
             'L': [1284. * u.MHz, 68./60. * u.deg],
             'S': [2625. * u.MHz, 33.3/60. * u.deg]}


def common_options(parser):
    """Add arguments that are common to all of the scripts"""

    parser.add_argument("--catalogue", type=str, default=os.path.join(static_dir, "AllSky.fits.gz"),
                        help=textwrap.dedent("""\
                                             Path to a FITS file containing a binary table with
                                             a source catalogue used to compute demerit values in its
                                             FIRST extension.

                                             The binary table MUST contain 3 columns:
                                             1:'RA(2000)'   # The J2000 Right Ascension
                                             2:'DEC(2000)'  # The J2000 Declination
                                             3:'PEAK'       # Peak flux density in Jy/beam

                                             The extension must further contain a header
                                             keyword 'REF_FREQ', specifying the frequency of
                                             PEAK in Hz.
                                             """))
    parser.add_argument("-b", "--band", type=str, default="L", choices=["UHF", "L", "S"],
                        help="Demerit band. Must be one of either 'UHF', 'L', 'S'. "
                             "(default: %(default)s)")
    parser.add_argument("-r", "--search-radius", type=int, default=5,
                        help="Width around each position in multiples of the primary beam FWHM \n"
                             "to search for sources in the provided catalogue. (default: %(default)s)")
    parser.add_argument("--pointing-rms", type=lambda x: float(x) << u.arcsec, default=RMS_POINTING,
                        help="RMS antenna pointing in arcminutes. (default: %(default)s)")
    parser.add_argument("--gain-rms", type=float, default=RMS_GAIN,
                        help="RMS antenna gain. (default: %(default)s)")
    return parser


def gaussian(rho, FWHM=1. * u.deg):
    """Return the attenuation of a Gaussian with given FWHM for a given offset rho."""
    const = -4. * np.log(2)
    fwhmratio = rho / FWHM
    exponent = const * fwhmratio * fwhmratio
    return np.exp(exponent)


def demerit_score_squared(beamfwhm, allattenflux, seps, pointing_rms, gain_error):
    beamfwhm = beamfwhm << u.deg
    pconst = 8. * np.log(2.) / (beamfwhm * beamfwhm)
    pointing_rms = pointing_rms << beamfwhm.unit
    # Pointing flux error
    dsp = pconst * allattenflux * seps * pointing_rms
    # Gain error
    dsg = allattenflux * gain_error
    ds_squared = (dsp * dsp) + (dsg * dsg)
    return ds_squared


def calculate_demerit_array(positions, catalogue, flux, beamfwhm, radius,
                            pointing_rms=RMS_POINTING, gain_error=RMS_GAIN):
    """ Calculate the demerit score D for positions using sources from given catalogue"""

    radius = radius * beamfwhm
    # Trap nan positions
    demerit = np.full(positions.shape, np.nan, dtype=np.float32) * u.Jy / u.beam
    true_position_idx = np.where(np.isfinite(positions.ra) & np.isfinite(positions.dec))[0]
    true_positions = positions[true_position_idx]
    true_demerit = np.zeros(true_positions.shape, dtype=np.float32) * u.Jy / u.beam

    idx1, idx2, seps, _ = catalogue.search_around_sky(true_positions, radius)

    if len(idx1) == 0:
        return demerit

    # Get indices of targets with matches (these will have demerit>0)
    demerit_indices = np.unique(idx1)

    allatten = gaussian(seps, FWHM=beamfwhm)
    allflux = flux[idx2]
    allattenflux = allflux * allatten
    # Where does idx1 move to its next target
    idx1change = np.where(idx1[:-1] != idx1[1:])[0] + 1
    idx1change = np.insert(idx1change, 0, 0)
    ds_squared = demerit_score_squared(beamfwhm, allattenflux, seps, pointing_rms, gain_error)
    true_demerit[demerit_indices] = np.sqrt(np.add.reduceat(ds_squared, idx1change, dtype=np.float32))
    # Only fill demerit array with non-nan positions
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
    flux_scale = (freq / cat_freq)**-0.7
    flux = flux_scale * ffdata['PEAK'] * (u.Jy / u.beam)
    catalogue = c.SkyCoord(ffdata['RA(2000)'] * u.deg, ffdata['DEC(2000)'] * u.deg)

    return catalogue, flux
