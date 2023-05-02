#! /usr/bin/env python
#
# Generate a 'demerit' map in FITS format.
# 
# The demerit score is computed at each pixel in an output map by summing
# the individual demerit contributions of each source in the input
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


