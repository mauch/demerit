demerit
=======

MeerKAT demeirt map utilities

# Installation

Best to clone from github and use pip install from a setup virtual environment to install the scripts and their dependencies:

$ git clone https://github.com/mauch/demerit
$ pip install ./demerit

# Scripts

There are 3 scripts (in all cases use --help for further information):

1. make-demerit-map.py will generate a FITS file with an all-sky demerit map for a selected band. Beware - 
	this takes over 24 hours to run for an all-sky demerit map.

2. get-demerit.py will print information about the Demerit score and bright sources at a given input position and band.

3. plot-cum-demerit.py will take a FITS input from make-demerit-map.py and make png plots of the cumulative demerit score
	and the demerit map itself. It will also produce .npy files containt the cumulative distribution which can replace
	these contained in the static dir of the package. These npy files are used to derive the % of sky value reported by
	the get-demerit.py script.
