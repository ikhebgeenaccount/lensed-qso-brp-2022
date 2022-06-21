from src.lensed_qso import LensedQSO
import pandas as pd


gal = 'J0806+2006'
"""
A LensedQSO uses a folder 'data/[gal]'. In this folder, it stores all files
pertaining to this galaxy, such as the SED in 'sed.csv'.
"""
lqso = LensedQSO(gal)


"""
SED construction
=============================================================
Several ways to construct an SED exist:
    - Manually, by adding fluxes, etc, to 'data/[gal]/sed.csv'
    - By converting magnitudes to fluxes
    - By reading a NED bar-separated table
"""

"""
Converting magnitudes to fluxes is done using filter specifications. These
specifications are stored in 'data/filter.csv'. If a filter is not defined there,
the magnitude cannot be converted, as it defines how it is converted.
"""
# Converting the magnitudes stored in 'mags.csv' to fluxes
from src.sed_compilation.mags_to_fluxes import mags_to_fluxes
mags_to_fluxes(lqso)

"""
Fluxes can also be split using magnitudes obtained from measurements for which
the used filter, the zeropoint or the used magnitude system is unknown.

This function automatically looks for the closest filter in the SED and splits
the total measured flux in that filter using the given magnitudes from SOURCE.
"""
from src.sed_compilation.mags_to_fluxes import mag_ratio_split_total_flux  # TODO: different name
mag_ratio_split_total_flux(lqso, ratio_source='SOURCE')
# SOURCE is the value in the 'source' column in mags.csv for the magnitude that is used to split a flux.

"""
The NASA Extragalactic Database (NED) is a very useful source for SEDs. After
searching for your target, the resulting photometry can be exported as a bar-
separated table. These bar-separated tables can be read and imported into the SED.

The table should be saved in 'data/[gal]/ned.txt'.
"""
from src.sed_compilation.ned_to_sed import ned_table_to_sed
ned_table_to_sed(lqso, ned_file='ned.txt')


"""
Lensing galaxy subtraction
=============================================================
After having constructed the SED of the lensed system, the lensing galaxy
can be subtracted. Model spectra from the Brown atlas (Brown et al. 2014?) are
fitted to the lensing galaxy SED (component G). An average of the best fitting
model spectra is then subtracted from the total.

"""
from src.lens_subtraction.model_subtraction import model_subtraction
from src.lens_subtraction.model_sed import fit

# These are the number of models used for the subtraction, they can be specified in properties.csv. Default=5
N = 5 if pd.isnull(lqso.props['no_models'].values[0]) else lqso.props['no_models'].values[0]

# This is the morphology of the foreground galaxy, as specified in properties.csv
morph = 'all' if pd.isnull(lqso.props.lens_type.values[0]) else lqso.props.lens_type.values[0]

# This function returns the averaged out fit as performed in model_sed
w, f, fe = fit(lqso, morph=morph, N=N)

# After the fit has been made, this function subtracts it from the SED in a new column
model_subtraction(lqso, w, f, fe)


"""
AGNFITTER
=============================================================
Sometimes new filters need to be added before running AGNfitter, there are
some functions that turn these into the proper files. Take the XML transmission profiles
for each filter from the csv filter profile service and save them to the folder FILTERPROFILES.
TODO: describe automated running of AGNfitter
"""

# Turns xml files into the proper format
from src.agnfitter.xml_to_txt import xml_to_txt
# xml_to_txt(Gemini_K.xml,Gemini_K.txt)


# It could be that the csv filter profile service lists transmission in percentage tp fraction, this function fixes it
from src.agnfitter.percent_to_fraction import percent_to_fraction
# percent_to_fraction(WIYN.U_HARRIS.txt, WIYN.U_HARRIS_fraction.txt)

# For e.g. radio wavelengths, no exact filter is known. Therefore we assume a tophat type filter. If it is given in angstrom or Ghz, 
# enter the central wavelength and the bandwidth. 
# If it is entered as an energy in keV, let the parameter central be the lower limit and bandwidth be the upper limit
# Output is always in angstrom
from src.agnfitter.tophat import tophat
# tophat(central,bandwidth,filtername)


