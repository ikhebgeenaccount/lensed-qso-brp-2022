from src.lensed_qso import LensedQSO


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
model spectra is then subtracted from the total
"""
from src.lens_subtraction.model_subtraction import model_subtraction
model_subtraction(lqso)


# AGNfitter

