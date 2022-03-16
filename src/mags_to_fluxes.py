import os.path

import numpy as np
import pandas as pd

from src.lensed_qso import LensedQSO

from src.filters import FILTER_PROPERTIES, get_wavelength

import warnings


def mags_to_fluxes(lqso, components=None):
    """
    Reads the mags.csv file of galaxy and converts the magnitudes to fluxes, using the appropriate conversion formulas.
    Saves these fluxes to sed.csv.
    :param galaxy:
    :param sed_file:
    :param mags_file:
    :return: DataFrame with fluxes
    """
    if not hasattr(lqso, 'mags'):
        warnings.warn(lqso.name + ' has no mags file')
        return

    if components is None:
        components = ['_G', '_A', '_B', '_C', '_D', '']

    for i, row in lqso.mags.iterrows():

        try:
            # Find the wavelength of the filter
            wavelength = get_wavelength(row.telescope, row['filter'])
        except IndexError:
            print(f'No conversion found for {row.telescope} {row["filter"]} filter in data/filter.csv.')
            continue

        # Flag that is set to 1 if there is a G component, then it is assumed that all components are listed and
        # a total flux can be calculated
        # Set to 2 if there is a total magnitude
        total = 0

        # Used to store fluxes and flux_errs for all components
        fs = []
        fes = []

        f_fields = []
        fe_fields = []

        # Calculate fluxes and flux errors for each component
        for comp in components:
            # magnitude 0 means no magnitude given
            if row[f'mag{comp}'] == 0:
                continue

            f, fe = mag_to_flux(row.telescope, row['filter'], row[f'mag{comp}'], row[f'mag{comp}_err'])

            # If a G (foreground galaxy) component magnitude is listed, assume all components are and produce a total
            if comp == '_G':
                total = 1

            fs.append(f)
            fes.append(fe)

            f_fields.append(f'flux{comp}' if comp != '' else 'flux_total')
            fe_fields.append(f'flux{comp}_err')

        if total == 1:
            fs.append(sum(fs))
            fes.append(np.linalg.norm(fes))

            if 'flux_total' in f_fields:
                raise ValueError('mags file has both all components and a total mag, this cannot be handled.')

            f_fields.append('flux_total')
            fe_fields.append('flux_err')

        # If this combination of filter and wavelength does not exist yet in lqso.sed
        if lqso.sed[(lqso.sed['filter'] == row['filter']) & (lqso.sed.source == row.source)].empty:
            # Add it
            lqso.sed.loc[len(lqso.sed.index), ['filter', 'wavelength', 'source', 'upper_limit'] + f_fields + fe_fields] =\
                                                        [row['filter'], wavelength, row.source, row.lower_limit] + fs + fes
        else:
            # It exists, overwrite
            # Find the index of the row with same filter and source
            index = lqso.sed.index[(lqso.sed['filter'] == row['filter']) & (lqso.sed.source == row.source)]
            lqso.sed.loc[index, ['filter', 'wavelength', 'source', 'upper_limit'] + f_fields + fe_fields] =\
                                                        [row['filter'], wavelength, row.source, row.lower_limit] + fs + fes

    lqso.save_sed()


def mag_to_flux(telescope, filter, mag, mag_err):
    """
    Converts the given magnitude to flux using the appropriate convertion formula for given telescope filter combination.
    :param telescope:
    :param filter:
    :param mag:
    :param mag_err:
    :return: flux, flux error
    """
    conversion = get_conversion(telescope, filter)

    zp = FILTER_PROPERTIES[(FILTER_PROPERTIES.telescope == telescope) * (FILTER_PROPERTIES.filtername == filter)].zeropoint.values[0]

    return conversion(filter, mag, mag_err, zeropoint=zp)


def get_conversion(telescope, filter):
    """
    Returns the function that correctly converts magnitude to flux as given by data/filter.csv.
    :param telescope:
    :param filter:
    :return:
    """
    conv = FILTER_PROPERTIES[(FILTER_PROPERTIES.telescope == telescope) * (FILTER_PROPERTIES.filtername == filter)].conversion.values[0]
    # eval() evaluates the string for which function it is, then returns that function
    return eval(f'conversion_{conv}')


def conversion_vega(filter, mag, mag_err, zeropoint=3650e3):
    """
    Converts magnitudes in the Vega system to fluxes.
    :param filter:
    :param mag:
    :param mag_err:
    :param zeropoint: zeropoint of magnitude system
    :return:
    """
    flux = zeropoint * np.power(10., - mag / 2.5)
    flux_err = np.abs(- np.log(10) / 2.5 * flux * mag_err)
    return flux, flux_err


def conversion_AB(filter, mag, mag_err, zeropoint=3631e3):
    return conversion_vega(filter, mag, mag_err, zeropoint=zeropoint)


SDSS_b = {
    'u': 1.4e-10,
    'g': 9e-11,
    'r': 1.2e-10,
    'i': 1.8e-10,
    'z': 7.4e-10
}


def conversion_SDSS(filter, mag, mag_err, zeropoint=3631e3):
    flux = 2. * SDSS_b[filter] * np.sinh(- mag * np.log(10.) / 2.5 - np.log(SDSS_b[filter])) * zeropoint
    flux_err = np.abs(2. * SDSS_b[filter] * np.cosh(- mag * np.log(10.) / 2.5 - np.log(SDSS_b[filter])) * - np.log(10.) / 2.5 * mag_err)
    return flux, flux_err
