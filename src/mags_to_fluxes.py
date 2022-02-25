import os.path

import numpy as np
import pandas as pd

from src.lensed_qso import LensedQSO


QSO_PROPERTIES = pd.read_csv(os.path.join('data', 'filter.csv'))


def mags_to_fluxes(galaxy, mags_file='mags.csv', sed_file='sed.csv'):
    """
    Reads the mags.csv file of galaxy and converts the magnitudes to fluxes, using the appropriate conversion formulas.
    Saves these fluxes to sed.csv.
    :param galaxy:
    :param sed_file:
    :param mags_file:
    :return: DataFrame with fluxes
    """
    mags_csv = pd.read_csv(os.path.join('data', galaxy, mags_file))
    mags_csv.fillna(0, inplace=True)

    lqso = LensedQSO(galaxy)

    for i, row in mags_csv.iterrows():
        try:
            conversion = get_conversion(row.telescope, row['filter'])

            # Find the wavelength of the filter
            wavelength = QSO_PROPERTIES[(QSO_PROPERTIES.telescope == row.telescope) * (QSO_PROPERTIES.filtername == row['filter'])].central_wavelength.values[0]
        except IndexError:
            print(f'No conversion found for {row.telescope} {row["filter"]} filter in data/filter.csv.')
            continue

        # Flag that is set to true if there is a G component, then it is assumed that all components are listed and
        # a total flux can be calculated
        total = False

        # Used to store fluxes and flux_errs for all components
        fs = []
        fes = []
        # Calculate fluxes and flux errors for each component
        for comp in ['_G', '_A', '_B', '_C', '_D']:
            # magnitude 0 means no magnitude given
            if row[f'mag{comp}'] == 0:
                fs.append(0.)
                fes.append(0.)
                continue

            # If a G (foreground galaxy) component magnitude is listed, assume all components are and produce a total
            if comp == '_G':
                total = True

            f, fe = conversion(row.filter, row[f'mag{comp}'], row[f'mag{comp}_err'])

            fs.append(f)
            fes.append(fe)

        total_f = sum(fs) * total
        total_fe = sum(fes) * total

        # If this combination of filter and wavelength does not exist yet in lqso.sed
        if lqso.sed[(lqso.sed['filter'] == row['filter']) & (lqso.sed.source == row.source)].empty:
            # Add it
            lqso.sed.loc[len(lqso.sed.index)] = [row['filter'], wavelength, total_f, total_fe,
                                                 fs[0], fes[0], fs[1], fes[1], fs[2], fes[2], fs[3], fes[3], fs[4], fes[4],
                                                 row.source, row.lower_limit]
        else:
            # It exists, overwrite
            # Find the index of the row with same filter and source
            index = lqso.sed.index[(lqso.sed['filter'] == row['filter']) & (lqso.sed.source == row.source)]
            lqso.sed.loc[index] = [row['filter'], wavelength, total_f, total_fe,
                                                 fs[0], fes[0], fs[1], fes[1], fs[2], fes[2], fs[3], fes[3], fs[4], fes[4],
                                                 row.source, row.lower_limit]

    lqso.sed.to_csv(os.path.join('data', galaxy, sed_file), index=False)


def get_conversion(telescope, filter):
    """
    Returns the function that correctly converts magnitude to flux as given by data/filter.csv.
    :param telescope:
    :param filter:
    :return:
    """
    conv = QSO_PROPERTIES[(QSO_PROPERTIES.telescope == telescope) * (QSO_PROPERTIES.filtername == filter)].conversion.values[0]
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


def conversion_AB(filter, mag, mag_err):
    return conversion_vega(filter, mag, mag_err, zeropoint=3631e3)


SDSS_b = {
    'u': 1.4e-10,
    'g': 9e-11,
    'r': 1.2e-10,
    'i': 1.8e-10,
    'z': 7.4e-10
}


def conversion_SDSS(filter, mag, mag_err):
    zeropoint = 3631e3  # mJy
    flux = 2. * SDSS_b[filter] * np.sinh(- mag * np.log(10.) / 2.5 - np.log(SDSS_b[filter])) * zeropoint
    flux_err = np.abs(2. * SDSS_b[filter] * np.cosh(- mag * np.log(10.) / 2.5 - np.log(SDSS_b[filter])) * - np.log(10.) / 2.5 * mag_err)
    return flux, flux_err
