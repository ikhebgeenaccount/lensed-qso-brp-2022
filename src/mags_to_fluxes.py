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
        except (IndexError, NameError):
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

            try:
                f, fe = mag_to_flux(row.telescope, row['filter'], row[f'mag{comp}'], row[f'mag{comp}_err'])
            except NameError:
                print(f'No conversion found for {row.telescope} {row["filter"]} filter in data/filter.csv.')
                continue

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


def mag_ratio_split_total_flux(lqso, ratio_source, max_filter_dist=1e3, components=None, overwrite=False):
    """
    Calculates the flux ratios of the components based on the components magnitudes.
    Applies this ratio to the closest filter (based on wavelength) it can find, if that wavelength is not split into components already.
    :param lqso: LensedQSO object
    :param ratio_source: indicates which source to use to find the ratio.
    :param max_filter_dist: maximum distance between filter that ratios are taken from and filter of which total will be split
    """
    if components is None:
        components = ['_G', '_A', '_B', '_C', '_D']
    source_is = lqso.mags.loc[lqso.mags.source == ratio_source].index

    if source_is.shape[0] > lqso.filter_sed().shape[0]:
        raise ValueError('Found more mag filters than filtered SED rows, which cannot be handled.')
        # See note below, at while loop

    n_comps = 0
    for c in components:
        n_comps += lqso.mags.loc[source_is[0], f'mag{c}'] > 0.

    print(lqso.name, f'has {n_comps} components.')

    # subbeds saves the array subbed for each index in source_is
    # Is used to compare the different closest filters and find wavelengths that have the same closest filter
    subbeds = []

    # Only find closest wavelength for allowed sources and flux_total > 0
    # So use lqso.filter_sed() to find wls
    fsed = lqso.filter_sed()

    for i in source_is:
        # Select mag based on source, calculate ratios
        r = lqso.mags.loc[i]

        # Find the wavelength of the used filter
        wlf = get_wavelength(r['telescope'], r['filter'])

        subbed = fsed['wavelength'] - wlf

        subbeds.append(np.abs(subbed))

        closest_wl = fsed['wavelength'].iloc[np.argmin(np.abs(subbed))]
        dist = subbed[np.argmin(np.abs(subbed))]

        print(f'For {r["telescope"]}, {r["filter"]} filter:\n\tclosest in SED is {fsed["filter"].iloc[np.argmin(np.abs(subbed))]} at {closest_wl} with dist {dist}\n\tSource: {fsed["source"].iloc[np.argmin(np.abs(subbed))]}')

        if dist > max_filter_dist:
            print('\tDistance between these filters is too large: skipping')

    # Find the closest index for each subbed in subbeds
    closest_is = np.argmin(subbeds, axis=1)

    # Array used to translate from closest_is to sed index
    fsed_lids = fsed.index

    print('SED rows for each mag row before duplicate check:')
    print(closest_is)
    ks = np.ones(len(closest_is), dtype=int)

    # While there are any duplicate values
    # Note: this breaks if there are more mag ratio source rows than sed rows
    while len(np.unique(closest_is)) != len(closest_is):
        us, ins, counts = np.unique(closest_is, return_index=True, return_counts=True)

        # Find duplicate values in closest_is
        # Take closest_is at indices for which the counts > 1
        dups = np.array(closest_is)[ins[counts > 1]]

        for i in dups:
            # Find which of the duplicates is closest to the target wavelength
            fis = np.where(closest_is == i)[0]  # [0] Since np.where returns tuple for dimensions, but only 1D here
            sm = np.argmin(np.array(subbeds)[fis][:,i])

            # Turn filter indices fis into a list so we can remove the lowest value as that one will remain the same
            fis = list(fis)
            fis.remove(sm)

            # Now find the kth best for the remaining cases
            for cid in fis:
                closest_is[cid] = np.where(subbeds[cid] == np.partition(subbeds[cid], ks[cid])[ks[cid]])[0][0]
                # Add one to k for next attempt at finding best
                ks[cid] += 1

    print('SED rows for each mag row after duplicate check:')
    print(closest_is)
    print(f'Finding these indices took {ks - 1} attempts.')

    for ci, ri in enumerate(source_is):
        # ri is magnitude row
        # ci is corresponding index in closest_is, so row in SED that we need to fill

        # Check dist ok
        if subbeds[ci][closest_is[ci]] > max_filter_dist:
            print(f'For mag row {ri}, closest row in SED {fsed_lids[closest_is[ci]]} is furhter than allowed max_filter_dist {max_filter_dist}, skipping.')
            continue

        # Check no _G value yet
        if fsed['flux_G'].iloc[closest_is[ci]] > 0.:
            print(f'For mag row {ri}, closest row in SED {fsed_lids[closest_is[ci]]} already has foreground value')
            if not overwrite:
                print('\tNot overwriting, if you want to, set overwriting=True')
                continue
            print('\tOverwriting.')

        # Flux ratios are given by m_1 - m_2 = -2.5 log10(f_1/f_2)
        for c in range(n_comps):
            comp = components[c]

            ratios = [1.]

            # Calculate flux ratios
            for cs in range(n_comps):
                # Flux ratio with itself is 1, already in starting value of div
                if c == cs:
                    continue
                print(f'flux ratio {components[cs]}/{comp}={np.power(10., (lqso.mags[f"mag{components[cs]}"].loc[ri] - lqso.mags[f"mag{components[c]}"].loc[ri]) / -2.5)}')

                ratios.append(np.power(10., (lqso.mags[f'mag{components[cs]}'].loc[ri] - lqso.mags[f'mag{components[c]}'].loc[ri]) / -2.5))

            sum_ratios = np.sum(ratios)

            # Calculate error
            stdsq = 0.

            # df_c / dF_T, c is current component comp
            fcft = 1. / sum_ratios
            stdsq += np.power(fcft * lqso.sed.loc[fsed_lids[closest_is[ci]], 'flux_err'], 2.)

            # df_c / dm_c
            fcmc = np.log(10.) / 2.5 * lqso.sed.loc[fsed_lids[closest_is[ci]], 'flux_total'] * (1. - sum_ratios) / np.power(sum_ratios, 2.)
            stdsq += np.power(fcmc * lqso.mags[f'mag{comp}_err'].loc[ri], 2.)

            for cs in range(n_comps):
                if c == cs:
                    continue
                fbmi = np.log(10.) / 2.5 * lqso.sed.loc[fsed_lids[closest_is[ci]], 'flux_total'] * ratios[cs] / np.power(sum_ratios, 2.)
                stdsq += np.power(fbmi * lqso.mags[f'mag{components[c]}_err'].loc[ri], 2.)

            print(f'flux{comp} =', fsed['flux_total'].iloc[closest_is[ci]] / sum_ratios, '+/-', np.sqrt(stdsq))
            lqso.sed.loc[fsed_lids[closest_is[ci]], [f'flux{comp}', f'flux{comp}_err']] = [fsed['flux_total'].iloc[closest_is[ci]] / sum_ratios, np.sqrt(stdsq)]

    lqso.plot_spectrum(component='_G')
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
