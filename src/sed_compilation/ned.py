import pandas as pd
import numpy as np
import os.path
import warnings

from src.app import App
from src.sed_compilation import ned_api


def update_sed(lqso):
    """
    Adds to the SED of lqso all data points found in NED.
    :param lqso:
    :return:
    """
    # Try all aliases until we find something
    for alias in lqso.aliases:
        ned_sed = ned_api.access_sed(alias)

        if ned_sed.shape[0] > 0:
            break

    # TODO: qualifiers
    # TODO: unit conversions (using astropy?)
    return ned_sed


def read_table_file(lqso, ned_file='ned.txt', wavelength_conversion=1e4, flux_conversion=1e3, qualifier=None, allowed_sources=None):
    """
    Reads a ned.txt file of a NED bar-separated table and enters it into the SED of the galaxy.
    :param galaxy:
    :param ned_file:
    :param wavelength_conversion: default: um to Angstrom
    :param flux_conversion: default: Jy to mJy
    :param qualifier:
    :param allowed_sources: List of sources to use, e.g. ['Chandra']
    :return:
    """
    ned_df = pd.read_csv(os.path.join(App.config().get(section='GENERAL', option='data_dir'), lqso.name, ned_file), delimiter='|')

    if allowed_sources is None:
        raise ValueError('No allowed sources, enter a list of allowed sources using ned_table_to_sed(..., allowed_sources=[...])')

    # Filter out non-measurement
    fs = np.sum(ned_df['Flux Density'].isna())
    if fs > 0:
        warnings.warn('Filtered out ' + str(fs) + ' measurements, check if this is correct!')
        ned_df = ned_df[ned_df['Flux Density'].notna()]

    # Give fields we want proper name, convertFalse
    ned_df['wavelength'] = ned_df['Wavelength'] * wavelength_conversion
    ned_df['flux_total'] = ned_df['Flux Density'] * flux_conversion
    ned_df['flux_err'] = ned_df['Upper limit of uncertainty'] * flux_conversion

    ned_df['filter'] = ned_df['Observed Passband'].apply(lambda v: v.strip())
    ned_df['Qualifiers'] = ned_df['Qualifiers'].apply(lambda v: v.strip())

    qualis = ned_df['Qualifiers'].unique()

    # Check if qualifier exists and exists in qualis
    qi = None
    if qualifier is None:
        qi = input(f'Found qualifiers: {qualis}, input index of which one to use (zero-indexed):')
    elif qualifier not in qualis:
        qi = input(f'Given qualifier {qualifier} not found in qualifiers {qualis}, input index of which one to use (zero-indexed):')

    if qi is not None:
        qualifier = qualis[int(qi)]

    if qualifier not in qualis:
        raise ValueError('Qualifier ' + qualifier + ' not found in qualifiers of given table.')

    # Filter on qualifier
    q_sel = ned_df.loc[(ned_df['Qualifiers'] == qualifier)]
    wls = q_sel['wavelength'].unique()

    # Combine measurements for each unique wavelength into single value and error
    for wl in wls:
        sel = q_sel.loc[(q_sel['wavelength'] == wl)]

        # Check for proper size of selection, if not exactly one warn the user and continue to next wavelength
        if sel.shape[0] > 1:
            warnings.warn('Multiple values found for wavelength ' + str(wl) +', qualifier ' + qualifier)
            continue

        if sel.shape[0] == 0:
            warnings.warn('No values found for wavelength ' + str(wl) +', qualifier ' + qualifier)
            continue

        # Get the required data
        flux_total = sel['flux_total'].values[0]
        flux_err = sel['flux_err'].values[0]
        filter = sel['filter'].values[0]

        source = ''

        # Check if source is in sources to skip
        # Default: skip
        skip = True
        # Check for every skip_sources if it occurs in the filter, if so, set skip to True
        if allowed_sources is not None:
            for als in allowed_sources:
                if als.lower() in filter.lower():
                    source = als  # Set source of this entry to the allowed source name
                    skip = False
                    break

        if skip:
            continue

        if not np.all(np.where(sel['Upper limit of uncertainty'] == sel['Lower limit of uncertainty'], 1, 0)):
            # print(ned_df[np.where(ned_df['Upper limit of uncertainty'] == ned_df['Lower limit of uncertainty'], 1, 0)][['Upper limit of uncertainty', 'Lower limit of uncertainty']])
            raise NotImplementedError('Upper and lower limit of uncertainty are different, this cannot be handled.')

        if np.all(np.where(sel['Upper limit of Flux Density'] == np.nan, 1, 0)) or np.all(np.where(sel['Lower limit of Flux Density'] == np.nan, 1, 0)):
            raise NotImplementedError('Table has upper or lower limit for flux density, this cannot be handled.')

        # Check if entry is not in lqso.sed yet
        if lqso.sed[(lqso.sed['wavelength'] == wl) & (lqso.sed['filter'] == filter)].empty:
            # Add it
            lqso.sed = lqso.sed.append({
                'wavelength': wl,
                'flux_total': flux_total,
                'flux_err': flux_err,
                'filter': filter,
                'source': source
            }, ignore_index=True)
        else:
            # Overwrite it
            index = lqso.sed.index[(lqso.sed['wavelength'] == wl) & (lqso.sed['filter'] == filter)]
            lqso.sed.loc[index, ['wavelength', 'flux_total', 'flux_err', 'filter', 'source']] = \
                [wl, flux_total, flux_err, filter, source]

    # Save to SED
    lqso.save_sed()
