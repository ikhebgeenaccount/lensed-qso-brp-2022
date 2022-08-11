import re

import pandas as pd
import numpy as np
import os.path
import warnings

from src.app import App
from src.lensed_qso import LensedQSO
from src.sed_compilation import ned_api


def update_sed(lqso: LensedQSO, ignore_radio=False, ignore_xray=False):
    """
    Adds to the SED of lqso all data points found in NED.
    :param ignore_xray:
    :param ignore_radio:
    :param lqso:
    :return:
    """
    if len(lqso.aliases) == 0:
        raise ValueError('lqso does not have a name.')

    # Try all aliases until we find something
    for alias in lqso.aliases:
        ned_sed = ned_api.access_sed(alias)

        if ned_sed.shape[0] > 0:
            break

    filtered_ned_sed = _filter_df_with_qualifiers(ned_sed, ignore_radio, ignore_xray)

    lqso.update_sed(filtered_ned_sed[['filter', 'wavelength', 'flux_total', 'flux_err', 'bibcode']])
    return filtered_ned_sed[['filter', 'wavelength', 'flux_total', 'flux_err', 'bibcode']]


def read_table_file(lqso: LensedQSO, ned_file='ned.txt', wavelength_conversion=1e4, flux_conversion=1e3, qualifier=None, allowed_sources=None):
    """
    Reads a ned.txt file of a NED bar-separated table and enters it into the SED of the galaxy.
    :param lqso:
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


QUALIFIERS = {}


def _filter_df_with_qualifiers(df, ignore_radio, ignore_xray):
    new_df = pd.DataFrame()

    df['added'] = ['-'] * df.shape[0]

    # For every unique filter:
    #   If already in QUALIFIERS:
    #       Use qualifier from QUALIFIERS if qualifer exists in current df
    #   else:
    for f in np.unique(df['filter'].values):
        # Check if wavelength falls within set bounds
        if ignore_radio and df[df['filter'] == f]['wavelength'].values[0] > App.config().getfloat(section='GENERAL', option='radio_cutoff'):
            _update_column_selection(df, 'added', 'n - radio', selection=df['filter'] == f)
            continue
        if ignore_xray and df[df['filter'] == f]['wavelength'].values[0] < App.config().getfloat(section='GENERAL', option='xray_cutoff'):
            _update_column_selection(df, 'added', 'n - xray', selection=df['filter'] == f)
            continue

        # Check if we need to use a qualifier
        if sum(df['filter'] == f) > 1:
            # Print the complete dataframe
            # The 'added' columns gets updated each iteration
            os.system('cls||clear')  # Clear screen
            print(df)

            # Check if we found this filter previously and stored qualifier applies
            if f in QUALIFIERS and QUALIFIERS[f] in df['qualifiers'][df['filter'] == f].values:
                new_df = new_df.append(df[(df['filter'] == f) & (df['qualifiers'] == QUALIFIERS[f])])

                _update_column_selection(df, 'added', 'y - set quali',
                                         selection=(df['filter'] == f) & (df['qualifiers'] == QUALIFIERS[f]))
                _update_column_selection(df, 'added', 'n - set quali',
                                         selection=(df['filter'] == f) & (df['qualifiers'] != QUALIFIERS[f]))
            else:
                print(df[df['filter'] == f])
                # We need to ask user which entry to use
                print(f'Options:'
                      f'\n\t\'skip\' to skip this filter'
                      f'\n\t\'sum\' to sum all entries'
                      f'\n\t\'sum {df[df["filter"] == f].index[0]} {df[df["filter"] == f].index[1]}\' to sum a subset of entries'
                      f'\n\t\'sel {df[df["filter"] == f].index[0]}\' to select an entry'
                      f'\n\t\'set {df[df["filter"] == f].index[0]}\' to select an entry and set its qualifier as preferred qualifier for this filter')

                valid_response = False
                while not valid_response:
                    response = input(f'Enter index or option for filter {f}: ')

                    re_match = re.match('(sel|set|sum|skip) *([0-9 ]*)', response)

                    valid_response = re_match is not None

                res_type = re_match.group(1)
                entry_ids = re_match.group(2).split(' ')

                if res_type == 'skip':
                    # Add skip as reason
                    _update_column_selection(df, 'added', 'n - skip', selection=df['filter'] == f)
                elif res_type == 'sum':
                    # Sum all entries or listed entries
                    raise NotImplementedError('sum not implemented yet')
                elif res_type == 'sel':
                    # Use values of listed entry
                    new_df = new_df.append(df.iloc[int(entry_ids[0])])
                elif res_type == 'set':
                    # Use values of listed entry and update QUALIFIERS
                    new_df = new_df.append(df.iloc[int(entry_ids[0])])
                    QUALIFIERS[f] = df['qualifiers'].iloc[int(entry_ids[0])]

                print(new_df)
                input()
        else:
            # Only one qualifier, add data
            new_df = new_df.append(df[df['filter'] == f])

    return new_df


def _update_column_selection(df, column, new_value, selection=None, ids=None):
    if selection is not None:
        ids = df[selection].index
    elif ids is not None:
        pass
    else:
        raise ValueError('One of selection or ids needs to be set.')

    for id in ids:
        df.at[id, column] = new_value
    return df
