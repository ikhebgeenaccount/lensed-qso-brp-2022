import warnings

import pandas as pd
import requests


def extract_table(url):
    """
    Extracts all data from tables on an individual quasar's CASTLES page.
    Combines any tables found, notifies if values of magnitudes are outside the range [15,25], as these are probably fluxes.
    :param url: URL to individual quasar's CASTLES page
    :return: pandas DataFrame containinng filter column followed by mag_{component}, mag_{component}_err columns
    """
    res = requests.get(url)

    # Read all tables in the html
    dfs = pd.read_html(res.text)

    # Loop through each table and add data to lists
    data = {
        'filter': [],
    }
    might_be_flux = False
    for df in dfs:
        # Only take rows with fluxes (usually these are mags)
        rel = df[df[df.columns[0]] == 'fluxes']

        for i, c in enumerate(rel.columns):
            if i > 1:  # From 3rd column onward is data
                # Split column data with
                spl = rel[c].str.split('±')

                sed_name = f'mag_{c[1]}'
                sed_err_name = f'mag_{c[1]}_err'

                # Check if components are in data already
                # 0 is default value
                if sed_name not in data:
                    data[sed_name] = [0] * len(data['filter'])
                else:
                    data[sed_name] += [0] * rel.shape[0]

                if sed_err_name not in data:
                    data[sed_err_name] = [0] * len(data['filter'])
                else:
                    data[sed_err_name] += [0] * rel.shape[0]

                # Loop through all entries and add them to
                for j, content in enumerate(spl):
                    # Index of current entry is dependent on what went before, as we append multiple tables
                    ind = len(data['filter']) - rel.shape[0] + j

                    # Sometimes it will not have an error, or nothing at all
                    if len(content) > 1:
                        mag, err = content[0], content[1]
                    else:
                        mag = content[0]
                        err = 0

                    # Convert to floats, if mag is unable to convert, continue to next values as it's a string
                    # If error is unable to convert, set it to 0
                    try:
                        mag = float(mag)
                    except ValueError:
                        continue

                    try:
                        err = float(err)
                    except ValueError:
                        err = 0

                    # Update flag if outside of general range
                    if not 15 < mag < 25:
                        might_be_flux = True

                    data[sed_name][ind] = mag
                    data[sed_err_name][ind] = err

            elif i == 1:  # 2nd column contains filters
                data['filter'] += list(rel[c].values)

    # Make sure all arrays are same length
    length = 0
    for c in data:
        length = len(data[c]) if len(data[c]) > length else length

    for c in data:
        if len(data[c]) != length:
            data[c] += [0] * (length - len(data[c]))

    if might_be_flux:
        warnings.warn(f'{url} might contain fluxes, check it.')

    return pd.DataFrame(data)
