import pandas as pd
import requests


def _extract_table(url):
    res = requests.get(url)

    # Read all tables in the html
    dfs = pd.read_html(res.text)

    # Loop through each table and add data to lists
    data = {
        'filter': [],
    }
    for df in dfs:
        # Only take rows with fluxes (usually these are mags)
        rel = df[df['Observations'] == 'fluxes']  # TODO: this won't work as cols are ('Data from ...', 'Observations')

        # TODO: how to know if column is data column (component column)?
        # TODO: or filter column
        for c in rel.columns:
            pass
            # TODO: split column data with
            df[c].str.split('±')
            # then loop through this array to save magnitudes

        # TODO: Check range of values, if fluxes between [15,25] usually mags, otherwise warn that they might be fluxes


_extract_table('https://lweb.cfa.harvard.edu/castles/Individual/CTQ414.html')
_extract_table('https://lweb.cfa.harvard.edu/castles/Individual/B1608.html')
_extract_table('https://lweb.cfa.harvard.edu/castles/noimages.html')
