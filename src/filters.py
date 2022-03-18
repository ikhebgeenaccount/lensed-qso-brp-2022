import glob
import os
import re

import pandas as pd

FILTER_PROPS_PATH = os.path.join('data', 'filter.csv')
FILTER_PROPERTIES = pd.read_csv(FILTER_PROPS_PATH)


def get_wavelength(telescope, tfilter):
    return FILTER_PROPERTIES[(FILTER_PROPERTIES.telescope == telescope) * (FILTER_PROPERTIES.filtername == tfilter)].central_wavelength.values[0]


def populate_filter_profile_path_column(profiles_dir=None):
    if profiles_dir is None:
        profiles_dir = os.path.join('data', 'Filterprofiles', 'TXT')

    filters = found = 0

    for i, row in FILTER_PROPERTIES.iterrows():
        filters += 1
        fp_pattern = f'{row["telescope"]}.{row["filtername"]}.[a-z]*'.replace('/', '.').replace("'", 'prime')

        matches = [f for f in os.listdir(profiles_dir) if re.search(fp_pattern, f)]

        if len(matches) == 1:
            found += 1
            FILTER_PROPERTIES.loc[i, 'file'] = matches[0]

    print(found, filters)

    print(FILTER_PROPERTIES)

    FILTER_PROPERTIES.to_csv(FILTER_PROPS_PATH, index=False)