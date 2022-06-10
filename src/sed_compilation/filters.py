import glob
import os
import re

import pandas as pd

from src.app import App

FILTER_PROPERTIES = pd.read_csv(App.config().get(section='GENERAL', option='filters_properties_file'))


def get_wavelength(telescope, tfilter):
    return FILTER_PROPERTIES[(FILTER_PROPERTIES.telescope == telescope) * (FILTER_PROPERTIES.filtername == tfilter)].central_wavelength.values[0]


def get_filename(telescope, tfilter):
    try:
        return FILTER_PROPERTIES[(FILTER_PROPERTIES.telescope == telescope) * (FILTER_PROPERTIES.filtername == tfilter)].file.values[0]
    except IndexError:
        print(telescope, tfilter)


def populate_filter_profile_path_column(profiles_dir=None):
    if profiles_dir is None:
        profiles_dir = os.path.join(App.config().get(section='GENERAL', option='data_dir'), 'Filterprofiles', 'TXT')

    filters = found = 0

    for i, row in FILTER_PROPERTIES.iterrows():
        filters += 1
        fp_pattern = f'{row["telescope"]}.{row["filtername"]}.[a-z]*'.replace('/', '.').replace("'", 'prime')

        matches = [f for f in os.listdir(profiles_dir) if re.search(fp_pattern, f)]

        if len(matches) == 1:
            found += 1
            FILTER_PROPERTIES.loc[i, 'file'] = matches[0]

    FILTER_PROPERTIES.to_csv(App.config().get(section='GENERAL', option='filters_properties_file'), index=False)
