import os

import pandas as pd

FILTER_PROPERTIES = pd.read_csv(os.path.join('data', 'filter.csv'))


def get_wavelength(telescope, filter):
    return FILTER_PROPERTIES[(FILTER_PROPERTIES.telescope == telescope) * (FILTER_PROPERTIES.filtername == filter)].central_wavelength.values[0]
