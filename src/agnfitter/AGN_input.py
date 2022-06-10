# -*- coding: utf-8 -*-
"""
To write the AGN input as a print
"""
import pandas as pd
import os

from src.app import App

FILTER_PROPS_PATH = os.path.join(App.config().get(section='GENERAL', option='data_dir'), 'filter.csv')
FILTER_PROPERTIES = pd.read_csv(FILTER_PROPS_PATH)


def format_filter_name(filter_name):
    return filter_name.replace("'", 'prime').replace('-', 'dash').replace('/', 'slash').replace('.', 'dot')


def format_telescope_name(tel_name):
    return tel_name.replace(' ', 'space')


def get_agnf_filter_path(tel, fil):
    TEL = tel.upper()
    file = FILTER_PROPERTIES.loc[(FILTER_PROPERTIES['telescope'] == tel) & (FILTER_PROPERTIES['filtername'] == fil)].file.values[0]
    return f'models/FILTERS/{TEL}/{file}'


def AGN_input_1():
    for i, row in FILTER_PROPERTIES.iterrows():

        tel = FILTER_PROPERTIES.telescope.values[i]
        TEL = tel.upper()
        tel = tel.replace(' ', 'space')
        fil = FILTER_PROPERTIES.filtername.values[i].replace("'", 'prime').replace('-', 'dash').replace('/', 'slash').replace('.', 'dot')
        file = FILTER_PROPERTIES.file.values[i]

        if pd.isnull(file):
            continue
        # print(tel, fil, file)
        print(f"    _{tel}_{fil}_file = path + 'models/FILTERS/{TEL}/{file}'")
        print(
            f"    _{tel}_{fil}_lambda,_{tel}_{fil}_factor = np.loadtxt(_{tel}_{fil}_file, usecols=(0,1), unpack = True)")
        print("")
    print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAHHHHHHHHHHHH')


def AGN_input_2():
    for i, row in FILTER_PROPERTIES.iterrows():

        tel = FILTER_PROPERTIES.telescope.values[i]
        TEL = tel.upper()
        tel = tel.replace(' ', 'space')
        fil = FILTER_PROPERTIES.filtername.values[i].replace("'", 'prime').replace('-', 'dash').replace('/', 'slash').replace('.', 'dot')
        file = FILTER_PROPERTIES.file.values[i]

        if pd.isnull(file):
            continue

        print(f"        if filters['_{tel}_{fil}']:")
        print(f"            files.append(_{tel}_{fil}_file)")
        print(f"            lambdas.append(_{tel}_{fil}_lambda)")
        print(f"            factors.append(_{tel}_{fil}_factor)")
        print('')

