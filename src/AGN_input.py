# -*- coding: utf-8 -*-
"""
To write the AGN input as a print
"""
import pandas as pd
import os
from src.lensed_qso import LensedQSO
from src.filters import get_filename

FILTER_PROPS_PATH = os.path.join('data', 'filter.csv')
FILTER_PROPERTIES = pd.read_csv(FILTER_PROPS_PATH)


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
    print('BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB')


def AGN_input_3(galaxy=None):
    lqso = LensedQSO(galaxy)

    # Iterate over each row to get the dictionary settings
    for i, row in FILTER_PROPERTIES.iterrows():

        tel = FILTER_PROPERTIES.telescope.values[i]
        fil = FILTER_PROPERTIES.filtername.values[i]
        file = FILTER_PROPERTIES.file.values[i]

        # als je geen galaxy geeft zet ze allemaal op nul
        if galaxy is None:
            if pd.isnull(file):
                continue

            print(f'    filters[_{tel}_{fil}]=True')

        # Als je wel een glaxy geeft zet je hem op true wanneer hij er in zit, false als niet
        else:
            # if [tel,fil] in lqso.sed[['telescope', 'filter']].values:
            if lqso.filter_sed(component='_sub').loc[(lqso.filter_sed(component='_sub')['telescope'] == tel) * (
                    lqso.filter_sed(component='_sub')['filter'] == fil)].shape[0] > 0:
                f = fil.replace("'", 'prime').replace('-', 'dash').replace('/', 'slash').replace('.', 'dot')
                print(f"    filters['_{tel.replace(' ', 'space')}_{f}']=True")
            else:
                f = fil.replace("'", 'prime').replace('-', 'dash').replace('/', 'slash').replace('.', 'dot')
                print(f"    filters['_{tel.replace(' ', 'space')}_{f}']=False")
