# -*- coding: utf-8 -*-
"""
To write the AGN input as a print
"""
import pandas as pd
import os
def AGN_input():
    FILTER_PROPS_PATH = os.path.join('data', 'filter.csv')
    FILTER_PROPERTIES = pd.read_csv(FILTER_PROPS_PATH)
    
    for i, row in FILTER_PROPERTIES.iterrows():
        
        tel=FILTER_PROPERTIES.telescope.values[i]
        TEL=tel.upper()
        fil=FILTER_PROPERTIES.filtername.values[i]
        file=FILTER_PROPERTIES.file.values[i]
        
        if pd.isnull(file):
            continue
        #print(tel, fil, file)
        print(f"    {tel}_{fil}_file = path + 'models/FILTERS/{TEL}/{file}'")
        print(f"    {tel}_{fil}_lambda,{tel}_{fil}_factor = np.loadtxt({tel}_{fil}_file, usecols=(0,1), unpack = True)")
        print("")
    print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAHHHHHHHHHHHH')
    
    for i, row in FILTER_PROPERTIES.iterrows():
        
        tel=FILTER_PROPERTIES.telescope.values[i]
        TEL=tel.upper()
        fil=FILTER_PROPERTIES.filtername.values[i]
        file=FILTER_PROPERTIES.file.values[i]
        
        print(f'    if filters[{tel}_{fil}]:')
        print(f'        files.append({tel}_{fil}_file)')
        print(f'        lambdas.append({tel}_{fil}_lambda)')
        print(f'        factors.append({tel}_{fil}_factor)')
        print('')
    print('BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB')
    
    for i, row in FILTER_PROPERTIES.iterrows():
        
        tel=FILTER_PROPERTIES.telescope.values[i]
        TEL=tel.upper()
        fil=FILTER_PROPERTIES.filtername.values[i]
        file=FILTER_PROPERTIES.file.values[i]
        print(f'    filters[{tel}_{fil}]=True')