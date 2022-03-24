# -*- coding: utf-8 -*-
"""
To write the AGN input as a print
"""
import pandas as pd
import os
from src.lensed_qso import LensedQSO
from src.filters import get_filename

def AGN_input_1():
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

def AGN_input_2():    
    for i, row in FILTER_PROPERTIES.iterrows():
        
        tel=FILTER_PROPERTIES.telescope.values[i]
        TEL=tel.upper()
        fil=FILTER_PROPERTIES.filtername.values[i]
        file=FILTER_PROPERTIES.file.values[i]
        
        if pd.isnull(file):
            continue
        
        print(f"    if filters[{tel}_{fil}]:")
        print(f"        files.append({tel}_{fil}_file)")
        print(f"        lambdas.append({tel}_{fil}_lambda)")
        print(f"        factors.append({tel}_{fil}_factor)")
        print('')
    print('BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB')
    
def AGN_input_3(galaxy=None):
    lqso = LensedQSO(galaxy)
    print(lqso.sed[['telescope', 'filter']].values)
    for i, row in FILTER_PROPERTIES.iterrows():
        
        #Iterate over each row to get the dictionary settings
        
        tel=FILTER_PROPERTIES.telescope.values[i]
        TEL=tel.upper()
        fil=FILTER_PROPERTIES.filtername.values[i]
        file=FILTER_PROPERTIES.file.values[i]
        
        
        #als je geen galaxy geeft zet ze allemaal op nul
        if pd.isnull(galaxy):     
            if pd.isnull(file):
                continue
            
            print(f'    filters[{tel}_{fil}]=True')
            
            
            
         #Als je wel een glaxy geeft zet je hem op true wanneer hij er in zit, false als niet   
        else:
            # if [tel,fil] in lqso.sed[['telescope', 'filter']].values:
            if lqso.filter_sed().loc[(lqso.filter_sed()['telescope'] == tel) * (lqso.filter_sed()['filter'] == fil)].shape[0] > 0:
                print(f"    filters['{tel}_{fil}']=True")
            else:
                print(f"    filters['{tel}_{fil}]=False'") 
            
            
            
            
            
            
            
            