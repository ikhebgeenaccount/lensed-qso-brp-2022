# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:03:40 2022

@author: Femke
"""
import os
import pandas as pd

def percent_to_fraction(filtername, newfiltername):
    file = pd.read_csv (os.path.join(App.config().get(section='GENERAL', option='data_dir'),'Filterprofiles','TXT',filtername), delim_whitespace=True, header=None,usecols=[ 0,1], names=['wavelength', 'transmission'])
    file['transmission'] = file['transmission']/100
    file.to_csv(os.path.join(App.config().get(section='GENERAL', option='data_dir'),'Filterprofiles','TXT',newfiltername))
    
