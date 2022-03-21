# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 10:12:08 2022
model subtraction shizzle
@author: Femke
"""

import src.model_sed as mod
from src.model_sed import fit
from src.model_sed import interp_fluxes
import matplotlib.pyplot as plt
import numpy as np
from src.filters import get_filename
import pandas as pd
import os

def model_subtraction(lqso):
    """
    lqso is the class of the galaxy on which we are working
    model is the name of the model that we wish to subtract
    scalar is by which number this model needs to be multiplied
    """
    #Fitting gives the proper model and the scalar
    model_name, scalar, error = fit(lqso, save_plots=False)
    #The wavelengths of the model
    x_model = mod.MODELS[model_name]['wavelength']* (1 + lqso.props['z_lens'].values[0])
    #The fluxes of the model
    y_model = scalar * mod.MODELS[model_name]['flux']
    
    #how to plot it
    #fig, ax = lqso.plot_spectrum(loglog=True) 
    #ax.plot(x_model, y_model, label = model_name)
    
    #an empty list to store the fluxes
    lqso.sed['flux_sub'] = 0.
    lqso.sed['flux_sub_err'] = 0.
    list_sub=[]
    list_sub_err=[]
    
    #iterating over the sed file
    for i, row in lqso.sed.iterrows(): 
        
        #if there is an upper limit, print it. Calculated the same way now
        if  row['upper_limit']:
            print(f'upper limit for row {i}')
        
        #exclude too large lambdas
        if row['wavelength'] >= np.max(x_model):
            list_sub.append( row['flux_total'])
            list_sub_err.append(row['flux_err'])
            continue
            
        #exclude too small lambdas
        elif row['wavelength'] <= np.min(x_model):  
            list_sub.append( row['flux_total'])
            list_sub_err.append(row['flux_err'])
            continue
            
        
            
        #if the entire row is empty
        elif row['flux_G'] == 0. and row['flux_A'] == 0. and row['flux_B'] == 0.\
                                and row['flux_C'] == 0. and row['flux_D'] == 0. and row['flux_total']==0:
            list_sub.append(0)
            list_sub_err.append(0)
            continue
            
        #if only a total flux is known, we are going to be subtracting the model
        elif row['flux_G'] == 0. and row['flux_A'] == 0. and row['flux_B'] == 0.\
                                and row['flux_C'] == 0. and row['flux_D'] == 0. and row['flux_total']!=0:
                                    
            #check if there is a telescope mentioned, otherwise just take closest wavelength value
            if row['telescope']==0:
                print(f'no telescope in sed file for row {i}')
                list_sub.append( float(row['flux_total']) - scalar * np.interp(row['wavelength'], xp=x_model, fp=y_model))
                list_sub_err.append(np.sqrt(row['flux_err'] ** 2 + (scalar * error) ** 2)) 
                continue
                
            #check if the telescope has a name in filters.csv
            #if not, just take the closest value
            filename = get_filename(row['telescope'], row['filter'])
            if pd.isnull(filename):
                print( 'filename not in filters.csv for', row['telescope'])
                list_sub.append( float(row['flux_total'])- scalar * np.interp(row['wavelength'], xp=x_model, fp=y_model))
                list_sub_err.append(np.sqrt(row['flux_err']**2 + (scalar * error) **2 ))
                continue
            
            #read in the filterprofile
            filterprofile = pd.read_csv (os.path.join('data','Filterprofiles','TXT',filename), \
                delim_whitespace=True, header=None,usecols=[ 0,1], names=['wavelength', 'transmission'])
            
            #middel over het filterprofiel: neem het stukje model op je filterprofiel
            x_model_range=x_model[(x_model >= min(filterprofile['wavelength'])) * (x_model <= max(filterprofile['wavelength']))]
            y_model_range=y_model[(x_model >= min(filterprofile['wavelength'])) * (x_model <= max(filterprofile['wavelength']))]
            #Neem de weights = de waarden van je filterprofiel op de range van je model
            weights_filter = np.interp(x_model_range, xp=filterprofile['wavelength'], fp=filterprofile['transmission'])
            #Neem het weighted average = de waarden van je model op de filterrange
            average_model = np.average (y_model_range, weights=weights_filter)
            
            list_sub.append( float(row['flux_total'])- scalar * average_model)
            list_sub_err.append(np.sqrt(row['flux_err']**2 + (error * scalar)**2))
            continue
            
            
        #if the componentwise data is available, use that instead of subtracting the data    
        else:
            list_sub.append( float(row['flux_A']) + float(row['flux_B']) + float(row['flux_C']) + float(row['flux_D']))
            list_sub_err.append(np.sqrt(row['flux_A_err']**2 +row['flux_B_err']**2 + row['flux_C_err']**2 + row['flux_D_err']**2 ))
    
    #write to the sed
    lqso.sed['flux_sub']=list_sub  
    lqso.sed['flux_sub_err']=list_sub_err
    
    
    lqso.plot_spectrum(component='_sub')
    plt.vlines(np.max(x_model), 0.9*np.min(list_sub), np.max(lqso.sed['flux_total']), alpha=0.5)
    plt.vlines(np.min(x_model), 0.9*np.min(list_sub), np.max(lqso.sed['flux_total']), alpha=0.5, label='model boundary')
    plt.legend()
    lqso.save_sed()
