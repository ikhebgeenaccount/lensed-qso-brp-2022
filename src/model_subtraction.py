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
    lqso.sed['flux_sub']=0.
    lqso.sed['flux_sub_err'] = 0.
    list_sub=[]
    list_sub_err=[]
    
    #select on whether there is component-wise data or only total flux
    for i, row in lqso.sed.iterrows(): #exclude radio and such
        if row['wavelength'] >= np.max(x_model):
            list_sub.append( row['flux_total'])
            list_sub_err.append(row['flux_err'])
        elif row['wavelength'] <= np.min(x_model): #exclude too small lambdas 
            list_sub.append( row['flux_total'])
            list_sub_err.append(row['flux_err'])
        #if only a total flux is known, subtract the model
        elif row['flux_G'] == 0. and row['flux_A'] == 0. and row['flux_B'] == 0.\
                                and row['flux_C'] == 0. and row['flux_D'] == 0.:
            list_sub.append( row['flux_total']- scalar * np.interp(row['wavelength'], xp=x_model, fp=y_model))
            list_sub_err.append(np.sqrt(row['flux_err']**2+error**2))
        #if the componentwise data is available, use that instead of subtracting the data    
        else:
            list_sub.append( row['flux_A'] + row['flux_B'] + row['flux_C'] + row['flux_D'])
            list_sub_err.append(np.sqrt(row['flux_A_err']**2 +row['flux_B_err']**2 + row['flux_C_err']**2 + row['flux_D_err']**2 ))
    #write to the sed
    lqso.sed['flux_sub']=list_sub  
    lqso.sed['flux_sub_err']=list_sub_err
    
    
    lqso.plot_spectrum(component='_sub')
    plt.vlines(np.max(x_model), 0.9*np.min(list_sub), np.max(lqso.sed['flux_total']), alpha=0.5)
    plt.vlines(np.min(x_model), 0.9*np.min(list_sub), np.max(lqso.sed['flux_total']), alpha=0.5, label='model boundary')
    plt.legend()
    lqso.save_sed()
