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
    model_name, scalar = fit(lqso, save_plots=False)
    #The wavelengths of the model
    x_model = mod.MODELS[model_name]['wavelength']* (1 + lqso.props['z_lens'].values[0])
    #The fluxes of the model
    y_model = scalar * mod.MODELS[model_name]['flux']
    
    #how to plot it
    #fig, ax = lqso.plot_spectrum(loglog=True) 
    #ax.plot(x_model, y_model, label = model_name)
    
    #make a new column in the sed file 
    list_sub=[]
    #now we want to access componentwise the sed file, later write it to sed_sub file
    for i, row in lqso.sed.iterrows():
        if row['wavelength'] >= np.max(x_model):
            list_sub.append( 'nope')
        elif row['flux_G'] == 0. and row['flux_A'] == 0.:
            list_sub.append( row['flux_total']- scalar * np.interp(row['wavelength'], xp=x_model, fp=y_model))
        else:
            list_sub.append( row['flux_A'] + row['flux_B'] + row['flux_C'] + row['flux_D'])
    lqso.sed['flux_sub']=list_sub    
    
    lqso.save_sed()
    
