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
    morph = 'all' if pd.isnull(lqso.props.lens_type.values[0]) else lqso.props.lens_type.values[0]

    #Fitting gives the model> see model_sed for details
    wavelength, flux, flux_error = fit(lqso, morph = morph)

    #The wavelengths of the model are returned already in restframe by model_sed

    #an empty list to store the fluxes, and new column in the sed
    lqso.sed['flux_sub'] = 0.
    lqso.sed['flux_sub_err'] = 0.
    list_sub=[]
    list_sub_err=[]

    #iterating over the sed file
    for i, row in lqso.filter_sed(component=None).iterrows():
        #if there is an upper limit, print it. Calculated the same way now
        if  row['upper_limit']:
            print(f'upper limit for row {i}')

        #exclude too large lambdas
        if row['wavelength'] >= np.max(wavelength):
            list_sub.append( row['flux_total'])
            list_sub_err.append(row['flux_err'])
            continue

        #exclude too small lambdas
        elif row['wavelength'] <= np.min(wavelength):
            list_sub.append( row['flux_total'])
            list_sub_err.append(row['flux_err'])
            continue

        #if the entire row is empty, subtracted is 0
        elif row['flux_G'] == 0. and row['flux_A'] == 0. and row['flux_B'] == 0.\
                                and row['flux_C'] == 0. and row['flux_D'] == 0. and row['flux_total']==0:
            list_sub.append(0)
            list_sub_err.append(0)
            continue

        #if only a total flux is known, we are going to be subtracting the model
        elif row['flux_G'] == 0. and row['flux_A'] == 0. and row['flux_B'] == 0.\
                                and row['flux_C'] == 0. and row['flux_D'] == 0. and row['flux_total']!=0:

            #check if there is a telescope mentioned, otherwise just take closest wavelength value
            if pd.isnull(row['telescope']):
                print(f'no telescope in sed file for row {i}')
                list_sub.append( float(row['flux_total']) - np.interp(row['wavelength'], xp=wavelength, fp=flux))
                list_sub_err.append(np.sqrt(row['flux_err'] ** 2 + np.interp(row['wavelength'], xp=wavelength, fp=flux_error) ** 2))
                continue

            #check if the telescope has a name in filters.csv
            #if not, just take the closest value
            filename = get_filename(row['telescope'], row['filter'])
            if pd.isnull(filename):
                print( 'filename not in filters.csv for', row['telescope'])
                list_sub.append( float(row['flux_total'])- np.interp(row['wavelength'], xp=wavelength, fp=flux))
                list_sub_err.append(np.sqrt(row['flux_err']**2 + np.interp(row['wavelength'], xp=wavelength, fp=flux_error) **2 ))
                continue

            #read in the filterprofile
            filterprofile = pd.read_csv (os.path.join('data','Filterprofiles','TXT',filename), \
                delim_whitespace=True, header=None,usecols=[ 0,1], names=['wavelength', 'transmission'])

            #middel over het filterprofiel: selecteer eerst het stukje model op je filterprofiel
            x_model_range=wavelength[(wavelength >= min(filterprofile['wavelength'])) & (wavelength <= max(filterprofile['wavelength']))]
            y_model_range=flux[(wavelength >= min(filterprofile['wavelength'])) & (wavelength <= max(filterprofile['wavelength']))]
            error_range = flux_error[(wavelength >= min(filterprofile['wavelength'])) & (wavelength <= max(filterprofile['wavelength']))]

            #Neem de weights = de waarden van je filterprofiel op de range van je model
            weights_filter = np.interp(x_model_range, xp=filterprofile['wavelength'], fp=filterprofile['transmission'])

            #Neem het weighted average = de waarden van je model op de filterrange
            average_model = np.average (y_model_range, weights=weights_filter)
            #average_error = np.average (error_range, weights=weights_filter)

            average_error = (1 / len(weights_filter)*np.sum(weights_filter)) * \
                np.sqrt (np.sum (np.square(weights_filter * error_range)))


            list_sub.append( float(row['flux_total']) -  average_model)
            list_sub_err.append(np.sqrt(row['flux_err']**2 + (average_error)**2))

        #if the componentwise data is available, use that instead of subtracting the data
        else:
            list_sub.append( float(row['flux_A']) + float(row['flux_B']) + float(row['flux_C']) + float(row['flux_D']))
            list_sub_err.append(np.sqrt(row['flux_A_err']**2 +row['flux_B_err']**2 + row['flux_C_err']**2 + row['flux_D_err']**2 ))

    #write to the sed
    lqso.sed.loc[lqso.filter_sed(component=None).index, 'flux_sub']=list_sub
    lqso.sed.loc[lqso.filter_sed(component=None).index, 'flux_sub_err']=list_sub_err

    # Remove magnification
    # Divide by magnification
    lqso.sed['flux_sub_demag'] = lqso.sed['flux_sub'] / lqso.props['magnification'].values[0]
    # Calculate new error
    lqso.sed['flux_sub_demag_err'] = np.sqrt(
        np.power(1. / lqso.props['magnification'].values[0] * lqso.sed['flux_sub_err'], 2.) +
        np.power(lqso.sed['flux_sub'] / np.power(lqso.props['magnification'].values[0], 2.) * lqso.props['magn_err'].values[0], 2.)
    )

    # Change upper limits with errors to upper limit = upper limit + error
    lqso.sed['flux_sub_demag'].loc[lqso.sed['upper_limit'] == 1] += lqso.sed['flux_sub_demag_err'].loc[lqso.sed['upper_limit'] == 1]
    lqso.sed['flux_sub_demag_err'].loc[lqso.sed['upper_limit'] == 1] = 0

    fig, ax = lqso.plot_spectrum(component='_sub')
    fig.savefig(os.path.join('plots', lqso.name, f'{lqso.name}_sub.jpg'))
    fig.savefig(os.path.join('plots', lqso.name, f'{lqso.name}_sub.pdf'))
    plt.vlines(np.max(wavelength), 0.9*np.min(list_sub), np.max(lqso.sed['flux_total']), alpha=0.5)
    plt.vlines(np.min(wavelength), 0.9*np.min(list_sub), np.max(lqso.sed['flux_total']), alpha=0.5, label='model boundary')
    lqso.save_sed()
