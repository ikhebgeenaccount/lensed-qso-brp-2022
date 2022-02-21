# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 11:35:07 2022

@author: Femke
"""

#Import the required packages
import numpy as np
import matplotlib.pyplot as plt

#read in the data
inputFile = "B1152+200.csv"
paper_name = np.loadtxt( inputFile, delimiter=';', usecols=0, skiprows=1 , dtype='str')
filter_name = np.loadtxt( inputFile, delimiter=';', usecols=1, skiprows=1 , dtype='str')
wavelength = np.loadtxt( inputFile, delimiter=';', usecols=2, skiprows=1 )  
flux = np.loadtxt( inputFile, delimiter=';', usecols=3, skiprows=1 ) 
flux_error = np.loadtxt( inputFile, delimiter=';', usecols=4, skiprows=1 )
flux_A = np.loadtxt( inputFile, delimiter=';', usecols=5, skiprows=1)
flux_B = np.loadtxt( inputFile, delimiter=';', usecols=6, skiprows=1 )

#loglog without zero fluxes
plt.figure(0, figsize=(10,8)) 
for i in range(len(flux[flux>0])):
    plt.errorbar(wavelength[flux>0][i], flux[flux>0][i], yerr=flux_error[flux>0][i], fmt='.', label=paper_name[flux>0][i])
    plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.xlabel('log wavelength [Angstrom]')
plt.ylabel('log flux [mJy]')
plt.title('log-log plot, zero fluxes left out')
plt.savefig('loglog_B1152+200.pdf', bboxinches='tight')


#loglog with zero fluxes
plt.figure(1, figsize=(10,8)) 
plt.errorbar(wavelength[flux>=0], flux[flux>=0], yerr=flux_error[flux>=0], fmt='.')
plt.xscale("log")
plt.xlabel('log wavelength [Angstrom]')
plt.ylabel('log flux [mJy]')
plt.title('log-log plot, zero fluxes left in')
plt.savefig('log_B1152+200.pdf', bboxinches='tight')