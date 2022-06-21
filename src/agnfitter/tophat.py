# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 10:35:23 2022

@author: Femke
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def tophat(central,bandwidth,filtername,freq_Ghz=False,energy_Kev=False,length=3000):
    """
    This makes a tophat-like profile function for the radio wavelengths, 
    with transmission =1 within the bandwidth and 0 outside it.
    """
        
    #wavelengths in angstrom or frequencies in hertz
    x_array = np.linspace (float(central) - float(bandwidth), float(central) + float(bandwidth), length)
    
    #corresponding transmission
    transmission_array = np.zeros (length)
    transmission_array[(x_array >= (float(central) - (float(bandwidth) /2)))\
                * (x_array <= (float(central) + (float(bandwidth) /2)))] = 1
                       
    #If given in Ghz, transform to Angstrom
    if freq_Ghz:
        wavelength_array = 3e8 * (10 ** 10) / ((x_array)*(10**9) )
        x_array = wavelength_array 
        
    #if given in kev, let central be the lower limit and bandwidth be the upper limit
    if energy_Kev:
        transmission_array = np.zeros(length)
        transmission_array[(x_array >= float(central)) * (x_array <= float(bandwidth))] = 1
        wavelength_array = 12398 / (x_array * 10 ** 3)
        x_array = wavelength_array
        
    plt.plot(x_array, transmission_array)
    plt.xlabel('wavelength in Angstrom')
    plt.ylabel('transmission')
    
    #writing to the file
    save_string = ''
    for i in range(3000):
        save_string += f'{x_array[i]} '
        save_string += f'{transmission_array[i]}\n'
        
    file = open(os.path.join(App.config().get(section='GENERAL', option='data_dir'),'Filterprofiles','TXT',filtername), "w")
    file.write(save_string)
    file.close()
