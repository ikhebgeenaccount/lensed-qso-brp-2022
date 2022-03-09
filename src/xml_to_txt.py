# -*- coding: utf-8 -*-
"""
from xml to txt for filter profiles
"""
import os
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

def xml_to_txt(xmlname,txtname):
    """
    The xml files taken from he filter profile service website follow the 
    same structure, this fil reads the entries in from the xmlname file and 
    writes them to the txtname file, assuming both the files are in the 
    FIlterprofiles folder. 
    """
    #Reading in the xml file
    tree = ET.parse(os.path.join('data','Filterprofiles',xmlname))
    root = tree.getroot()
    
    #For plotting the filterprofile
    transmission = []
    wavelength = []
    
    #For writing to the new file
    save_string = ''
    
    #Seperate into wavelength and transmission
    for TD in root.iter("TD"):
        if float(TD.text) <= 1 :
            transmission.append(float(TD.text))
            save_string += f' {TD.text}\n'
        else:
            wavelength.append(float(TD.text))
            save_string += f'{TD.text}'
    #plotting to check everything
    plt.plot(wavelength, transmission)
    
    #writing to the new file
    file = open(os.path.join('data','Filterprofiles',txtname), "w")
    file.write(save_string)
    file.close()
    
    

            