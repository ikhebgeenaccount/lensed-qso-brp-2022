# -*- coding: utf-8 -*-
"""
from xml to txt for filter profiles
"""
import os
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

def xml_to_txt(filename):
    """
    The xml files taken from he filter profile service website follow a certain
    structure, this is printed for all the datapoints (wavelength, transmission)
    """
    tree = ET.parse(os.path.join('data','Filterprofiles',filename))
    root = tree.getroot()
    transmission=[]
    wavelength=[]
    for TD in root.iter("TD"):
        if float(TD.text) < 1 :
            transmission.append(float(TD.text))
        else:
            wavelength.append(float(TD.text))
            
    print(wavelength)
    print(transmission)
    plt.plot(wavelength,transmission)
            