import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes
from src.ned_to_sed import ned_table_to_sed
from src.xml_to_txt import xml_to_txt

if __name__ == '__main__':
    #xml_to_txt('HST.ACS_HRC.F330W.xml', 'HST.ACS_HRC.F330W.txt')
    
    galaxy = 'J1330+1810'
    lqso = LensedQSO(galaxy)
    mags_to_fluxes(lqso)
    #ned_table_to_sed(lqso,'ned_galex', allowed_sources=['Chandra', 'WISE', '2MASS', 'Galex'])
    lqso.plot_spectrum(loglog=True)
    
    plt.show()
