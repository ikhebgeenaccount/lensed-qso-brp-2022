import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes
from src.ned_to_sed import ned_table_to_sed
from src.xml_to_txt import xml_to_txt
from src.tophat import tophat

if __name__ == '__main__':
    #xml_to_txt('HST.ACS_HRC.F330W.xml', 'HST.ACS_HRC.F330W.txt')
    
    galaxy = 'B1152+200'
    lqso = LensedQSO(galaxy)
    #mags_to_fluxes(lqso)
    #ned_table_to_sed(lqso,'ned_galex_wise_2mass', allowed_sources=['Chandra', 'WISE', '2MASS', 'Galex'])
    lqso.plot_spectrum(loglog=True, component='_A')
    tophat(6.08,2, 'ATCA_6.08.txt',True)
    
    plt.show()
