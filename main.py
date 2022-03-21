import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes
from src.ned_to_sed import ned_table_to_sed
from src.xml_to_txt import xml_to_txt
from src.tophat import tophat
from src.model_subtraction import model_subtraction
from src.percent_to_fraction import percent_to_fraction
from src.filters import populate_filter_profile_path_column

if __name__ == '__main__':
    #filterprofiles
    #xml_to_txt('WIYN.B_HARRIS.xml', 'WIYN.B_HARRIS.txt')
    #tophat(0.2,0.5, 'CHANDRA_0.2_0.5.txt',freq_Ghz=False, energy_Kev=True)
    #percent_to_fraction('WIYN.U_HARRIS.txt','WIYN.U_HARRIS_fraction.txt')
    
    #photometry
    galaxy = 'B1152+200'
    lqso = LensedQSO(galaxy)
    #mags_to_fluxes(lqso)
    #ned_table_to_sed(lqso,'ned_galex_wise_2mass', allowed_sources=['Chandra', 'WISE', '2MASS', 'Galex'])
    lqso.plot_spectrum(loglog=True)
    
    #lqso.sed['telescope']=0
    #lqso.save_sed()

    model_subtraction(lqso)

    
    
    plt.show()
