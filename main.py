import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes, mag_ratio_split_total_flux
from src.ned_to_sed import ned_table_to_sed
from src.xml_to_txt import xml_to_txt
from src.tophat import tophat
from src.model_subtraction import model_subtraction
from src.percent_to_fraction import percent_to_fraction
from src.filters import populate_filter_profile_path_column
from src.AGN_input import AGN_input

if __name__ == '__main__':
    #filterprofiles
    xml_to_txt('VLT_CONICA_H.xml', 'VLT_CONICA_H.txt')
    #tophat(230.609583,0.560, 'IRAM_1.3mm.txt',freq_Ghz=True, energy_Kev=False)
    #percent_to_fraction('WIYN.U_HARRIS.txt','WIYN.U_HARRIS_fraction.txt')
    
    #photometry
    #galaxy = 'B1608+656'
    #lqso = LensedQSO(galaxy)
    #mags_to_fluxes(lqso)
    #ned_table_to_sed(lqso,'ned_galex_wise_2mass', allowed_sources=['Chandra', 'WISE', '2MASS', 'Galex'])
    #lqso.plot_spectrum(loglog=True)

    #model subtraction
    #model_subtraction(lqso)
    
    
    #AGN_input('J0806+2006')
    
    
    
    GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']
    PLOTS_SAVE = 'plots'

    def all_galaxies():
        for g in GALAXIES:#['J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']:
            lqso = LensedQSO(g)
            lqso.plot_spectrum(loglog=True)
            model_subtraction(lqso)
    
        

    plt.show()
