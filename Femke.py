import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
import numpy as np
from src.mags_to_fluxes import mags_to_fluxes, mag_ratio_split_total_flux
from src.ned_to_sed import ned_table_to_sed
from src.xml_to_txt import xml_to_txt
from src.tophat import tophat
from src.model_subtraction import model_subtraction
from src.plots import plot_lqso_in_speagle
from src.percent_to_fraction import percent_to_fraction
from src.filters import populate_filter_profile_path_column

if __name__ == '__main__':
    
    def single():
        #general
        galaxy = 'J0806+2006' 
        lqso = LensedQSO(galaxy)
        
        #photometry
        #ned_table_to_sed(lqso,'ned_galex_wise_2mass', allowed_sources=['Chandra', 'WISE', '2MASS', 'Galex'])
        lqso.plot_spectrum(loglog=True, component='_sub')
        
        #filterprofiles
        #xml_to_txt('VLT_CONICA_H.xml', 'VLT_CONICA_H.txt')
        #tophat(230.609583,0.560, 'IRAM_1.3mm.txt',freq_Ghz=True, energy_Kev=False)
    
        #model subtraction
        #model_subtraction(lqso)
        
        #AGN input
        #print(lqso.sed_to_agn_fitter(rX=False))
        #print(lqso.agn_settings(rX=False))
        
        #speagle
        #plot_lqso_in_speagle(lqso)
        
        #agn_output
        #compare_test(lqso)
    #single()
    
    #running all galaxies
    GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']
    PLOTS_SAVE = 'plots'
    cols = ['tau', 'age', 'Nh', 'irlum', 'SB', 'BB', 'GA', 'TO', 'EBVbbb', 'EBVgal', 'logMstar', 'SFR_opt', 'LIR(8-1000)', 'Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)', 'Lsb(1-30)', 'SFR_IR', '-ln_like']
    cols_simple = ['SFR_IR']
    
    def all_galaxies():
        ax = None
        for g in GALAXIES:#['J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']:
            #general
            lqso = LensedQSO(g)
            #lqso.plot_spectrum(loglog=True)
            lqso.find_best_run(run_times=10, verbose=False)#, sub_folder='5runs_100nwalkers')
            
            #AGnfitter checking
            print(g)       
            for variable in cols_simple:
                    sub_output = lqso.get_agnf_output_field(variable)[0]
                    print(f'best value of {variable}', np.log10(sub_output))
            
            #lqso.plot_error_percentage() #how much of the sub fluxes errors they are in percentages
            
            #speagle
            # if lqso.agn_fitter_output(copy=True) is not None:
            #     if ax is None:
            #          fig, ax = plot_lqso_in_speagle(lqso)
            #     else:
            #          plot_lqso_in_speagle(lqso, fig=fig, ax=ax)
            
    
    all_galaxies()

    plt.show()
