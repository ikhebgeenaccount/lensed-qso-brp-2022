import matplotlib.pyplot as plt
import pandas as pd


from src.lensed_qso import LensedQSO
import numpy as np
from src.mags_to_fluxes import mags_to_fluxes, mag_ratio_split_total_flux
from src.ned_to_sed import ned_table_to_sed
from src.xml_to_txt import xml_to_txt
from src.tophat import tophat
from src.model_subtraction import model_subtraction
from src.plots import plot_lqso_in_speagle
from src.plots import residual_plot, plot_evolution
from src.percent_to_fraction import percent_to_fraction
from src.filters import populate_filter_profile_path_column
from src.model_sed import fit

if __name__ == '__main__':
    
    def single():
        #general
        galaxy = 'B1152+200'  
        lqso = LensedQSO(galaxy)
        lqso.agn_fitter_output(check_git=True, run_time=5)
        
        #lqso.plot_spectrum(loglog=True)#, component='_sub')
        #plot_lqso_in_speagle(lqso)
        #compare_test(lqso)
        #residual_plot(lqso, errors=True)
        plot_evolution(lqso)
        
    #single()
    
    #running all galaxies
    GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']
    PLOTS_SAVE = 'plots'
    cols = ['tau', 'age', 'Nh', 'irlum', 'SB', 'BB', 'GA', 'TO', 'EBVbbb', 'EBVgal', 'logMstar', 'SFR_opt', 'LIR(8-1000)', 'Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)', 'Lsb(1-30)', 'SFR_IR', '-ln_like']
    cols_simple = ['-ln_like']
    run_times = [9,6,5,8,2,1,2,2,0,5]
    
    def all_galaxies():
        ax = None
        ax2=None
        for i in range(10):#['J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']:
            
            g=GALAXIES[i]
            print(g)
            lqso = LensedQSO(g)
            lqso.agn_fitter_output(check_git=True, run_time=run_times[i])
    
            #if '1330' in g:
                #continue
                
            #lqso.plot_spectrum(loglog=True)
            #lqso.plot_error_percentage() #how much of the sub fluxes errors they are in percentages
            #residual_plot(lqso, errors=True)
            
            #speagle
#            if ax is None:
#                 fig, ax = plot_lqso_in_speagle(lqso)
#            else:
#                 plot_lqso_in_speagle(lqso, fig=fig, ax=ax)
            
            #plot evolution all galaxies
            if ax2 is None:
                fig2, ax2 = plot_evolution(lqso)
            else:
               plot_evolution(lqso, fig=fig2, ax=ax2) 
               
            plot_evolution(lqso, single=True)
        plt.show()
    all_galaxies()
    
    

    plt.show()
