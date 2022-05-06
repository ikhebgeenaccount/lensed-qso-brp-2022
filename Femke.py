import matplotlib.pyplot as plt
import pandas as pd


from src.lensed_qso import LensedQSO
import numpy as np
from src.mags_to_fluxes import mags_to_fluxes, mag_ratio_split_total_flux
from src.ned_to_sed import ned_table_to_sed
from src.xml_to_txt import xml_to_txt
from src.tophat import tophat
from src.model_subtraction import model_subtraction
from src.plots import plot_lqsos_in_speagle, plot_agnf_output, plot_n_runs_pars, plot_lqsos_vs_stacey, residual_plot, plot_speagle_residual, plot_evolution
from src.percent_to_fraction import percent_to_fraction
from src.filters import populate_filter_profile_path_column
from src.model_sed import fit
import src.ms_data_reader
import os

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
        ax2 = None
        ax3 = None
        lqsos=[]
        for i in range(10):#['J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']:
            
            g=GALAXIES[i]
            print(g)
            lqso = LensedQSO(g)
            lqso.agn_fitter_output(check_git=True, run_time=run_times[i])
            lqsos.append(lqso)
    
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
            
            #plot evolution all galaxies, including the other data
            if ax2 is None:
                fig2, ax2 = plot_evolution(lqso)
            else:
               plot_evolution(lqso, fig=fig2, ax=ax2) 
               
            
            #plot evolution excluding all other data
            if ax3 is None:
                fig3, ax3 = plot_evolution(lqso)
            else:
                plot_evolution(lqso, fig=fig3, ax=ax3) 
               
            #plot_evolution(lqso, single=True)
            
        
        from astropy.cosmology import LambdaCDM
        LCDM = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        for label, df in src.ms_data_reader.FILES.items():
            if 'Birkin' in label:
                continue
            ax2.scatter(LCDM.age(df[df['redshift'] > 0]['redshift']) * 1e9, np.power(10., df[df['redshift'] > 0]['logMstar']), label=label,
                 zorder=50, s=50, alpha=.7)
        ax2.legend()
#        ax2.set_ylim(ymax=1e12)
        ax2.set_xlim(xmin=3.5e9)
        
        fig2.savefig(os.path.join('plots', 'total_evolution_withdata.pdf'))

        #obscuration plot
        fig4, ax4 = plot_agnf_output(lqsos, 'EBVbbb', 'Nh', color_scale_field='SFR_IR', component='_sub')
        # Add Type1/2 AGN separation line as found in AGNfitter paper
        ax4.vlines(0.2, ymin=21.5, ymax=25, color='black', ls='--')
        ax4.hlines(21.5, xmin=0.2, xmax=1, color='black', ls='--')
        fig4.savefig(os.path.join('plots', 'EBVbbb_Nh.pdf'))
            
        plt.show()
    all_galaxies()
    
    

    plt.show()
