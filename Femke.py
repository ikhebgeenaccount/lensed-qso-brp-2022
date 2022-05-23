import matplotlib.pyplot as plt
import pandas as pd


from src.lensed_qso import LensedQSO
import numpy as np
from src.mags_to_fluxes import mags_to_fluxes, mag_ratio_split_total_flux
from src.ned_to_sed import ned_table_to_sed
from src.xml_to_txt import xml_to_txt
from src.tophat import tophat
from src.model_subtraction import model_subtraction
from src.plots import plot_lqsos_in_speagle, plot_agnf_output, plot_n_runs_pars, plot_lqsos_vs_stacey, residual_plot, plot_speagle_residual, plot_evolution, hist_stellarmass, plot_evolution_df
from src.percent_to_fraction import percent_to_fraction
from src.filters import populate_filter_profile_path_column
from src.model_sed import fit
import src.ms_data_reader
import os




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
    

    
#running all galaxies
GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']
PLOTS_SAVE = 'plots'
cols = ['tau', 'age', 'Nh', 'irlum', 'SB', 'BB', 'GA', 'TO', 'EBVbbb', 'EBVgal', 'logMstar', 'SFR_opt', 'LIR(8-1000)', 'Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)', 'Lsb(1-30)', 'SFR_IR', '-ln_like']
cols_simple = ['-ln_like']
run_times = [2,1,1,0,8,9,2,9,2,6]
    
def all_galaxies():
    ax = None
    ax2 = None
    ax3 = None
    lqsos = []
    
    # Create dict for all lqsos
    lqsos_dict = {
        'name': [],
        'redshift': [],
        'magnification': [],
        'magn_err': [],
        'stacey_sfr': [],
        'stacey_sfr_me': [],
        'stacey_sfr_pe': [],
        'mu_m_gas':[],
    }
    
# Add lists for all AGNfitter output fields
    for f in src.lensed_qso.AGNFITTER_FIELDS:
        lqsos_dict[f] = []
        lqsos_dict[f'{f}_pe'] = []
        lqsos_dict[f'{f}_me'] = []

    #Loop over all galaxies
    for i in range(10):#['J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']:
        
        g=GALAXIES[i]
        print(g)
        #find the latest output
        lqso = LensedQSO(g)
        lqso.agn_fitter_output(check_git=True, run_time=run_times[i])
        lqsos.append(lqso)
        #model_subtraction(lqso)
        #lqso.plot_spectrum(loglog=True)
        #lqso.plot_error_percentage() #how much of the sub fluxes errors they are in percentages
        #residual_plot(lqso, errors=True)
        
        


          
        #plot_evolution(lqso, single=True)
        
        #Make the dictionary
        lqsos_dict['name'].append(lqso.name)
        lqsos_dict['redshift'].append(lqso.props['z_qso'].values[0])
        lqsos_dict['magnification'].append(lqso.props['magnification'].values[0])
        lqsos_dict['magn_err'].append(lqso.props['magn_err'].values[0])
        lqsos_dict['stacey_sfr'].append(lqso.props['stacey_sfr'].values[0])
        lqsos_dict['stacey_sfr_me'].append(lqso.props['stacey_sfr_me'].values[0])
        lqsos_dict['stacey_sfr_pe'].append(lqso.props['stacey_sfr_pe'].values[0])
        lqsos_dict['mu_m_gas'].append(lqso.props['mu_m_gas'].values[0])
        
        for f in src.lensed_qso.AGNFITTER_FIELDS:
        # We set demag=True for every field, since demag only happens for logMstar, SFR_IR, SFR_opt
            val, pe, me = lqso.get_agnf_output_field(f, demag=True)
            lqsos_dict[f].append(val)
            lqsos_dict[f'{f}_pe'].append(pe)
            lqsos_dict[f'{f}_me'].append(me)

            
    
    #make the dataframe        
    lqsos_df = pd.DataFrame(lqsos_dict)
            
    # Calculate logSFRs (_IR and _opt) for all lqsos
    lqsos_df['logSFR_IR'] = np.log10(lqsos_df['SFR_IR'])
    lqsos_df['logSFR_IR_pe'] = lqsos_df['SFR_IR_pe'] / (lqsos_df['SFR_IR'] * np.log(10.))
    lqsos_df['logSFR_IR_me'] = lqsos_df['SFR_IR_me'] / (lqsos_df['SFR_IR'] * np.log(10.))

    lqsos_df['logSFR_opt'] = np.log10(lqsos_df['SFR_opt'])
    lqsos_df['logSFR_opt_pe'] = lqsos_df['SFR_opt_pe'] / (lqsos_df['SFR_opt'] * np.log(10.))
    lqsos_df['logSFR_opt_me'] = lqsos_df['SFR_opt_me'] / (lqsos_df['SFR_opt'] * np.log(10.))

    # Calculate total logSFR
    lqsos_df['SFR'] = lqsos_df['SFR_IR'] + lqsos_df['SFR_opt']
    lqsos_df['SFR_pe'] = np.sqrt(np.square(lqsos_df['SFR_IR_pe']) + np.square(lqsos_df['SFR_opt_pe']))
    lqsos_df['SFR_me'] = np.sqrt(np.square(lqsos_df['SFR_IR_me']) + np.square(lqsos_df['SFR_opt_me']))

    lqsos_df['logSFR'] = np.log10(lqsos_df['SFR'])
    lqsos_df['logSFR_pe'] = lqsos_df['SFR_pe'] / (lqsos_df['SFR'] * np.log(10.))
    lqsos_df['logSFR_me'] = lqsos_df['SFR_me'] / (lqsos_df['SFR'] * np.log(10.))
    
    lqsos_df['Mgas'] = lqsos_df['mu_m_gas'] / lqsos_df['magnification']
    lqsos_df['Mgas_err'] = lqsos_df['mu_m_gas'] / np.square(lqsos_df['magnification']) * lqsos_df['magn_err']
    lqsos_df['f_gas'] = lqsos_df['Mgas'] / (np.power(10., lqsos_df['logMstar']) + lqsos_df['Mgas'])
    lqsos_df['f_gas_pe'] = np.sqrt((np.square(np.power(10., lqsos_df['logMstar']) * lqsos_df['Mgas_err']) +\
                           np.square(lqsos_df['Mgas'] * np.log(10.) * np.power(10., lqsos_df['logMstar']) * lqsos_df['logMstar_pe'])) /\
                           np.power(np.power(10., lqsos_df['logMstar']) + lqsos_df['Mgas'], 4.))
    lqsos_df['f_gas_me'] = np.sqrt((np.square(np.power(10., lqsos_df['logMstar']) * lqsos_df['Mgas_err']) +\
                           np.square(lqsos_df['Mgas'] * np.log(10.) * np.power(10., lqsos_df['logMstar']) * lqsos_df['logMstar_me'])) /\
                           np.power(np.power(10., lqsos_df['logMstar']) + lqsos_df['Mgas'], 4.))
    print(lqsos_df['f_gas'])
    # Make plots
    plot_lqsos_in_speagle(lqsos_df, label=lqsos_df['name'], group=False)
    plot_lqsos_in_speagle(lqsos_df, label=lqsos_df['name'], group=False, sfr_type='logSFR_opt', save_name='speagle_opt')
    plot_lqsos_in_speagle(lqsos_df, label=lqsos_df['name'], group=False, sfr_type='logSFR_IR', save_name='speagle_IR')

    # Make speagle plot for every galaxy of all runs
    # for gal, df in lqsos_all_runs_df.items():
    #     plot_lqsos_in_speagle(df, label=df['name'], group=False, save_name=f'{gal}_speagle', errorbar_kwargs={'alpha': .6})

#    
    #obscuration plot# Add Type1/2 AGN separation line as found in AGNfitter paper
#    fig, ax = plot_agnf_output(lqsos, 'EBVbbb', 'Nh', color_scale_field='SFR_IR', component='_sub')
#    ax.vlines(0.2, ymin=21.5, ymax=25, color='black', ls='--')
#    ax.hlines(21.5, xmin=0.2, xmax=1, color='black', ls='--')
#    fig.savefig(os.path.join('plots', 'EBVbbb_Nh.pdf'))
    
    #plot_lqsos_vs_stacey(lqsos)
    #plot_agnf_output(lqsos, 'SFR_IR', 'SFR_opt', color_scale_field='log age', equals_line=True, logx=True, logy=True)
    
    f , a =plt.subplots(figsize=(10,8))
    f4, a4 = None, None
    hist_stellarmass(lqsos_df, f, a, label = 'our sample', zorder=10)
    for label, df in src.ms_data_reader.FILES.items():
        hist_stellarmass(df, f, a, label = label)
        f4, a4 = plot_lqsos_in_speagle(df, label=label, fig=f4, ax=a4, group=True, errorbar_kwargs={'markersize': 3, 'alpha':.5}, save_name='speagle_comp')
    f4, a4 = plot_lqsos_in_speagle(lqsos_df, label='This work', fig=f4, ax=a4, group=True, errorbar_kwargs={'zorder': 200, 'markersize': 10, 'alpha': 1, 'color': 'black'}, save_name='speagle_comp')
    a.legend()
    a.set_xlabel('Log M_star', fontsize=14)
    a.set_ylabel('normalised number density', fontsize=14)
    
    #plot evolution all galaxies, including the other data
    fig2, ax2 = plot_evolution_df(lqsos_df)
    

    
    #plot evolution excluding all other data

    fig3, ax3 = plot_evolution_df(lqsos_df)
    
    #Adding the other data to evolution plot
    from astropy.cosmology import LambdaCDM
    LCDM = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    for label, df in src.ms_data_reader.FILES.items():
        if 'Birkin' in label:
            continue
        if 'Sun' in label:
            continue
        ax2.scatter(LCDM.age(df[df['redshift'] > 0]['redshift']) * 1e9, np.power(10., df[df['redshift'] > 0]['logMstar']), label=label,
            zorder=50, s=50, alpha=.7)
    ax2.legend()
#   ax2.set_ylim(ymax=1e12)
    #ax2.set_xlim(xmin=3.5e9)
    
    fig2.savefig(os.path.join('plots', 'total_evolution_withdata.pdf'))
    
    
if __name__ == '__main__':
    #single()
    all_galaxies()
    
    
    
    
    
    
    
    



