import os.path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.agn_fitter_automated import run_agn_fitter
from src.lensed_qso import LensedQSO
from src.latex_output import dataframe_to_latex_table, sed_to_latex_table, plots_in_subfigures,agnf_output_table
from src.mags_to_fluxes import mags_to_fluxes, mag_ratio_split_total_flux
from src.ned_to_sed import ned_table_to_sed
from src.filters import populate_filter_profile_path_column
from src.model_subtraction import model_subtraction
from src.plots import plot_lqsos_in_speagle, plot_agnf_output, plot_n_runs_pars, plot_lqsos_vs_stacey, residual_plot, plot_speagle_residual, hist_stellarmass, plot_lqsos_in_speagle_z_scaled, plot_evolution_df
import src.ms_data_reader

import src.model_sed

import os

plt.style.use('brp.mplstyle')

GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']

PLOTS_SAVE = 'plots'


def _init_lqsos_dict():
    # Create dict for all lqsos
    lqsos_dict = {
        'name': [],
        'redshift': [],
        'magnification': [],
        'magn_err': [],
        'stacey_sfr': [],
        'stacey_sfr_me': [],
        'stacey_sfr_pe': [],
        'mu_m_gas': [],
    }
    # Add lists for all AGNfitter output fields
    for f in src.lensed_qso.AGNFITTER_FIELDS:
        lqsos_dict[f] = []
        lqsos_dict[f'{f}_pe'] = []
        lqsos_dict[f'{f}_me'] = []

        # Also add magnified version
        if f in ['SFR_IR', 'SFR_opt', 'logMstar']:
            lqsos_dict[f'mu_{f}'] = []
            lqsos_dict[f'mu_{f}_pe'] = []
            lqsos_dict[f'mu_{f}_me'] = []

    return lqsos_dict


def _update_lqsos_dict(lqsos_dict, lqso, name=None):
    if name is None:
        lqsos_dict['name'].append(lqso.name)
    else:
        lqsos_dict['name'].append(name)

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

        # Also add magnified version
        if f in ['SFR_IR', 'SFR_opt', 'logMstar']:
            val, pe, me = lqso.get_agnf_output_field(f, demag=False)
            lqsos_dict[f'mu_{f}'].append(val)
            lqsos_dict[f'mu_{f}_pe'].append(pe)
            lqsos_dict[f'mu_{f}_me'].append(me)


def _lqsos_dict_to_df(lqsos_dict):
    # Turn lqsos into a dataframe
    lqsos_df = pd.DataFrame(lqsos_dict)

    # Calculate logSFRs (_IR and _opt) for all lqsos
    lqsos_df['logSFR_IR'] = np.log10(lqsos_df['SFR_IR'])
    lqsos_df['logSFR_IR_pe'] = lqsos_df['SFR_IR_pe'] / (lqsos_df['SFR_IR'] * np.log(10.))
    lqsos_df['logSFR_IR_me'] = lqsos_df['SFR_IR_me'] / (lqsos_df['SFR_IR'] * np.log(10.))

    # Calculate logmu_SFR_IR for comparison with Stacey
    lqsos_df['logmu_SFR_IR'] = np.log10(lqsos_df['mu_SFR_IR'])
    lqsos_df['logmu_SFR_IR_pe'] = lqsos_df['mu_SFR_IR_pe'] / (lqsos_df['mu_SFR_IR'] * np.log(10.))
    lqsos_df['logmu_SFR_IR_me'] = lqsos_df['mu_SFR_IR_me'] / (lqsos_df['mu_SFR_IR'] * np.log(10.))

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

    # Calculate gas fraction
    alpha_CO = 3.6
    lqsos_df['Mgas'] = alpha_CO * lqsos_df['mu_m_gas'] / lqsos_df['magnification']
    lqsos_df['Mgas_err'] = alpha_CO * lqsos_df['mu_m_gas'] / np.square(lqsos_df['magnification']) * lqsos_df['magn_err']
    lqsos_df['f_gas'] = lqsos_df['Mgas'] / (np.power(10., lqsos_df['logMstar']) + lqsos_df['Mgas'])
    lqsos_df['f_gas_pe'] = np.sqrt((np.square(np.power(10., lqsos_df['logMstar']) * lqsos_df['Mgas_err']) +\
                           np.square(lqsos_df['Mgas'] * np.log(10.) * np.power(10., lqsos_df['logMstar']) * lqsos_df['logMstar_pe'])) /\
                           np.power(np.power(10., lqsos_df['logMstar']) + lqsos_df['Mgas'], 4.))
    lqsos_df['f_gas_me'] = np.sqrt((np.square(np.power(10., lqsos_df['logMstar']) * lqsos_df['Mgas_err']) +\
                           np.square(lqsos_df['Mgas'] * np.log(10.) * np.power(10., lqsos_df['logMstar']) * lqsos_df['logMstar_me'])) /\
                           np.power(np.power(10., lqsos_df['logMstar']) + lqsos_df['Mgas'], 4.))
    lqsos_df['t_dep'] = lqsos_df['Mgas'] / lqsos_df['SFR']

    return lqsos_df


def load_all_galaxies(n=10, sub_folder=None, generate_lqso_plots=False, from_file=False):
    if not from_file:
        lqsos_dict = _init_lqsos_dict()

        lqsos_all_runs_df = {}

        for g in GALAXIES:
            lqso = LensedQSO(g)
            print(g)
            # lqso.find_best_run(run_times=n, verbose=True, sub_folder=sub_folder, copy=False)
            # _update_lqsos_dict(lqsos_dict, lqso)

            if generate_lqso_plots:
                lqso.plot_spectrum(disallowed_sources=['chandra', 'luichies'] + src.lensed_qso.FILTERED_SOURCES_AGNFITTER[g])
                lqso.plot_spectrum(component='_sub', disallowed_sources=['chandra', 'luichies'] + src.lensed_qso.FILTERED_SOURCES_AGNFITTER[g])

                # Model subtraction also creates lqso.plot_spectrum but without above filtered sources, so messes up the saved plots
                # model_subtraction(lqso)

                # residual_plot(lqso, errors=True)

                #
                # # Fill DataFrame with every run's output
                # d = _init_lqsos_dict()
                # for i in range(n):
                #     lqso.agn_fitter_output(run_time=i, sub_folder=sub_folder)
                #     _update_lqsos_dict(d, lqso, name=g + str(i))
                # lqsos_all_runs_df[g] = _lqsos_dict_to_df(d)
                #
                # # Plot AGNfitter output stuff
                # plot_n_runs_pars(lqso, sub_folder=sub_folder, n=n)  # when running this one, have to run lqso.find_best_run afterwards again, otherwise stuck on last run
                #
                # # Reload best run
                # lqso.find_best_run(run_times=n, verbose=False, sub_folder=sub_folder)

        # lqsos_df = _lqsos_dict_to_df(lqsos_dict)
        #
        # lqsos_df.to_csv(os.path.join('data', 'final_output.csv'), index=False)

    else:
        lqsos_df = pd.read_csv(os.path.join('data', 'final_output.csv'))
        lqsos_all_runs_df = {}

    return lqsos_df, lqsos_all_runs_df


def generate_context_plots(lqsos_df, lqsos_all_runs_df):

    pd.set_option('display.max_columns', None)
    print(lqsos_df[['name', 'Mgas', 'Mgas_err', 'f_gas', 'f_gas_pe', 'f_gas_me', 't_dep']])
    print(lqsos_df.columns)
    print(np.mean(lqsos_df['f_gas']))
    print(np.mean(lqsos_df['redshift']))

    # Make plots
    plot_lqsos_in_speagle(lqsos_df, label=lqsos_df['name'] + ', z=' + lqsos_df['redshift'].astype(str), group=False)
    plot_lqsos_in_speagle(lqsos_df, label=lqsos_df['name'] + ', z=' + lqsos_df['redshift'].astype(str), group=False, sfr_type='logSFR_opt', save_name='speagle_opt')
    plot_lqsos_in_speagle(lqsos_df, label=lqsos_df['name'] + ', z=' + lqsos_df['redshift'].astype(str), group=False, sfr_type='logSFR_IR', save_name='speagle_IR')

    # Make speagle plot for every galaxy of all runs
    for gal, df in lqsos_all_runs_df.items():
        plot_lqsos_in_speagle(df, label=df['name'] + ',$\ln L$=' + df['-ln_like'].map('{:.1f}'.format), group=False, save_name=f'{gal}_speagle', errorbar_kwargs={'alpha': .6})

    fig, ax = plot_agnf_output(lqsos_df, 'EBVbbb', 'Nh', unique_markers=True)

    plot_agnf_output(lqsos_df, 'SFR_IR', 'SFR_opt', equals_line=True, logx=True, logy=True, unique_markers=True)

    plot_agnf_output(lqsos_df, 'logMstar', 'logSFR', color_scale_field='Lbb(0.1-1)', unique_markers=True)
    plot_agnf_output(lqsos_df, 'logMstar', 'logSFR', color_scale_field='f_gas', unique_markers=True)
    plot_agnf_output(lqsos_df, 'logMstar', 'logSFR', color_scale_field='t_dep', unique_markers=True)
    plot_agnf_output(lqsos_df, 'redshift', 'f_gas', unique_markers=True)

    plot_lqsos_vs_stacey(lqsos_df[lqsos_df['stacey_sfr'] > 0])

    f, a = None, None
    fr, ar = None, None
    fr2, ar2 = None, None
    fh, ah = None, None

    fres, ares = None, None

    # Stellar mass hist figs and axs
    flz, alz = None, None
    flz, alz = hist_stellarmass(lqsos_df, flz, alz, label='This work')
    fhz, ahz = None, None
    fhz, ahz = hist_stellarmass(lqsos_df, fhz, ahz, label='This work')

    # SFR hist
    fhs, ahs = None, None

    fh, ah = hist_stellarmass(lqsos_df, fh, ah, label = 'This work', zorder=10)
    for label, df in src.ms_data_reader.FILES.items():
        if max(df['redshift']) > 0.5:
            fhz, ahz = hist_stellarmass(df, fhz, ahz, label)
        else:
            flz, alz = hist_stellarmass(df, flz, alz, label)

        fh, ah = hist_stellarmass(df, fh, ah, label = label)
        f, a = plot_lqsos_in_speagle(df, label=label, fig=f, ax=a, group=True, errorbar_kwargs={'markersize': 5, 'alpha':.7, 'capsize': 0, 'linewidth': 0}, save_name='speagle_comp')

        fr, ar = plot_speagle_residual(df, label=label, fig=fr, ax=ar, errorbar_kwargs={'markersize': 3, 'alpha':.7, 'capsize': 0, 'linewidth': 0}, save_name='speagle_res')
        fr2, ar2 = plot_speagle_residual(df, label=label, fig=fr2, ax=ar2, x_field='logMstar', x_label='log$M_*$', errorbar_kwargs={'markersize': 3, 'alpha':.7, 'capsize': 0, 'linewidth': 0}, save_name='speagle_res_logMstar')

        fres, ares = plot_lqsos_in_speagle_z_scaled(df, label=label, fig=fres, ax=ares, group=True, errorbar_kwargs={'markersize': 5, 'alpha':.7, 'capsize': 0, 'linewidth': 0}, save_name='speagle_comp_z_scaled')

        if 'SMG' in label:
            fhs, ahs = hist_stellarmass(df, fhs, ahs, label=label, field='logSFR')

    f, a = plot_lqsos_in_speagle(lqsos_df, label='This work', fig=f, ax=a, group=True, errorbar_kwargs={'zorder': 200, 'markersize': 10, 'alpha': 1, 'color': 'black'}, save_name='speagle_comp')
    fres, ares = plot_lqsos_in_speagle_z_scaled(lqsos_df, label='This work', fig=fres, ax=ares, group=True, errorbar_kwargs={'zorder': 200, 'markersize': 10, 'alpha': 1, 'color': 'black'}, save_name='speagle_comp_z_scaled')
    fr, ar = plot_speagle_residual(lqsos_df, label='This work', fig=fr, ax=ar, errorbar_kwargs={'zorder': 200, 'markersize': 10, 'alpha': 1, 'color': 'black'}, save_name='speagle_res')
    fr2, ar2 = plot_speagle_residual(lqsos_df, label='This work', fig=fr2, ax=ar2, x_field='logMstar', x_label='log$M_*$', errorbar_kwargs={'zorder': 200, 'markersize': 10, 'alpha': 1, 'color': 'black'}, save_name='speagle_res_logMstar')

    plot_evolution_df(lqsos_df, context=False)
    plot_evolution_df(lqsos_df)


def big_plot():
    fig, ax = plt.subplots()
    for g in GALAXIES:
        lqso = LensedQSO(g)

        sed = lqso.filter_sed().sort_values(by='wavelength')
        sed = sed.loc[sed.flux_total > 0]

        ax.errorbar(sed.wavelength, sed.flux_total, yerr=sed.flux_err, label=g, alpha=0.7)

    ax.legend()
    ax.set_yscale('log')
    ax.set_xscale('log')

    fig.savefig(os.path.join(PLOTS_SAVE, 'SED_all.pdf'))
    fig.savefig(os.path.join(PLOTS_SAVE, 'SED_all.png'))


def latex(lqsos_df):
    # Filter table to latex
    # print(dataframe_to_latex_table(src.filters.FILTER_PROPERTIES.loc[src.filters.FILTER_PROPERTIES['conversion'].notnull()],
    #                                usecols=['telescope', 'filtername', 'central_wavelength', 'conversion', 'zeropoint'],
    #                                header=['Telescope', 'Filter', 'Central wavelength ($\\unit{\\angstrom}$)', 'Conversion', 'Zeropoint ($\\unit{\milli\jansky}$)'],
    #                                label='table:filter_conv', caption='All filters for which magnitudes were found, with their respective conversion methods and zeropoints.'))

    # Galaxy SEDs to latex
    # for g in GALAXIES:
    #     lqso = LensedQSO(g)

    #     sed_to_latex_table(lqso)

    # SED plots to latex
    # plots_in_subfigures(GALAXIES, 'SED_total', label='sed')

    # Foreground fit plots
    # plots_in_subfigures(GALAXIES, 'G_model_fit', label='fg_fit')

    # AGNfitter residuals
    # plots_in_subfigures(GALAXIES, 'agnf_residuals', label='agnf_res')

    # AGNfitter MS plots 10 runs
    # plots_in_subfigures(GALAXIES, 'speagle', label='10runs')

    # chi^2 plots
    plots_in_subfigures(GALAXIES, 'models_chisq', label='chi_sq')

    # AGNfitter output table
    # print(lqsos_df.columns)
    # agnf_output_table(lqsos_df, cols=['name'] + src.lensed_qso.AGNFITTER_FIELDS[0:10] + ['-ln_like'], label='tab:agnf_output_pars')

    # agnf_output_table(lqsos_df, cols=['name', 'logSFR_IR', 'logSFR_opt', 'logSFR', 'logMstar', 'Mgas', 'f_gas'], label='tab:results_sfrs_ms')


def plot_ell_models():
    fig, ax = plt.subplots()
    models = [3265, '0855', 4621, 4660, 4458]

    for mn in models:
        model = f'NGC_{mn}'
        m = src.model_sed.MODELS[model]

        #ax.plot(m.wavelength, m.flux, label='Brown')
        ax.plot(m.wavelength, np.abs(m.flux), label=model)

    # ax.set_title(model)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()


if __name__ == '__main__':
    lqsos_df, lqsos_all_runs_df = load_all_galaxies(from_file=True, generate_lqso_plots=False)

    generate_context_plots(lqsos_df, lqsos_all_runs_df)
    # plot_ell_models()
    # fit_foreground()
    # fg_subtraction()
    # single_galaxy()
    # known_mag_gals()

    # latex(lqsos_df)

    plt.show()
