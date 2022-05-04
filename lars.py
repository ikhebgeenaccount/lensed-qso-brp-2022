import os.path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.agn_fitter_automated import run_agn_fitter
from src.lensed_qso import LensedQSO
from src.latex_output import dataframe_to_latex_table, sed_to_latex_table, all_seds_plot, plots_in_subfigures
from src.mags_to_fluxes import mags_to_fluxes, mag_ratio_split_total_flux
from src.ned_to_sed import ned_table_to_sed
from src.filters import populate_filter_profile_path_column
from src.model_subtraction import model_subtraction
from src.plots import plot_lqso_in_speagle, plot_agnf_output, plot_n_runs_pars, plot_lqsos_vs_stacey, residual_plot
import src.ms_data_reader

import src.model_sed

import os

plt.style.use('brp.mplstyle')

GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']

PLOTS_SAVE = 'plots'


def all_galaxies(n=10, sub_folder=None):
    lqsos = []
    fig = None
    ax = None
    figm, axm = None, None
    lqsos = []
    for g in GALAXIES:
        lqso = LensedQSO(g)
        lqsos.append(lqso)
        print(g)
        lqso.find_best_run(run_times=n, verbose=True, sub_folder=sub_folder, copy=False)

        # Plot in speagle but grouped
        figm, axm = plot_lqso_in_speagle(lqso, figm, axm, label='This work', save_name='speagle_comp.pdf',
                                         errorbar_kwargs={'color': 'black', 'zorder': 100, 'markersize': 6, 'alpha': .7})

        # Plot in standard speagle
        # fig, ax = plot_lqso_in_speagle(lqso, fig=fig, ax=ax)

        # Plot AGNfitter output stuff
        # plot_n_runs_pars(lqso, sub_folder=sub_folder, n=n)
        # residual_plot(lqso, errors=True)

        # figs, axs = None, None
        # for i in range(n):
        #     lqso.agn_fitter_output(run_time=i, sub_folder=sub_folder)
        #     figs, axs = plot_lqso_in_speagle(lqso, figs, axs, label=lqso.name + str(i), save_name=f'{lqso.name}_speagle')

        # lqso.find_best_run(run_times=n, verbose=False, sub_folder=sub_folder)

    fig, ax = plot_agnf_output(lqsos, 'EBVbbb', 'Nh', color_scale_field='SFR_IR', component='_sub')
    # Add Type1/2 AGN separation line as found in AGNfitter paper
    ax.vlines(0.2, ymin=21.5, ymax=25, color='black', ls='--')
    ax.hlines(21.5, xmin=0.2, xmax=1, color='black', ls='--')
    fig.savefig(os.path.join('plots', f'EBVbbb_Nh.pdf'))

    plot_agnf_output(lqsos, 'SFR_IR', 'SFR_opt', color_scale_field='age', equals_line=True, logx=True, logy=True)
    plot_lqsos_vs_stacey(lqsos)

    for label, df in src.ms_data_reader.FILES.items():
        axm.errorbar(df['logMstar'], df['logSFR'],
                     xerr=None if np.sum([df['logMstar_me'], df['logMstar_pe']]) == 0 else
                     np.reshape([df['logMstar_me'], df['logMstar_pe']], (2, len(df['logMstar_pe']))),
                     yerr=None if np.sum([df['logSFR_me'], df['logSFR_pe']]) == 0 else
                     np.reshape([df['logSFR_me'], df['logSFR_pe']], (2, len(df['logSFR_pe']))), label=label, fmt='o',
                     zorder=50, markersize=6, alpha=.7)

    axm.legend([axm.get_legend().get_texts()[0].get_text(), 'This work'] + list(src.ms_data_reader.FILES.keys()))
    figm.savefig(os.path.join('plots', 'speagle_comp.pdf'))


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


def single_galaxy():
    galaxy = 'J1650+4251'
    lqso = LensedQSO(galaxy)
    # lqso.find_best_run(run_times=5)

    # print(lqso.filter_sed(disallowed_sources=src.lensed_qso.FILTERED_SOURCES_AGNFITTER[lqso.name]).sort_values(by='wavelength')[['source', 'wavelength']])
    # print(lqso.sed_to_agn_fitter(component='_sub_demag'))

    # ned_table_to_sed(lqso, ned_file='ned_wise.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])
    # ned_table_to_sed(lqso, ned_file='ned_2mass.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])
    # ned_table_to_sed(lqso, ned_file='ned_chandra.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])

    # mags_to_fluxes(lqso, components=None if galaxy != 'B1608+656' else ['_G', '_G2', '_A', '_B', '_C', '_D', ''])
    # m = 'all' if pd.isnull(lqso.props.lens_type.values[0]) else lqso.props.lens_type.values[0]
    # mag_ratio_split_total_flux(lqso, 'Koopmans+2003', overwrite=True,  components=None if galaxy != 'B1608+656' else ['_G', '_G2', '_A', '_B', '_C', '_D', ''])

    # a = src.model_sed.fit(lqso, m)

    model_subtraction(lqso)

    lqso.plot_spectrum()
    # lqso.plot_spectrum(component='_sub')
    # lqso.plot_spectrum(component='_sub_demag')

    # catalog, length = lqso.sed_to_agn_fitter(rX=True)

    # print(catalog)
    # print(length)

    # lqso.agn_settings(rX=True)

    # lqso.load_agnf_output()
    # print(lqso.get_agnf_output_field('SFR_IR', component='_sub'))
    # print(lqso.get_agnf_output_field('SFR_IR', component='_sub_demag_test'))
    #print(lqso.agn_fitter_output()[['tau', 'age', 'LIR(8-1000)', 'SFR_IR', 'SFR_opt', 'logMstar']])
    # plot_lqso_in_speagle(lqso)

    # plot_agnf_output([galaxy], 'SFR_IR', 'SFR_opt', color_scale_field='age')

    # fig, ax = plot_agnf_output([galaxy], 'EBVbbb', 'Nh', color_scale_field='SFR_IR')

    # Add Type1/2 AGN separation line as found in AGNfitter paper
    # ax.vlines(0.2, ymin=21.5, ymax=25, color='black', ls='--')
    # ax.hlines(21.5, xmin=0.2, xmax=1, color='black', ls='--')


def known_mag_gals():
    known_mag_gals = ['J1524+4409', 'B1608+656' , 'J1455+1447']

    fig, ax = None, None
    for g in known_mag_gals:
        lqso = LensedQSO(g)
        fig, ax = plot_lqso_in_speagle(lqso, fig=fig, ax=ax)

    plot_agnf_output(known_mag_gals, 'SFR_IR', 'SFR_opt', color_scale_field='age')

    fig, ax = plot_agnf_output(known_mag_gals, 'EBVbbb', 'Nh', color_scale_field='SFR_IR')

    # Add Type1/2 AGN separation line as found in AGNfitter paper
    ax.vlines(0.2, ymin=21.5, ymax=25, color='black', ls='--')
    ax.hlines(21.5, xmin=0.2, xmax=1, color='black', ls='--')


def latex():
    # Filter table to latex
    # dataframe_to_latex_table(src.filters.FILTER_PROPERTIES.loc[src.filters.FILTER_PROPERTIES['conversion'].notnull()],
    #                                usecols=['telescope', 'filtername', 'central_wavelength', 'conversion', 'zeropoint'],
    #                                header=['Telescope', 'Filter', 'Central wavelength ($\\unit{\\angstrom}$)', 'Conversion', 'Zeropoint ($\\unit{\milli\jansky}$)'],
    #                                label='table:filter_conv', caption='All filters for which magnitudes were found, with their respective conversion methods and zeropoints.')

    # Galaxy SEDs to latex
    # for g in GALAXIES:
    #     lqso = LensedQSO(g)

    #     sed_to_latex_table(lqso)

    # SED plots to latex
    # all_seds_plot(GALAXIES)

    # Foreground fit plots
    plots_in_subfigures(GALAXIES, 'G_model_fit', 'Fits to foreground galaxy datapoints for every galaxy.')


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
    all_galaxies()
    # plot_ell_models()
    # fit_foreground()
    # fg_subtraction()
    # single_galaxy()
    # known_mag_gals()

    # latex()

    plt.show()
