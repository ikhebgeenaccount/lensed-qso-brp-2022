import os.path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.lensed_qso import LensedQSO
from src.latex_output import dataframe_to_latex_table
from src.mags_to_fluxes import mags_to_fluxes, mag_ratio_split_total_flux
from src.ned_to_sed import ned_table_to_sed
from src.filters import populate_filter_profile_path_column
from src.model_subtraction import model_subtraction
from src.speagle import plot_lqso_in_speagle

import src.model_sed

GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']

PLOTS_SAVE = 'plots'


def all_galaxies():
    ax = None
    for g in GALAXIES:#['J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']:
        lqso = LensedQSO(g)
        lqso.plot_spectrum(loglog=True)

        # mags_to_fluxes(lqso, components=None if g != 'B1608+656' else ['_G', '_G2', '_A', '_B', '_C', '_D', ''])

        count = lqso.filter_sed().loc[lqso.sed['flux_G'] > 0].shape[0]
        print(f'{g}, {count}')

        # telescope = 'Magellan'
        # if lqso.mags.loc[lqso.mags.telescope == telescope].shape[0] > 0:
            # print(g, 'has', telescope)

        model_subtraction(lqso)

        if lqso.agn_fitter_output(copy=False) is not None:
            if ax is None:
                fig, ax = plot_lqso_in_speagle(lqso)
            else:
                plot_lqso_in_speagle(lqso, fig=fig, ax=ax)


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
    galaxy = 'J1455+1447'
    lqso = LensedQSO(galaxy)

    # ned_table_to_sed(lqso, ned_file='ned_wise.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])
    # ned_table_to_sed(lqso, ned_file='ned_2mass.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])
    # ned_table_to_sed(lqso, ned_file='ned_chandra.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])

    # mags_to_fluxes(lqso, components=None if galaxy != 'B1608+656' else ['_G', '_G2', '_A', '_B', '_C', '_D', ''])
    m = 'all' if pd.isnull(lqso.props.lens_type.values[0]) else lqso.props.lens_type.values[0]
    # mag_ratio_split_total_flux(lqso, 'Koopmans+2003', overwrite=True,  components=None if galaxy != 'B1608+656' else ['_G', '_G2', '_A', '_B', '_C', '_D', ''])

    src.model_sed.fit(lqso, m)

    # lqso.plot_spectrum()

    # model_subtraction(lqso)

    #catalog, length = lqso.sed_to_agn_fitter()

    #print(catalog)

    #print(lqso.agn_settings())

    #print(lqso.agn_fitter_output())
    #print(lqso.agn_fitter_output()[['tau', 'age', 'LIR(8-1000)', 'SFR_IR', 'SFR_opt', 'logMstar']])
    plot_lqso_in_speagle(lqso)


def latex():
    dataframe_to_latex_table(src.filters.FILTER_PROPERTIES.loc[src.filters.FILTER_PROPERTIES['conversion'].notnull()],
                                   usecols=['telescope', 'filtername', 'central_wavelength', 'conversion', 'zeropoint'],
                                   header=['Telescope', 'Filter', 'Central wavelength ($\\unit{\\angstrom}$)', 'Conversion', 'Zeropoint ($\\unit{\milli\jansky}$)'],
                                   label='table:filter_conv', caption='All filters for which magnitudes were found, with their respective conversion methods and zeropoints.')


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
    # all_galaxies()
    # plot_ell_models()
    # fit_foreground()
    # fg_subtraction()
    # plot_single_model()
    single_galaxy()

    # latex()

    # populate_filter_profile_path_column()

    plt.show()
