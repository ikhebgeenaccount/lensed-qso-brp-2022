import os.path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes
from src.ned_to_sed import ned_table_to_sed

import src.model_sed

GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']

PLOTS_SAVE = 'plots'


def sdss_panstarrs_flux_discrepancy():
    SDSS_V_PANSTARRS_GS = ['B1152+200', 'B1600+434', 'J0806+2006', 'J1633+3134', 'J1455+1447', 'J1524+4409']
    flux_discrepancy = pd.DataFrame(columns=['galaxy', 'filter', 'sdss', 'panstarrs', 'ratio', 'diff'])

    ratios = []
    diffs = []

    panstarrs_mag = 'ap'

    for g in SDSS_V_PANSTARRS_GS:
        try:
            mags_to_fluxes(g)
        except FileNotFoundError:
            continue

        print(f'{g}, PanSTARRS mag choice: {panstarrs_mag}')

        lqso = LensedQSO(g)
        # lqso.plot_spectrum()
        lqso.plot_spectrum(loglog=True)
        lqso.plot_spectrum(mags=True)

        for i, r in lqso.sed.iterrows():
            if 'SDSS' not in r['source'] and panstarrs_mag not in r['source']:
                continue
            if r['filter'] in ['g', 'r', 'i', 'z']:
                source = 'sdss' if 'SDSS' in r.source else 'panstarrs'
                flux_discrepancy.loc[r['filter'], source] = r.flux_total

        flux_discrepancy['ratio'] = flux_discrepancy.sdss.div(flux_discrepancy.panstarrs.where(flux_discrepancy.panstarrs != 0, np.nan))
        # Make all ratios > 1 so we can compare['B1600+434']
        flux_discrepancy['ratio'] = flux_discrepancy['ratio'].apply(lambda r: 1. / r if r < 1  else r)
        flux_discrepancy['diff'] = flux_discrepancy['sdss'] - flux_discrepancy['panstarrs']
        print(flux_discrepancy)
        print()

        ratios.append(list(flux_discrepancy['ratio']))
        diffs.append(list(flux_discrepancy['diff']))

    ratios = np.array(ratios).flatten()
    diffs= np.array(diffs).flatten()

    print(f'ratio avg: {np.nanmean(ratios)}, ratio std: {np.nanstd(ratios)}')
    print(f'diff avg: {np.nanmean(diffs)}, diff std: {np.nanstd(diffs)}')


def all_galaxies():
    for g in GALAXIES:
        lqso = LensedQSO(g)
        lqso.plot_spectrum(loglog=True)

        # mags_to_fluxes(lqso, components=None if g != 'B1608+656' else ['_G', '_G2', '_A', '_B', '_C', '_D', ''])

        count = lqso.filter_sed().loc[lqso.sed['flux_G'] > 0].shape[0]
        print(f'{g}, {count}')

        # telescope = 'Magellan'
        # if lqso.mags.loc[lqso.mags.telescope == telescope].shape[0] > 0:
            # print(g, 'has', telescope)


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
    galaxy = 'J1330+1810'
    lqso = LensedQSO(galaxy)

    # ned_table_to_sed(lqso, ned_file='ned_wise.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])
    # ned_table_to_sed(lqso, ned_file='ned_2mass.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])
    # ned_table_to_sed(lqso, ned_file='ned_chandra.txt', allowed_sources=['Chandra', 'WISE', '2MASS'])

    # mags_to_fluxes(lqso, components=None if galaxy != 'B1608+656' else ['_G', '_G2', '_A', '_B', '_C', '_D', ''])

    lqso.plot_spectrum(loglog=True)


def fit_foreground():
    for g in ['B1600+434']:
        lqso = LensedQSO(g)

        m = 'all' if pd.isnull(lqso.props.lens_type.values[0]) else lqso.props.lens_type.values[0]

        if lqso.filter_sed(component='_G').shape[0] > 1:
             # src.model_sed.fit(lqso, method='minimize', morph=m)
            src.model_sed.fit(lqso, method='curve_fit', morph=m)
        else:
            print(f'{g} has {lqso.filter_sed(component="_G").shape[0]} foreground galaxy datapoints, can\'t be fitted.')


def plot_single_model():
    fig, ax = plt.subplots()

    m = src.model_sed.MODELS['CGCG_453-062_spec']
    ax.plot(m.wavelength, m.flux_cgs)
    ax.plot(m.wavelength, m.flux)


if __name__ == '__main__':
    # all_galaxies()
    fit_foreground()
    # single_galaxy()
    plt.show()
