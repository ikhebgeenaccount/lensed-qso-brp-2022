#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:47:55 2022

@author: abbo
"""
import glob
import re
import pandas

import numpy as np

# Load the models
MODELS = {}

NAME_MATCH = 'data\/brown_seds\/([A-Za-z0-9_\-+]*).dat'
NAME_MATCH_RE = re.compile(NAME_MATCH)

for sed_file in glob.glob('data/brown_seds/*.dat'):
    name = re.findall(NAME_MATCH_RE, sed_file)[0]
    MODELS[name] = pandas.read_csv(sed_file, delim_whitespace=True, comment='#', names=['wavelength', 'flux', 'observed_wavelength', 'source'])


def fit(lqso):
    sed = lqso.filter_sed()

    avg_scores = []
    model_mults = []

    for model in MODELS:
        mults = []

        for i, row in sed.iterrows():
            closest_wl = closest_wavelength(row.wavelength, model)
            mults.append(row.flux_total / MODELS[model].loc[MODELS[model].wavelength == closest_wl].flux.values[0])

        print(sed.flux_total.shape[0])
        print(MODELS[model].flux.shape[0])
        print(sed.flux_total - MODELS[model].flux)
        avg_scores.append(np.sum(np.power(sed.flux_total - MODELS[model].flux * np.mean(mults), 2.)))
        model_mults.append(mults)

    best_model = list(MODELS.keys())[np.argmin(avg_scores)]
    print(best_model)

    plot_fit(lqso, best_model, model_mults[np.argmin(avg_scores)])


def plot_fit(lqso, model, mults):
    fig, ax = lqso.plot_spectrum(loglog=True)

    for m in mults:
        ax.plot(MODELS[model].wavelength, MODELS[model].flux * m, alpha=0.3)

    ax.plot(MODELS[model].wavelength, MODELS[model].flux * np.mean(m), color='black', alpha=.8)


def closest_wavelength(wl, model):
    return MODELS[model].wavelength.loc[np.argmin(np.abs(MODELS[model].wavelength - wl))]
