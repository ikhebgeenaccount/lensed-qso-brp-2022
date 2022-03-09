import glob
import os
import re
import pandas

import numpy as np
from matplotlib import pyplot as plt

from scipy.optimize import curve_fit

# Load the models
MODELS = {}

NAME_MATCH = '([A-Za-z0-9_\-+]*).dat'
NAME_MATCH_RE = re.compile(NAME_MATCH)

for sed_file in glob.glob(os.path.join('data', 'brown_seds', '*.dat')):
    name = re.findall(NAME_MATCH_RE, sed_file)[0]
    MODELS[name] = pandas.read_csv(sed_file, delim_whitespace=True, comment='#', names=['wavelength', 'flux', 'observed_wavelength', 'source'])


def fit(lqso):
    """
    Fits a Brown SED to the foreground galaxy data points of given LensedQSO.
    :param lqso:
    :return:
    """
    sed = lqso.filter_sed(component='_G')
    sed.sort_values(by='wavelength')  # Doesn't matter

    # Arrays to store scores and multiplication factors per model
    scores = []
    model_mults = []

    for model in MODELS:

        # f is the function that will be used for curve_fit
        def f(x, mult):
            cwl = closest_wavelength(x, model)
            return mult * MODELS[model].loc[MODELS[model].wavelength.isin(cwl)].flux

        popt, pcov = curve_fit(f, sed.wavelength.values, sed.flux_G.values, sigma=sed.flux_G_err.values, p0=[1e12])

        # Keep track of scores and multiplication factors
        # Score is sum of absolute difference between fitted function f and data of foreground galaxy
        scores.append(np.sum(np.abs(f(sed.wavelength.values, popt[0]) - sed.flux_G.values)))
        model_mults.append(popt[0])

    # Lowest score is best model
    best_model = list(MODELS.keys())[np.argmin(scores)]
    best_mult = model_mults[np.argmin(scores)]
    print(best_model)
    print(best_mult)

    # Histogram of scores
    fig, ax = plt.subplots()
    ax.hist(scores, bins=20)

    plot_fit(lqso, best_model, best_mult)

    return best_model, best_mult


def plot_fit(lqso, model, mults):
    # Plot the model on just the foreground galaxy data
    fig, ax = lqso.plot_spectrum(loglog=True, component='_G')
    ax.plot(MODELS[model].wavelength, MODELS[model].flux * mults, color='black', alpha=.8)

    # Plot the model on total flux data
    fig, ax = lqso.plot_spectrum(loglog=True)
    ax.plot(MODELS[model].wavelength, MODELS[model].flux * mults, color='black', alpha=.8)


def closest_wavelength(wl, model):
    """
    Given an array of wavelenghts wl, returns an array containing the closest wavelenghts that are present in the model
    :param wl:
    :param model:
    :return:
    """
    # print(wl)
    # Reshape such that wl becomes a 2D array
    # print(wl.reshape(wl.size, 1))
    # Reshape the model wavelengths to fill one row
    # print(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size))
    # Now, subtracting causes the model wavelength array to be broadcast and have extra rows added, while the wl array
    # will add columns to conform to the model wavelength array shape. Then, each given wl wavelenght is subtracted from
    # each value of model wavelength.
    # print(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size) - wl.reshape(wl.size, 1))
    # Find the argmin of the abs after subtraction, i.e. the closest value of model wavelength for each given
    # wavelength, thus axis=1
    # print(np.argmin(np.abs(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size) - wl.reshape(wl.size, 1)), axis=1))
    # Use these argmin as indices for the model wavelenghts to find the closest wavelenghts to the given wavelenghts.
    # print(MODELS[model].wavelength.loc[np.argmin(np.abs(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size) - wl.reshape(wl.size, 1)), axis=1)].values)
    return MODELS[model].wavelength.loc[np.argmin(
        np.abs(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size) -
               wl.reshape(wl.size, 1)), axis=1)].values
