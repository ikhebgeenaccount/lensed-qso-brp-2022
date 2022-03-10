import glob
import os
import re
import pandas

import numpy as np
import tqdm as tqdm
from matplotlib import pyplot as plt

from scipy.optimize import curve_fit, minimize

# Load model properties
MODEL_PROPERTIES = pandas.read_csv(os.path.join('data', 'brown_seds', 'sed_properties.dat'), delimiter='|',
                               names=['name', 'ra_dec', 'morph', 'type', 'bpt_class'], usecols=[0, 1, 2, 3, 4])

# Remove all leading and trainling spaces
MODEL_PROPERTIES['name'] = MODEL_PROPERTIES['name'].apply(lambda n: n.strip().replace(' ', '_'))  # Replace spaces with underscores such that names are equal for model spec files
MODEL_PROPERTIES['ra_dec'] = MODEL_PROPERTIES['ra_dec'].apply(lambda n: n.strip())
MODEL_PROPERTIES['morph'] = MODEL_PROPERTIES['morph'].apply(lambda n: n.strip())
MODEL_PROPERTIES['type'] = MODEL_PROPERTIES['type'].apply(lambda n: n.strip())

# Load the models
MODELS = {}

NAME_MATCH = '([A-Za-z0-9_\-+]*)_spec.dat'
NAME_MATCH_RE = re.compile(NAME_MATCH)

for sed_file in glob.glob(os.path.join('data', 'brown_seds', '*.dat')):
    try:
        name = re.findall(NAME_MATCH_RE, sed_file)[0]
    except IndexError:
        continue
    MODELS[name] = pandas.read_csv(sed_file, delim_whitespace=True, comment='#', names=['wavelength', 'flux_cgs', 'observed_wavelength', 'source'])
    # Convert units to Jansky
    MODELS[name]['flux'] = MODELS[name]['flux_cgs'] * np.power(MODELS[name]['wavelength'], 2.) / 2.998e18 * 1e26

# Fitting
#
# Duncan+2017 (https://academic.oup.com/mnras/article/473/2/2655/4315948?login=true) uses EAZY
# (https://github.com/gbrammer/eazy-photoz/tree/master) on single template mode, EAZY is for photometric redshift though
#
# Leja+2017 (https://iopscience.iop.org/article/10.3847/1538-4357/aa5ffe/meta) uses scipy minimize combined with MCMC


def fit(lqso, morph=None, method='curve_fit'):
    """
    Fits a Brown SED to the foreground galaxy data points of given LensedQSO using scipy.optimize.curve_fit.
    :param lqso:
    :param morph: List of allowed morphologies
    :return:
    """
    sed = lqso.filter_sed(component='_G')
    sed.sort_values(by='wavelength')  # Doesn't matter

    # Arrays to store scores and multiplication factors per model
    scores = []
    covs = []
    model_mults = []

    if morph is None:
        model_set = MODEL_PROPERTIES
    else:
        model_set = MODEL_PROPERTIES.loc[MODEL_PROPERTIES['morph'].isin(morph)]

    # Chi squared necessities
    # Parameter space of mults
    mults_space = np.logspace(-7, -2, 10000)  # logspace, 10**-7 through 10**-2
    # Arrays to save chisq and mults values for plots
    chisqs = []

    for i in tqdm.tqdm(range(0, model_set.shape[0])):
        m = model_set.loc[[i]]
        model = m['name'].values[0]

        if method == 'curve_fit':
            ff = FitFunction(sed, model)

            popt, pcov = curve_fit(ff.f, sed.wavelength.values, sed.flux_G.values, sigma=sed.flux_G_err.values, p0=[1e-2])

            # Keep track of scores and multiplication factors
            # Score is sum of squared difference between fitted function f and data of foreground galaxy divided by flux errors
            scores.append(ff.chi_squared(popt[0]))
            covs.append(np.sqrt(np.diag(pcov))[0])
            model_mults.append(popt[0])
        elif method == 'minimize':
            ff = FitFunction(sed, model)

            res = minimize(ff.chi_squared, [1e-2], method='Powell')

            scores.append(res.fun)
            model_mults.append(res.x)
            covs.append(1)  # WONTFIX: we don't use minimize

        elif method == 'chi_squared':
            ff = FitFunction(sed, model)

            chisqs.append([])

            # TODO: optimize loop away
            for m in mults_space:
                chisqs[i].append(ff.chi_squared([m]))

            l_index = np.argmin(chisqs[i])

            scores.append(chisqs[i][l_index])
            model_mults.append(mults_space[l_index])

            covs.append(1)  # TODO: bootstrap?

    # Lowest score is best model
    bm_index = np.argmin(scores)
    best_model = model_set.loc[[bm_index]]['name'].values[0]  # FIXME: when using a subset of the models, this will not point to the correct one, I think
    best_mult = model_mults[bm_index]

    print(f'{lqso.name}, best model: {best_model}, mult: {best_mult:.4e}, std: {covs[bm_index]:.4e}, chisq: {scores[bm_index]}')

    if method == 'chi_squared':
        fig, ax = plt.subplots()
        ax.plot(mults_space, chisqs[bm_index])
        ax.set_xscale('log')
        ax.set_yscale('log')

    # Histogram of scores
    # fig, ax = plt.subplots()
    # ax.hist(scores, bins=20)

    plot_fit(lqso, best_model, best_mult)

    return best_model, best_mult


class FitFunction:

    def __init__(self, sed, model, interp=True, component='_G'):
        """
        Class with function to fit.
        :param sed:
        :param model:
        :param interp:
        :param component:
        """
        self.sed = sed
        self.model = model
        self.interp = interp
        self.component = component

    # f is the function that will be used for curve_fit
    def f(self, x, mult):
        model = self.model
        if self.interp:
            return mult * interp_fluxes(x, model)
        else:
            cwl = closest_wavelength(x, model)
            return mult * MODELS[model].loc[MODELS[model].wavelength.isin(cwl)].flux

    def chi_squared(self, mults):
        mult = mults[0]
        # TODO: calculate proper model wavelengths based on lens redshift

        # print(mult)
        # print('interp values', self.f(self.sed.wavelength.values, mult))
        # print('true values', self.sed[f'flux{self.component}'].values)
        # print('subtracted', self.f(self.sed.wavelength.values, mult) - self.sed[f'flux{self.component}'].values)
        # print('errors', self.sed[f'flux{self.component}_err'].values)
        # print('div', (self.f(self.sed.wavelength.values, mult) - self.sed[f'flux{self.component}'].values) /
        #                        self.sed[f'flux{self.component}_err'].values)
        return np.sum(np.power((self.f(self.sed.wavelength.values, mult) - self.sed[f'flux{self.component}'].values) /
                               self.sed[f'flux{self.component}_err'].values, 2))


def plot_fit(lqso, model, mults):
    # Plot the model on just the foreground galaxy data
    fig, ax = lqso.plot_spectrum(loglog=True, component='_G')
    ax.plot(MODELS[model].wavelength, MODELS[model].flux * mults, color='black', alpha=.6, label=model)
    ax.legend()

    # Plot the model on total flux data
    fig, ax = lqso.plot_spectrum(loglog=True)
    ax.plot(MODELS[model].wavelength, MODELS[model].flux * mults, color='black', alpha=.6, label=model)
    ax.legend()


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
    # print(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size) - wl.reshape(wl.size, 1))
    # Find the argmin of the abs after subtraction, i.e. the closest value of model wavelength for each given
    # wavelength, thus axis=1
    # print(np.argmin(np.abs(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size) - wl.reshape(wl.size, 1)), axis=1))
    # Use these argmin as indices for the model wavelenghts to find the closest wavelenghts to the given wavelenghts.
    # print(MODELS[model].wavelength.loc[np.argmin(np.abs(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size) - wl.reshape(wl.size, 1)), axis=1)].values)
    return MODELS[model].wavelength.loc[np.argmin(
        np.abs(MODELS[model].wavelength.values.reshape(1, MODELS[model].wavelength.values.size) -
               wl.reshape(wl.size, 1)), axis=1)].values


def interp_fluxes(wl, model):
    return np.interp(wl, xp=MODELS[model].wavelength.values, fp=MODELS[model].flux.values)
