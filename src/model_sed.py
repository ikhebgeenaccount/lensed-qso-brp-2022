import glob
import os
import re
import pandas

import numpy as np
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


def fit(lqso, morph='all', method='curve_fit', save_plots=True, save_location='plots'):
    """
    Fits a Brown SED to the foreground galaxy data points of given LensedQSO using scipy.optimize.curve_fit.
    :param lqso:
    :param morph: type of allowed morphologies, valid values are 'all', 'spiral', 'elliptical'
    :return:
    """
    sed = lqso.filter_sed(component='_G', allow_zero_error=False).copy()  # Copy such that changing the wavelength doesn't affect the original
    sed['wavelength'] = sed['wavelength'] / (1 + lqso.props['z_lens'].values[0])
    sed.sort_values(by='wavelength')  # Doesn't matter

    model_set_is = MODEL_PROPERTIES.index

    if morph == 'spiral':
        model_set_is = MODEL_PROPERTIES.loc[MODEL_PROPERTIES['morph'].str.contains('S')].index
    elif morph == 'elliptical':
        model_set_is = MODEL_PROPERTIES.loc[MODEL_PROPERTIES['morph'].str.contains('E')].index

    # Arrays to keep track of chi squared values, stds and mults
    chisqs = []
    stds = []
    mults = []

    for i, label_index in enumerate(model_set_is):
        m = MODEL_PROPERTIES.loc[label_index]
        model = m['name']

        if method == 'curve_fit':
            ff = FitFunction(sed, model)

            popt, pcov = curve_fit(ff.f, sed.wavelength.values, sed.flux_G.values, sigma=sed.flux_G_err.values, p0=[1e-2])

            # Keep track of scores and multiplication factors
            chisqs.append(ff.chi_squared(popt))
            stds.append(np.sqrt(np.diag(pcov))[0])
            mults.append(popt[0])
        elif method == 'minimize':
            ff = FitFunction(sed, model)

            res = minimize(ff.chi_squared, [1e-2], method='Powell')

            #model_set['score'].iloc[i] = res.fun
            #model_set['mult'].iloc[i] = res.x  # FIXME: x is value on laptop, array on strw pc, numpy version dependent i guess
            #model_set['std'].iloc[i] = 1  # TODO: determine error, perhaps bootstrap but only very few data points

    # Turn the resulting lists into a DataFrame for easy sorting
    result = {
        'name': MODEL_PROPERTIES['name'].iloc[model_set_is],
        'morph': MODEL_PROPERTIES['morph'].iloc[model_set_is],
        'type': MODEL_PROPERTIES['type'].iloc[model_set_is],
        'score': chisqs,
        'std': stds,
        'mult': mults
    }

    model_set = pandas.DataFrame(result)

    # Lowest score is best model
    model_set.sort_values(by='score', ascending=True, inplace=True)

    print(f'{lqso.name}, best model: {model_set.iloc[[0]]["name"].values[0]}, mult: {model_set.iloc[[0]]["mult"].values[0]:.4e}, std: {model_set.iloc[[0]]["std"].values[0]:.4e}, chisq: {model_set.iloc[[0]]["score"].values[0]:.4e}')
    print(f'Model properties: {model_set.iloc[[0]].morph.values[0]} {model_set.iloc[[0]].type.values[0]}')

    # print(model_set.head(10))

    # Histogram of scores
    fig, ax = plt.subplots()
    ax.hist(model_set['score'].head(25), bins=25)

    plot_fit(lqso, model_set, save_plots=save_plots, save_location=save_location)

    return model_set.iloc[[0]]["name"].values[0], model_set.iloc[[0]]["mult"].values[0], model_set.iloc[[0]]["std"].values[0]


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


def plot_fit(lqso, models, save_plots=True, save_location='plots', count=5):

    # Plot the model on just the foreground galaxy data
    fig, ax = lqso.plot_spectrum(loglog=True, component='_G')
    for i in range(count):
        ax.plot(MODELS[models['name'].iloc[i]].wavelength * (1 + lqso.props['z_lens'].values[0]), MODELS[models['name'].iloc[i]].flux * models['mult'].iloc[i], color='black' if i == 0 else None, alpha=.6 / count * (count - i), label=models['name'].iloc[i])

    ax.legend()

    if save_plots:
        fig.savefig(os.path.join(save_location, lqso.name, 'G_model_fit.pdf'))
        fig.savefig(os.path.join(save_location, lqso.name, 'G_model_fit.jpg'))

    # Plot the model on total flux data
    fig, ax = lqso.plot_spectrum(loglog=True)
    ax.plot(MODELS[models['name'].iloc[0]].wavelength * (1 + lqso.props['z_lens'].values[0]), MODELS[models['name'].iloc[0]].flux * models['mult'].iloc[0], color='black', alpha=.6, label=models['name'].iloc[0])
    ax.legend()

    if save_plots:
        fig.savefig(os.path.join(save_location, lqso.name, 'G_model_fit_full_sed.pdf'))
        fig.savefig(os.path.join(save_location, lqso.name, 'G_model_fit_full_sed.jpg'))


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
