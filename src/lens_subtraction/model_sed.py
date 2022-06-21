import glob
import os
import re
import pandas

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.mplot3d import Axes3D

from scipy.optimize import curve_fit, minimize


# LQSO_NO_MODELS = {
#     'J0806+2006': 2, #very clear
#     'J0924+0219': 2, #very clear
#     'B1152+200': 4, #semi clear
#     'J1330+1810': 5, #not clear, made guess
#     'J1455+1447': 5, #1 datapoint
#     'J1524+4409': 5, #not clear, made guess
#     'B1600+434': 4, #spiral>take best 5, changed to 4 since 4/5 best models have Herschel data
#     'B1608+656': 3, #semi clear  
#     'J1650+4251': 5 #1 datapoint
# }


from src.app import App



# Load model properties
MODEL_PROPERTIES = pandas.read_csv(os.path.join(App.config().get(section='GENERAL', option='data_dir'), 'brown_seds', 'sed_properties.dat'), delimiter='|',
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

for sed_file in glob.glob(os.path.join(App.config().get(section='GENERAL', option='data_dir'), 'brown_seds', '*.dat')):
    try:
        name = re.findall(NAME_MATCH_RE, sed_file)[0]
    except IndexError:
        continue
    MODELS[name] = pandas.read_csv(sed_file, delim_whitespace=True, comment='#', names=['wavelength', 'flux_cgs', 'observed_wavelength', 'source'])
    # Convert units to milliJansky
    MODELS[name]['flux'] = MODELS[name]['flux_cgs'] * np.power(MODELS[name]['wavelength'], 2.) / 2.998e18 * 1e26

    # For the one spiral we have, add Herschel data (flux in millijansky)
    #TODO: this is for only the best fitting models, either add all others or fix the boundary that it cuts off
    if name == 'UGC_12150':
        new_model=MODELS[name].append({'wavelength': 2.4476e6, 'flux' : 5.086e3, 'observed_wavelength':2.5e6,'source':4}, ignore_index = True)
        newest_model=new_model.append({'wavelength': 3.4257e6, 'flux' : 2.031e3, 'observed_wavelength':3.5e6,'source':4}, ignore_index = True)
        newerest_model=newest_model.append({'wavelength': 4.8953e6, 'flux' : 0.611e3, 'observed_wavelength':5e6,'source':4}, ignore_index = True)
        MODELS[name]=newerest_model

    if name == 'NGC_5104':
        new_model=MODELS[name].append({'wavelength': 2.4476e6, 'flux' : 5.266e3, 'observed_wavelength':2.5e6,'source':4}, ignore_index = True)
        newest_model=new_model.append({'wavelength': 3.4257e6, 'flux' : 2.0051e3, 'observed_wavelength':3.5e6,'source':4}, ignore_index = True)
        newerest_model=newest_model.append({'wavelength': 4.8953e6, 'flux' : 0.658e3, 'observed_wavelength':5e6,'source':4}, ignore_index = True)
        MODELS[name]=newerest_model

    if name == 'NGC_5033':
        new_model=MODELS[name].append({'wavelength': 2.4476e6, 'flux' : 40.78e3, 'observed_wavelength':2.5e6,'source':4}, ignore_index = True)
        newest_model=new_model.append({'wavelength': 3.4257e6, 'flux' : 18.221e3, 'observed_wavelength':3.5e6,'source':4}, ignore_index = True)
        newerest_model=newest_model.append({'wavelength': 4.8953e6, 'flux' : 6.474e3, 'observed_wavelength':5e6,'source':4}, ignore_index = True)
        MODELS[name]=newerest_model

    if name == 'NGC_4594':
        new_model=MODELS[name].append({'wavelength': 2.4476e6, 'flux' : 25.6e3, 'observed_wavelength':2.5e6,'source':4}, ignore_index = True)
        newest_model=new_model.append({'wavelength': 3.4257e6, 'flux' : 12.1e3, 'observed_wavelength':3.5e6,'source':4}, ignore_index = True)
        newerest_model=newest_model.append({'wavelength': 4.8953e6, 'flux' : 5.56e3, 'observed_wavelength':5e6,'source':4}, ignore_index = True)
        MODELS[name]=newerest_model


def add_model_data(model_name, wavelength, flux):
    MODELS[model_name].append(pandas.DataFrame({'wavelength': wavelength, 'flux': flux}))



# Fitting
#
# Duncan+2017 (https://academic.oup.com/mnras/article/473/2/2655/4315948?login=true) uses EAZY
# (https://github.com/gbrammer/eazy-photoz/tree/master) on single template mode, EAZY is for photometric redshift though
#
# Leja+2017 (https://iopscience.iop.org/article/10.3847/1538-4357/aa5ffe/meta) uses scipy minimize combined with MCMC


def fit(lqso, morph='all', method='curve_fit', save_plots=True, save_location='plots', verbose_plots=False, N=5):

    """
    Fits a Brown SED to the foreground galaxy data points of given LensedQSO using scipy.optimize.curve_fit.
    :param verbose_plots:
    :param save_plots:
    :param method:
    :param lqso:
    :param morph: type of allowed morphologies, valid values are 'all', 'spiral', 'elliptical'
    :return: three 1D arrays, wavelengths, fluxes, errors for each flux at each wavelength
    """
    sed = lqso.filter_sed(component='_G', allow_zero_error=lqso.name == 'J1650+4251').copy()  # Copy such that changing the wavelength doesn't affect the original
    sed['wavelength'] = sed['wavelength'] / (1 + lqso.props['z_lens'].values[0])
    sed.sort_values(by='wavelength')  # Doesn't matter

    if sed.flux_G.values.shape[0] == 0:
        print(f'Not enough foreground data points to fit {lqso.name}.')
        return -1, -1, -1

    model_set_is = MODEL_PROPERTIES.index

    if morph == 'spiral':
        model_set_is = MODEL_PROPERTIES.loc[MODEL_PROPERTIES['morph'].str.contains('S')].index
    elif morph == 'elliptical':
        model_set_is = MODEL_PROPERTIES.loc[MODEL_PROPERTIES['morph'].str.contains('E')].index

    # Arrays to keep track of chi squared values, stds and mults
    chisqs = []
    stds = []
    mults = []

    single_g = False        
    if lqso.name in ['J1650+4251', 'J1455+1447']:
        # Fit only selection of models
        models = ['NGC_3265', 'NGC_0855', 'NGC_4621', 'NGC_4660', 'NGC_4458']

        model_set_is = MODEL_PROPERTIES.loc[MODEL_PROPERTIES['name'].isin(models)].index

        single_g = True
        
        #TODO: fix this by detecting when there is a single datapoint

    elif lqso.name == 'B1608+656':
        # Fit only selection of models
        models = ['NGC_3265', 'NGC_4660', 'NGC_4458']

        model_set_is = MODEL_PROPERTIES.loc[MODEL_PROPERTIES['name'].isin(models)].index

        single_g = True
        #TODO: fix when only a selection of models, negative fluxes

    elif lqso.name == 'B1600+434':
        # Only fit good fitting spiral models with Herschel data
        models = ['UGC_12150', 'NGC_5104', 'NGC_5033', 'NGC_4594']

        model_set_is = MODEL_PROPERTIES.loc[MODEL_PROPERTIES['name'].isin(models)].index
        #TODO: fix when only spirals

    for i, label_index in enumerate(model_set_is):
        m = MODEL_PROPERTIES.loc[label_index]
        model = m['name']

        if method == 'curve_fit':
            ff = FitFunction(sed, model)
            popt, pcov = curve_fit(ff.f, sed.wavelength.values, sed.flux_G.values, sigma=None if lqso.name == 'J1650+4251' else sed.flux_G_err.values, p0=[1e-2])

            # Keep track of scores and multiplication factors
            chisqs.append(ff.chi_squared(popt))
            stds.append(np.sqrt(np.diag(pcov))[0])
            mults.append(popt[0])
        elif method == 'minimize':
            ff = FitFunction(sed, model)

            res = minimize(ff.chi_squared, [1e-3], method='Powell')


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
    # fig, ax = plt.subplots()
    # ax.hist(model_set['score'], bins=25)

    model_set['red_chi_sq'] = model_set['score'] / (sed.flux_G.shape[0] - 1)

    fig, ax = plt.subplots()
    x = range(model_set.shape[0])
    ax.scatter(x, model_set['red_chi_sq'])
    ax.set_xticks(x, model_set['name'].values, rotation=90)
    ax.set_ylabel('$\chi^2_\mathrm{red}$')
    # ax.set_title(f'Reduced $\chi^2$ values of models for {lqso.name}_G')

    if save_plots:
        fig.savefig(os.path.join(App.config().get(section='GENERAL', option='plots_dir'), f'{lqso.name}_models_chisq.pdf'), bbox_inches='tight')
        # fig.savefig(os.path.join(save_location, f'{lqso.name}_models_chisq.jpg'))

    print(model_set[['name', 'red_chi_sq', 'mult', 'std']].head(15))

    if model_set['red_chi_sq'].iloc[0] / model_set['red_chi_sq'].iloc[1] <= 0.5:
        print('Best model is twice as good as next best')

    # Combine N models
    #N = LQSO_NO_MODELS[lqso.name]



    if N != 0:

        models_wl = [list(MODELS[row['name']]['wavelength'].values) for i, row in model_set.head(N).iterrows()]
        models_flux = [MODELS[row['name']]['flux'].values * row['mult'] for i, row in model_set.head(N).iterrows()]

        # Models are not guaranteed to be same length
        # Take highest starting point (highest np.min of models) as starting point
        # Take lowest end point as end point
        start_wl = np.max([np.min(ws) for ws in models_wl])
        end_wl = np.min([np.max(ws) for ws in models_wl])

        all_wls = np.array([w for wls in models_wl for w in wls])

        # Get all the unique wavelenghts that the models have, so we can interp at those wavelenghts
        wls = np.unique(all_wls[(all_wls >= start_wl) & (all_wls <= end_wl)])
        wls = np.sort(wls)

        norm_weights = (1. / model_set['red_chi_sq'].head(N).values / np.sum(1. / model_set['red_chi_sq'].head(N).values)).reshape((N, 1))

        if np.isnan(norm_weights).any():
            norm_weights = 1

        interp_models_flux = np.stack([np.interp(wls, models_wl[i], models_flux[i]) for i in range(N)])
        avg_model = np.average(interp_models_flux, axis=0, weights=norm_weights.reshape(N) if not single_g else None)

        stds = model_set['std'].head(N).values.reshape((N, 1))

        # mult_err_prop = np.linalg.norm(interp_models_flux * stds, axis=0)
        models_std = np.sqrt(np.sum(norm_weights * np.power(interp_models_flux - avg_model, 2.), axis=0))

        avg_model_err = np.sqrt(np.sum(np.power(norm_weights * stds * interp_models_flux, 2.), axis=0))

        if np.isnan(avg_model_err).any() or np.isinf(avg_model_err).any():
            avg_model_err = 0

        avg_model_errs = np.sqrt(np.power(avg_model_err, 2.) + np.power(models_std, 2.))

        print(np.sum(np.abs(models_std - avg_model_errs)))

        print('Avg model red chi sq:', FitFunction(sed, wls=wls, fluxs=avg_model).chi_squared([1]) / (sed.flux_G.shape[0] - 1))

    if verbose_plots:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        fig = plt.figure()
        moax = fig.add_subplot(111, projection='3d')

        models_flux = [MODELS[row['name']]['flux'].values * row['mult'] for _, row in model_set.iterrows()]
        models_wl = [list(MODELS[row['name']]['wavelength'].values) for _, row in model_set.iterrows()]

        all_wls = np.array([w for wls in models_wl for w in wls])

        # Get all the unique wavelenghts that the models have, so we can interp at those wavelenghts
        mwls = np.unique(all_wls[(all_wls >= all_wls) & (all_wls <= all_wls)])
        mwls = np.sort(mwls)
        interp_models_flux = np.stack([np.interp(mwls, models_wl[j], models_flux[j]) for j in range(model_set.shape[0])])

        for i in range(model_set.shape[0]):
            t_model = np.average(interp_models_flux[0:i + 1], axis=0, weights=1. / model_set['red_chi_sq'].head(i + 1).values)

            moax.plot3D([i] * len(mwls), np.log10(mwls), np.log10(interp_models_flux[i]))
            ax.plot3D([i] * len(mwls), np.log10(mwls), np.log10(t_model))

        ax.set_xlabel('Model count')
        ax.set_ylabel('$\mathit{Wavelength}\ (\mathrm{\AA})$')
        ax.set_zlabel('$\mathit{Flux\ density}\ (\mathrm{mJy})$')

    if N == 0:
        print(f'N = 0 for {lqso.name}.')
        return None

    wls = wls * (1. + lqso.props.z_lens.values[0])

    plot_fit(lqso, model_set, avg_model=(wls, avg_model, avg_model_errs), save_plots=save_plots, count=N)

    return wls, avg_model, avg_model_errs


class FitFunction:

    def __init__(self, sed, model=None, wls=None, fluxs=None, interp=True, component='_G'):
        """
        Class with function to fit.
        :param sed:
        :param model:
        :param interp:
        :param component:
        """
        self.sed = sed
        self.interp = interp
        self.component = component
        if model is not None:
            self.model = model
            self.wls = MODELS[model]['wavelength'].values
            self.fluxs = MODELS[model]['flux'].values
        elif wls is not None and fluxs is not None:
            self.wls = wls
            self.fluxs = fluxs
        else:
            raise ValueError('Either model or both wls and fluxs must be given.')

    # f is the function that will be used for curve_fit
    def f(self, x, mult):
        if self.interp:
            return mult * interp_fluxes(x, wls=self.wls, fluxs=self.fluxs)
        else:
            model = self.model
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


def plot_fit(lqso, models, avg_model, save_plots=True, save_location='plots', count=5):
    avg_wls, avg_model, avg_err = avg_model

    # Plot the model on just the foreground galaxy data
    fig, ax, labels, legend_list = lqso.plot_spectrum(loglog=True, component='_G')
    for i in range(count):
        labels.append(models['name'].iloc[i])
        legend_list.append(ax.plot(MODELS[models['name'].iloc[i]].wavelength * (1 + lqso.props['z_lens'].values[0]), MODELS[models['name'].iloc[i]].flux * models['mult'].iloc[i], alpha=.6 / count * (count - i) + .4, label=models['name'].iloc[i])[0])

    # labels.append('Average')
    # legend_list.append(ax.plot(avg_wls, avg_model, label='Average', color='black')[0])
    avg = (ax.fill_between(avg_wls, avg_model - avg_err, avg_model + avg_err, color='grey', alpha=.5), ax.plot(avg_wls, avg_model, color='black')[0])
    legend_list.append(avg)
    labels.append('Avg. model $\pm 1\sigma$')

    ax.legend(legend_list, labels, handler_map={avg: HandlerTuple(), tuple: HandlerTuple(ndivide=None)})

    if save_plots:
        fig.savefig(os.path.join(App.config().get(section='GENERAL', option='plots_dir'), f'{lqso.name}_G_model_fit.pdf'))

    rfig, rax = fig, ax

    # Plot the model on total flux data
    fig, ax, labels, legend_list = lqso.plot_spectrum(loglog=True)
    # labels.append('Average')
    avg = (ax.fill_between(avg_wls, avg_model - avg_err, avg_model + avg_err, color='grey', alpha=.5), ax.plot(avg_wls, avg_model, color='black')[0])
    legend_list.append(avg)
    labels.append('Avg. model $\pm 1\sigma$')

    # TODO: just as a test for now, remove later
    #if lqso.name == 'B1600+434':
        # These numbers are flux densities for wavelenghts that are longer than model, since lensing galaxy is a spiral it can have radio contribution
        #rax.scatter([60e4, 2.14137e9], [8.17, 38.3], label='radio model', color='fuchsia')
        #rax.legend()
        #ax.scatter([60e4, 2.14137e9], [8.17, 38.3], label='radio model', color='fuchsia')

    ax.legend(legend_list, labels, handler_map={avg: HandlerTuple(), tuple: HandlerTuple(ndivide=None)})

    if save_plots:
        fig.savefig(os.path.join(App.config().get(section='GENERAL', option='plots_dir'), f'{lqso.name}_G_model_fit_full_sed.pdf'))


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


def interp_fluxes(wl, model=None, wls=None, fluxs=None):
    """
    If wls given, must also have fluxs given. model can be name of model in MODELS dict.
    """
    if wls is not None and fluxs is not None:
        return np.interp(wl, xp=wls, fp=fluxs)
    elif model is not None:
        return np.interp(wl, xp=MODELS[model].wavelength.values, fp=MODELS[model].flux.values)
    else:
        raise ValueError('Either model, or both wls and fluxs must be given.')
