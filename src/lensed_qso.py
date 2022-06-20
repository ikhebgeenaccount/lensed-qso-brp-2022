import os
import distutils.dir_util
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.legend_handler import HandlerTuple

from src.agnfitter.AGN_input import format_filter_name, format_telescope_name, get_agnf_filter_path
from src.app import App
from src.sed_compilation.filters import get_wavelength, FILTER_PROPERTIES

RADIO_CUTOFF = App.config().getfloat(section='GENERAL', option='radio_cutoff')  # wavelengths >2.5e7 Angstrom are classified as radio
XRAY_CUTOFF = App.config().getfloat(section='GENERAL', option='xray_cutoff')  # wavelengths < 300 Angstrom are classified as Xray

# Dictionary that tracks colors of sources used in .plot_spectrum,
# to make sure colours are consistent between plots of different LensedQSOs
SOURCES_COLORS = {}

COLORS = plt.get_cmap('tab10').colors + plt.get_cmap('Dark2').colors + plt.get_cmap('Accent').colors

# TODO: what is AGNfitter rX-version output?
AGNFITTER_FIELDS = ['tau', 'log age', 'Nh', 'irlum', 'SB', 'BB', 'GA', 'TO', 'EBVbbb', 'EBVgal', 'logMstar', 'SFR_opt', 'LIR(8-1000)', 'Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)', 'Lsb(1-30)', 'SFR_IR', '-ln_like']

PROPERTIES = pd.read_csv(os.path.join(App.config().get(section='GENERAL', option='data_dir'), 'properties.csv'), skiprows=0)


class LensedQSO:

    def __init__(self, name, sed_source='sed.csv',  mags_source='mags.csv', save_all_plots=True, save_location='plots'):
        """
        :param name: Name of the lensed qso.
        :param sed_source: source of the spectral energy density (SED). Must be located in 'data/[name]'. Default is
        'data/[name]/sed.csv'.
        """
        self.name = name

        # Save location for plots
        self.save_all_plots = save_all_plots
        self.save_location = os.path.join(save_location)

        # Check if location exists
        # Base plots folder
        if not os.path.isdir(save_location):
            os.mkdir(save_location)
        # This specific galaxy
        if not os.path.isdir(self.save_location):
            os.mkdir(self.save_location)

        # Read SED
        self.sed_file = sed_source
        self.sed = pd.read_csv(os.path.join(App.config().get(section='GENERAL', option='data_dir'), name, sed_source))
        self.sed.fillna(0., inplace=True)

        try:
            # Read mags
            self.mags = pd.read_csv(os.path.join(App.config().get(section='GENERAL', option='data_dir'), name, mags_source))
            self.mags.fillna(0., inplace=True)
        except FileNotFoundError:
            print('No mags file found for galaxy ' + self.name)

        self.props = PROPERTIES.loc[PROPERTIES.galaxy == self.name]

        # filtered_sed only selects entries that have a wavelength and a flux_total
        self.filtered_sed = self.sed[(self.sed.wavelength > 0) & (self.sed.flux_total > 0)].copy()

        # self.load_agnf_output()

    def filter_sed(self, disallowed_sources=None, component='_total', allow_zero_error=True, rX=True):
        """
        Returns a DataFrame that contains only those rows that are allowed through LensedQSO.allowed_sources
        :param disallowed_sources: passed to LensedQSO.allowed_sources
        :return:
        """
        allowed_sources = self.allowed_sources(disallowed_sources)

        compfilter = np.ones(self.sed.shape[0], dtype=bool)
        if component is not None:
            compfilter = self.sed[f'flux{component}'] > 0

        rX_filter = np.ones(self.sed.shape[0], dtype=bool)
        if not rX:
            rX_filter = (self.sed['wavelength'] >= XRAY_CUTOFF) & (self.sed['wavelength'] <= RADIO_CUTOFF)

        errorfilter = np.ones(self.sed.shape[0], dtype=bool)
        if not allow_zero_error:
            errorfilter= self.sed[f'flux{"" if component == "_total" else component}_err'] > 0

        return self.sed.loc[self.sed.source.isin(allowed_sources) * compfilter * errorfilter * rX_filter]

    def allowed_sources(self, disallowed_sources=None):
        """
        Returns a list of allowed sources based on FILTERED_SOURCES and disallowed_sources.
        :param disallowed_sources:
        :return:
        """
        if disallowed_sources is None:
            disallowed_sources = []

        if App.config().has_option(section='FILTERED_SOURCES', option=self.name):
            # Add standard filtered sources to disallowed sources
            disallowed_sources += App.config().getlist(section='FILTERED_SOURCES', option=self.name)

            # No need to warn about PanSTARRS, otherwise do warn
            if len(disallowed_sources) > 0:
                pass # warnings.warn(self.name + ': Filtering sources ' + str(disallowed_sources))

        u_sources = list(self.sed.source.unique())
        to_remove = []

        # Check for occurrences of disallowed sources in the unique sources
        for ds in disallowed_sources:
            for s in u_sources:
                if ds.lower() in s.lower() and s not in to_remove:
                    to_remove.append(s)

        # Remove all found disallowed sources
        for r in to_remove:
            u_sources.remove(r)

        return u_sources

    def plot_spectrum(self, loglog=True, mags=False, disallowed_sources=None, component=None, **kwargs):
        if mags:
            raise NotImplementedError('I broke mags plot so if you want it let me know.')

        # Fill with NaNs in case something was added
        self.sed.fillna(0, inplace=True)

        fig, ax = plt.subplots(figsize=(10, 8))
        legend_list = []

        # upper_limits = ()

        data = None
        data_type = None
        data_err = None

        if mags:
            data = self.mags
            data_type = 'mag'
            data_err = 'mag_err'
            limit = 'lower_limit'

            def wl(row):
                return get_wavelength(row.telescope, row['filter'])

            data['wavelength'] = data.apply(lambda row: wl(row), axis=1)

        else:
            data = self.sed.copy()
            data['source'] = data['source'].str.replace('_filter', '', regex=False)
            data_type = 'flux_total' if component is None else f'flux{component}'
            data_err = 'flux_err' if component is None else f'flux{component}_err'
            limit = 'upper_limit'

        # For every unique source, add their data separately
        u_sources = self.allowed_sources(disallowed_sources)

        plotted_sources = []

        for l in u_sources:
            # Filter based on source
            # Only take those that have both a wavelength and a total flux
            sel = data[(data.source == l) & (data.wavelength > 0) & (data[data_type] > 0)]

            # If there are no entries after filtering, continue to next source
            if len(sel) == 0:
                continue

            plotted_sources.append(l)

            # Separate upper limits from regular data points
            sel_upper_limit = sel[sel[limit] == 1]
            sel_reg = sel[sel[limit] == 0]

            if l in SOURCES_COLORS:
                color = SOURCES_COLORS[l]
            else:
                color = COLORS[len(SOURCES_COLORS) + 1]
                SOURCES_COLORS[l] = color

            # Plot regular data points and upper limits separately, upper limits with special marker
            le_1, _, _ = ax.errorbar(sel_reg.wavelength, sel_reg[data_type], sel_reg[data_err], fmt='o', label=l, color=color, **kwargs)

            if len(sel_upper_limit) > 0:
                # le_2, _, _ = ax.errorbar(sel_upper_limit.wavelength, sel_upper_limit[data_type], sel_upper_limit[data_err],
                #                          fmt='v', label=l, color=le_1.get_color())
                le_2 = ax.scatter(sel_upper_limit.wavelength, sel_upper_limit[data_type],# sel_upper_limit[data_err],
                                         marker='v', label=l, color=le_1.get_color())
                legend_list.append((le_1, le_2))
            else:
                legend_list.append(le_1)

            # upper_limits += (le_2, )

        # ax.legend(legend_list + [upper_limits], list(self.sed.source.unique()) + ['upper limit'], loc='upper left', handler_map={tuple: HandlerTuple(ndivide=None)})
        ax.legend(legend_list, plotted_sources, handler_map={tuple: HandlerTuple(ndivide=None)})
        ax.set_xscale('log')

        if loglog:
            ax.set_yscale('log')

        # ax.set_title(f'{self.name} SED' if component is None else f'{self.name}{component} SED')
        ax.set_xlabel('$\mathit{Wavelength}\ [\mathrm{\AA}]$')
        if mags:
            ax.set_ylabel('$\mathit{mag}$')
        else:
            ax.set_ylabel('$\mathit{Flux\ density}\ [\mathrm{mJy}]$')

        fig.tight_layout()
        if self.save_all_plots:
            fig.savefig(os.path.join(App.config().get(section='GENERAL', option='plots_dir'), f'{self.name}_SED{component if component is not None else "_total"}.pdf'))
            # fig.savefig(os.path.join(self.save_location, f'{self.name}_SED{component if component is not None else "_total"}.png'))

        return fig, ax, plotted_sources, legend_list

    def plot_error_percentage(self, loglog=True):
        fig, ax = plt.subplots(figsize=(10, 8))
        ratio=self.sed['flux_sub_err']/self.sed['flux_sub'] *100
        ax.scatter(self.sed['wavelength'][ratio<=100],ratio[ratio<=100], label='sub')
        ratio2=self.sed['flux_err']/self.sed['flux_total'] *100
        ax.scatter(self.sed['wavelength'][ratio2<=100],ratio2[ratio2<=100], label='total')
        ax.set_xscale('log')
        ax.set_ylabel('errors of sub in percentage', fontsize=15)
        ax.set_title(f'{self.name} percentage error')

        return fig,ax

    def save_sed(self):
        self.sed.to_csv(os.path.join(App.config().get(section='GENERAL', option='data_dir'), self.name, self.sed_file), index=False)

    def sed_to_agn_fitter(self, rX=False, component='_sub', run_times=1):
        """
        Translates the SED to a catalog that is compatible with AGNfitter.
        :return: str
        """
        id = self.agn_fitter_id(component=component)

        header = '# ID redshift [wavelength_angstrom flux_mJy flux_error_mJy]\n'

        l = 2
        catalog = header
        for j in range(max(2, run_times)):
            catalog_line = f'{str(id) + ("" if j == 0 else str(j))} {self.props.z_qso.values[0]} '
            for i, row in self.filter_sed(component=component, rX=rX,
                                          disallowed_sources=App.config().getlist(section='FILTERED_SOURCES_AGNFITTER', option=self.name)
                                          if App.config().has_option(section='FILTERED_SOURCES_AGNFITTER', option=self.name) else None).iterrows():

                if not row.upper_limit:
                    catalog_line += f'{row.wavelength} {row[f"flux{component}"]} {row[f"flux{component}_err"]} '
                else:
                    # Upper limit has error -99 as flag for AGNfitter
                    catalog_line += f'{row.wavelength} {row[f"flux{component}"]} -99 '

                if j== 0:
                    l += 3
            catalog += catalog_line + '\n'

        return catalog, l

    def agn_fitter_id(self):
        t = self.name.replace('B', '').replace('J', '').replace('+', '')
        return t if t[0] != '0' else t[1:]

    def agn_fitter_settings(self, rX=False):
        # First load default settings to make sure every setting is set
        settings_use = dict()

        if rX:
            with open(App.config().get(section='AGNFITTER-RX', option='template_path_rX'), 'r') as template_file:
                template = template_file.read()

            filters, filternames, filterfilenames = self.agn_fitter_input_filters(rX)

            hasxray = np.sum(self.filter_sed(component='_sub', rX=True)['wavelength'] <= XRAY_CUTOFF) > 0
            hasradio = np.sum(self.filter_sed(component='_sub', rX=True)['wavelength'] >= RADIO_CUTOFF) > 0
            settings_use.update({
                'length': self.sed_to_agn_fitter()[1],
                'length + 1': self.sed_to_agn_fitter()[1] + 1,
                'length + 2': self.sed_to_agn_fitter()[1] + 2,
                'name': self.name,
                'redshift': self.props.z_qso.values[0],
                'filters': filters,
                'filternames': filternames,
                'filterfilenames': filterfilenames,
                'hasxray': hasxray,
                'hasradio': hasradio
            })
        else:
            with open(App.config().get(section='AGNFITTER', option='template_path'), 'r') as template_file:
                template = template_file.read()
            settings_use.update({
                'length': self.sed_to_agn_fitter()[1],
                'length + 1': self.sed_to_agn_fitter()[1] + 1,
                'length + 2': self.sed_to_agn_fitter()[1] + 2,
                'name': self.name,
                'redshift': self.props.z_qso.values[0],
                'filters': self.agn_fitter_input_filters()
            })

        # Update default settings with settings that have been set by user
        settings_use.update(dict(App.config()['AGNFITTER']))
        if rX:
            settings_use.update(dict(App.config()['AGNFITTER-RX']))

        return template.format(**settings_use)

    def agn_fitter_input_filters(self, rX=False):
        """
        Returns string of filters for regular AGNfitter settings.
        Returns triplet of (string of filters, list of filternames, list of filterfilenames) for rX AGNfitter version
        """
        outstr = ''

        # rX reqs
        # filters = self.filter_sed(component='_sub')['filter'].values
        filternames = []
        filterfilenames = []

        # Iterate over each row to get the dictionary settings
        for i, row in FILTER_PROPERTIES.iterrows():

            tel = FILTER_PROPERTIES.telescope.values[i]
            fil = FILTER_PROPERTIES.filtername.values[i]

            fsed = self.filter_sed(component='_sub', rX=rX,
                                   disallowed_sources=App.config().getlist(section='FILTERED_SOURCES_AGNFITTER', option=self.name)
                                   if App.config().has_option(section='FILTERED_SOURCES_AGNFITTER', option=self.name) else None)

            exists_in_sed = fsed.loc[(fsed['telescope'] == tel) & (fsed['filter'] == fil)].shape[0] > 0

            if rX:
                filterfilenames += [get_agnf_filter_path(tel, fil)]
                filternames += [fil]

            outstr += f"    filters['_{format_telescope_name(tel)}_{format_filter_name(fil)}']={exists_in_sed}\n"

        if rX:
            return outstr, filternames, filterfilenames
        else:
            return outstr

    def load_agnf_output(self, rX=False, copy=False):
        self.agn_fitter_output(rX=rX, copy=copy, check_git=False)

    def find_best_run(self, run_times=1, rX=False, component='_sub', verbose=True, copy=False, sub_folder=None):
        """
        Finds the best AGNfitter run, based on the log likelihood. Highest log likelihood corresponds to best fit.
        Loads that best AGNfitter run as output.
        :param run_times:
        :param rX:
        :param component:
        :return:
        """
        lls = []
        index = []
        self.agnf_output = [0]

        if verbose:
            print('Run\t-ll')

        for i in range(run_times):
            output = self.agn_fitter_output(rX=rX, agnf_id=self.agn_fitter_id() + str(i if i != 0 else ''), sub_folder=sub_folder)

            if output is None:
                continue

            lls.append(output['-ln_like'].values[2])
            index.append(i)

            if verbose:
                print(f'{i}\t{lls[-1]}')

        if verbose:
            print(f'Best run: {np.argmax(lls)}')

        return self.agn_fitter_output(rX=rX, agnf_id=self.agn_fitter_id() + str(index[np.argmax(lls)] if index[np.argmax(lls)] != 0 else ''), copy=copy, sub_folder=sub_folder)

    def agn_fitter_output(self, rX=False, agnf_id=None, copy=False, check_git=False, run_time=0, sub_folder=''):
        # TODO: agnfitter output reading
        if agnf_id is None:
            agnf_id = self.agn_fitter_id() + (str(run_time) if run_time != 0 else '')
        if rX:
            raise NotImplementedError('rX path not yet')
        else:
            if sub_folder is None:
                path = os.path.join(os.pardir, 'AGNfitter', 'OUTPUT', str(agnf_id))
            else:
                path = os.path.join(os.pardir, 'AGNfitter', 'OUTPUT', sub_folder, str(agnf_id))
        repo_path = os.path.join(App.config().get(section='GENERAL', option='data_dir'), f'{self.name}', 'agnfitter')

        # Three cases:
        # path exists: points to AGNfitter output
        # repo_path exists: points to data in repo
        # neither exist: sucks
        if os.path.isdir(path):
            par_values_file = f'parameter_outvalues_{agnf_id}.txt'
            if copy:
                distutils.dir_util.copy_tree(path, os.path.join(App.config().get(section='GENERAL', option='data_dir'), self.name, 'agnfitter'))
        elif os.path.isdir(repo_path) and check_git:
            # TODO: fix
            # print(re.escape(os.path.join(repo_path, f'parameter_outvalues_{self.agn_fitter_id(component=component)}')).replace('/', '\\/') + '[0-9]?\.txt')
            # par_values_file = glob.glob(re.escape(os.path.join(repo_path, f'parameter_outvalues_{self.agn_fitter_id(component=component)}')).replace('/', '\\/') + '[0-9]?\.txt')
            # print(par_values_file)
            par_values_file = f'parameter_outvalues_{agnf_id}.txt'
            path = repo_path
        else:
            # print('Are you working on vdesk or strw? If not, then this output has not been copied to our github repo, or doesn\'t exist')
            return

        # Read agnfitter sed_data and sed_realizations output
        try:
            self.agnf_sed_data = pd.read_csv(os.path.join(path, 'sed_data.csv'))
            self.agnf_sed_realizations = pd.read_csv(os.path.join(path, 'sed_realizations.csv'))
        except:
            print('No SED output found in AGNfitter output for' + self.name)

        # Read parameter outvalues
        output = pd.read_csv(os.path.join(path, par_values_file),
                             delim_whitespace=True, skiprows=4, header=None, names=AGNFITTER_FIELDS)

        if not hasattr(self, 'agnf_output'):
            self.agnf_output = [0]

        self.agnf_output[cid] = {}
        for c in AGNFITTER_FIELDS:
            self.agnf_output[c] = []

            self.agnf_output[c].append(output[c].iloc[2])
            self.agnf_output[c].append(output[c].iloc[3] - self.agnf_output[c][0])
            self.agnf_output[c].append(self.agnf_output[c][0] - output[c].iloc[1])

        return output

    def get_agnf_output_field(self, field, demag, component='_sub'):
        if hasattr(self, 'agnf_output'):

            if demag and component == '_sub':
                # Save values in easy vars
                value, value_pe, value_me = self.agnf_output[field]
                mu = self.props['magnification'].values[0]
                mu_err = self.props['magn_err'].values[0]

                if field == 'SFR_IR' or field == 'SFR_opt':
                    new_value = value / mu

                    new_pe = np.sqrt(np.power(value_pe / mu, 2.) + np.power((value / np.power(mu, 2.)) * mu_err, 2.))
                    new_me = np.sqrt(np.power(value_me / mu, 2.) + np.power((value / np.power(mu, 2.)) * mu_err, 2.))

                    return new_value, new_pe, new_me
                elif field == 'logMstar':
                    new_value = value - np.log10(mu)

                    new_pe = np.sqrt(np.power(value_pe, 2.) + np.power(mu_err / (np.log(10) * mu), 2.))
                    new_me = np.sqrt(np.power(value_me, 2.) + np.power(mu_err / (np.log(10) * mu), 2.))

                    return new_value, new_pe, new_me
                else:
                    pass
                    # print(f'No demagnification specified for {field}')

            elif demag:
                raise ValueError('Don\'t demag a different component than _sub')
            return self.agnf_output[field]
        raise AttributeError('Call load_agnf_output() first to load AGNfitter output.')
