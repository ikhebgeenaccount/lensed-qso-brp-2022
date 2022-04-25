import re
import os
import distutils.dir_util
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.legend_handler import HandlerTuple

from src.AGN_input import format_filter_name, format_telescope_name, get_agnf_filter_path
from src.filters import get_wavelength, FILTER_PROPERTIES

import warnings


FILTERED_SOURCES = {
    'B1152+200': ['panstarrs'],
    'B1600+434': ['panstarrs'],
    'B1608+656': [ 'luichies'],#['Koopmans+2003' ],
    'J0806+2006': ['panstarrs'],
    'J0924+0219': ['panstarrs', 'faure'],
    'J1330+1810': ['panstarrs'],
    'J1455+1447': ['panstarrs'],
    'J1524+4409': ['panstarrs'],
    'J1633+3134': ['panstarrs'],
    'J1650+4251': ['panstarrs']
}


FILTERED_SOURCES_AGNFITTER = {
    'B1152+200': ['toft'],
    'B1600+434': ['munoz'],
    'B1608+656': [ 'luichies'],#['Koopmans+2003' ],
    'J0806+2006': [],
    'J0924+0219': [],
    'J1330+1810': [],
    'J1455+1447': [],
    'J1524+4409': ['Oguri'],
    'J1633+3134': ['Morgan'],
    'J1650+4251': ['Morgan']
}


RADIO_CUTOFF = 1e8  # wavelengths >1e8 Angstrom are classified as radio
XRAY_CUTOFF = 300  # wavelengths < 300 Angstrom are classified as Xray

DEFAULT_AGNFITTER_SETTINGS = {
    'nwalkers': 100,
    'nburnsets': 2,
    'nburn': 5000,
    'nmcmc': 15000,
    'iprint': 1000,
    'plot_tracesburn-in': True,
    'plot_tracesmcmc': True,
    'plot_posteriortriangle': False,
    'realizations2int': 100,
    'realizations2plot': 10,
    'plotSEDrealizations': True
}

DEFAULT_AGNFITTER_SETTINSG_RX = {
    'plot_residuals': True,
    'modelset': 'modelsv1',
    'GALAXY': 'BC03_metal',
    'STARBURST': 'S17_radio',
    'BBB': 'SN12',
    'TORUS': 'SKIRTORM',
    'PRIOR_energy_balance': 'Flexible',
    'PRIOR_AGNfraction': True,
    'PRIOR_midIR_UV': False,
    'PRIOR_galaxy_only': False
}


COMPONENT_ID = {
    '_sub': 0,
    '_sub_demag': 1,
    '_sub_demag_test': 2
}


class LensedQSO:

    def __init__(self, name, sed_source='sed.csv',  mags_source='mags.csv', properties='properties.txt', save_all_plots=True, save_location='plots'):
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
        self.sed = pd.read_csv(os.path.join('data', name, sed_source))
        self.sed.fillna(0., inplace=True)

        try:
            # Read mags
            self.mags = pd.read_csv(os.path.join('data', name, mags_source))
            self.mags.fillna(0., inplace=True)
        except FileNotFoundError:
            print('No mags file found for galaxy ' + self.name)

        self.props = pd.read_csv(os.path.join('data', properties), skiprows=1)
        self.props = self.props.loc[self.props.galaxy == self.name]

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

        if self.name in FILTERED_SOURCES:
            # Add standard filtered sources to disallowed sources
            disallowed_sources += FILTERED_SOURCES[self.name]

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
            data = self.sed
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

            # Plot regular data points and upper limits separately, upper limits with special marker
            # TODO: consistent colours between all plots for same sources (SDSS, PanSTARRS, etc)
            le_1, _, _ = ax.errorbar(sel_reg.wavelength, sel_reg[data_type], sel_reg[data_err], fmt='o', label=l, **kwargs)

            if len(sel_upper_limit) > 0:
                le_2, _, _ = ax.errorbar(sel_upper_limit.wavelength, sel_upper_limit[data_type], sel_upper_limit[data_err],
                                         fmt='v', label=l, color=le_1.get_color(), **kwargs)
                legend_list.append((le_1, le_2))
            else:
                legend_list.append(le_1)

            # upper_limits += (le_2, )

        # ax.legend(legend_list + [upper_limits], list(self.sed.source.unique()) + ['upper limit'], loc='upper left', handler_map={tuple: HandlerTuple(ndivide=None)})
        ax.legend(legend_list, plotted_sources, loc='upper left', handler_map={tuple: HandlerTuple(ndivide=None)})
        ax.set_xscale('log')

        if loglog:
            ax.set_yscale('log')

        ax.set_title(f'{self.name} SED' if component is None else f'{self.name}{component} SED')
        ax.set_xlabel('$\mathit{Wavelength}\ (\mathrm{\AA})$')
        if mags:
            ax.set_ylabel('$\mathit{mag}$')
        else:
            ax.set_ylabel('$\mathit{Flux\ density}\ (\mathrm{mJy})$')

        if self.save_all_plots:
            fig.savefig(os.path.join(self.save_location, f'{self.name}_SED{component if component is not None else "_total"}.pdf'))
            # fig.savefig(os.path.join(self.save_location, f'{self.name}_SED{component if component is not None else "_total"}.png'))

        return fig, ax, plotted_sources, legend_list

    def plot_error_percentage(self, loglog=True):
        fig, ax = plt.subplots(figsize=(10, 8))
        ratio=self.sed['flux_sub_err']/self.sed['flux_sub'] *100
        ax.scatter(self.sed['wavelength'][ratio<=100],ratio[ratio<=100])
        ax.set_xscale('log')
        ax.set_ylabel('errors of sub in percentage', fontsize=15)
        ax.set_title(f'{self.name} percentage error')

        return fig,ax

    def save_sed(self):
        self.sed.to_csv(os.path.join('data', self.name, self.sed_file), index=False)

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
            for i, row in self.filter_sed(component=component, rX=rX, disallowed_sources=FILTERED_SOURCES_AGNFITTER[self.name]).iterrows():
                if row[f'flux{component}'] <= 0:
                    print('Skipping SED row', i)
                    continue

                if not row.upper_limit:
                    catalog_line += f'{row.wavelength} {row[f"flux{component}"]} {row[f"flux{component}_err"]} '
                else:
                    # Upper limit has error -99 as flag for AGNfitter
                    catalog_line += f'{row.wavelength} {row[f"flux{component}"]} -99 '

                if j== 0:
                    l += 3
            catalog += catalog_line + '\n'

        return catalog, l

    def agn_fitter_id(self, component='_sub'):
        t = self.name.replace('B', '').replace('J', '').replace('+', '') + str(COMPONENT_ID[component])
        return t if t[0] != '0' else t[1:]

    def agn_settings(self, rX=False, settings=None):
        # First load default settings to make sure every setting is set
        settings_use = DEFAULT_AGNFITTER_SETTINGS.copy()
        template = None

        if rX:
            settings_use.update(DEFAULT_AGNFITTER_SETTINSG_RX)

            template = settings_template_rX
            filters, filternames, filterfilenames = self.agn_fitter_input_filters(rX)

            hasxray = sum(self.filter_sed(component='_sub', rX=True)['wavelength'] <= XRAY_CUTOFF) > 0
            hasradio = sum(self.filter_sed(component='_sub', rX=True)['wavelength'] >= RADIO_CUTOFF) > 0
            print(hasxray, hasradio)
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
            template = settings_template
            settings_use.update({
                'length': self.sed_to_agn_fitter()[1],
                'length + 1': self.sed_to_agn_fitter()[1] + 1,
                'length + 2': self.sed_to_agn_fitter()[1] + 2,
                'name': self.name,
                'redshift': self.props.z_qso.values[0],
                'filters': self.agn_fitter_input_filters()
            })

        # Update default settings with settings that have been passed by user
        if settings is not None:
            settings_use.update(settings)

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

            exists_in_sed =self.filter_sed(component='_sub', rX=rX, disallowed_sources=FILTERED_SOURCES_AGNFITTER[self.name]).loc[(self.filter_sed(component='_sub', rX=rX, disallowed_sources=FILTERED_SOURCES_AGNFITTER[self.name])['telescope'] == tel) & (
                    self.filter_sed(component='_sub', rX=rX, disallowed_sources=FILTERED_SOURCES_AGNFITTER[self.name])['filter'] == fil)].shape[0] > 0

            if rX:
                filterfilenames += [get_agnf_filter_path(tel, fil)]
                filternames += [fil]

            outstr += f"    filters['_{format_telescope_name(tel)}_{format_filter_name(fil)}']={exists_in_sed}\n"

        if rX:
            return outstr, filternames, filterfilenames
        else:
            return outstr

    def load_agnf_output(self, rX=False, copy=False):
        self.agnf_output = [0] * len(COMPONENT_ID)
        for c in COMPONENT_ID.keys():
            self.agn_fitter_output(rX=rX, copy=copy, component=c, check_git=False)

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
        self.agnf_output = [0] * len(COMPONENT_ID)

        if verbose:
            print('Run\t-ll')

        for i in range(run_times):
            output = self.agn_fitter_output(rX=rX, agnf_id=self.agn_fitter_id(component=component) + str(i if i != 0 else ''), sub_folder=sub_folder)

            if output is None:
                continue

            lls.append(output['-ln_like'].values[2])
            index.append(i)

            if verbose:
                print(f'{i}\t{lls[-1]}')

        if verbose:
            print(f'Best run: {np.argmax(lls)}')

        return self.agn_fitter_output(rX=rX, agnf_id=self.agn_fitter_id(component=component) + str(index[np.argmax(lls)] if index[np.argmax(lls)] != 0 else ''), copy=copy, sub_folder=sub_folder)

    def agn_fitter_output(self, rX=False, agnf_id=None, copy=False, check_git=False, component='_sub', run_time=0, sub_folder=''):
        if agnf_id is None:
            agnf_id = self.agn_fitter_id(component=component) + (str(run_time) if run_time != 0 else '')
        if rX:
            raise NotImplementedError('rX path not yet')
        else:
            if sub_folder is None:
                path = os.path.join(os.pardir, 'AGNfitter', 'OUTPUT', str(agnf_id))
            else:
                path = os.path.join(os.pardir, 'AGNfitter', 'OUTPUT', sub_folder, str(agnf_id))
        repo_path = os.path.join('data', f'{self.name}', 'agnfitter')

        # Three cases:
        # path exists: points to AGNfitter output
        # repo_path exists: points to data in repo
        # neither exist: sucks
        if os.path.isdir(path):
            par_values_file = f'parameter_outvalues_{agnf_id}.txt'
            if copy:
                distutils.dir_util.copy_tree(path, os.path.join('data', self.name, 'agnfitter'))
        elif os.path.isdir(repo_path) and check_git:
            # print(re.escape(os.path.join(repo_path, f'parameter_outvalues_{self.agn_fitter_id(component=component)}')).replace('/', '\\/') + '[0-9]?\.txt')
            # par_values_file = glob.glob(re.escape(os.path.join(repo_path, f'parameter_outvalues_{self.agn_fitter_id(component=component)}')).replace('/', '\\/') + '[0-9]?\.txt')
            # print(par_values_file)
            par_values_file = f'parameter_outvalues_{agnf_id}.txt'
            path = repo_path
        else:
            # print('Are you working on vdesk or strw? If not, then this output has not been copied to our github repo, or doesn\'t exist')
            return

        cols = ['tau', 'age', 'Nh', 'irlum', 'SB', 'BB', 'GA', 'TO', 'EBVbbb', 'EBVgal', 'logMstar', 'SFR_opt', 'LIR(8-1000)', 'Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)', 'Lsb(1-30)', 'SFR_IR', '-ln_like']

        output = pd.read_csv(os.path.join(path, par_values_file),
                             delim_whitespace=True, skiprows=4, header=None, names=cols)

        cid = COMPONENT_ID[component]
        self.agnf_output[cid] = {}
        for c in cols:
            self.agnf_output[cid][c] = []

            self.agnf_output[cid][c].append(output[c].iloc[2])
            self.agnf_output[cid][c].append(output[c].iloc[3] - self.agnf_output[cid][c][0])
            self.agnf_output[cid][c].append(self.agnf_output[cid][c][0] - output[c].iloc[1])

        return output

    def get_agnf_output_field(self, field, component='_sub', demag=False):
        if hasattr(self, 'agnf_output'):

            if demag and component == '_sub':
                # Save values in easy vars
                value, value_pe, value_me = self.agnf_output[COMPONENT_ID[component]][field]
                mu = self.props['magnification'].values[0]
                mu_err = self.props['magn_err'].values[0]

                if field == 'SFR_IR' or field == 'SFR_opt':
                    new_value = value / mu

                    new_pe = np.sqrt(np.power(value_pe / mu, 2.) + np.power(value / np.power(mu, 2.) * mu_err, 2.))
                    new_me = np.sqrt(np.power(value_me / mu, 2.) + np.power(value / np.power(mu, 2.) * mu_err, 2.))

                    return new_value, new_pe, new_me
                elif field == 'logMstar':
                    new_value = value - np.log10(mu)

                    new_pe = np.sqrt(np.power(value_pe, 2.) + np.power(mu_err / (np.log(10) * mu), 2.))
                    new_me = np.sqrt(np.power(value_me, 2.) + np.power(mu_err / (np.log(10) * mu), 2.))

                    return new_value, new_pe, new_me
                else:
                    print(f'No demagnification specified for {field}')

            elif demag:
                raise ValueError('Don\'t demag a different component than _sub')
            return self.agnf_output[COMPONENT_ID[component]][field]
        raise AttributeError('Call load_agnf_output() first to load AGNfitter output.')




settings_template_rX = "'''\n" +\
"AGNfitter setting file:\n" +\
"required:\n" +\
"CATALOG_settings\n" +\
"FILTERS_settings\n" +\
"MCMC_settings\n" +\
"OUTPUT_settings\n" +\
"For default use (test example with 2 redshifts and default filter set)\n" +\
"Change only the functions which state \n" +\
"***USER INPUT NEEDED***.\n" +\
"'''\n" +\
"import os\n" +\
"\n" +\
"def CATALOG_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    ***USER INPUT NEEDED***\n" +\
"    Set the right values to be able to read your catalog's format.\n" +\
"    FITS option is not available yet.\n" +\
"    =================================='''\n" +\
"\n" +\
"\n" +\
"    cat = dict()\n" +\
"\n" +\
"\n" +\
"    ##GENERAL\n" +\
"    cat['path'] = os.getcwd().replace('lensed-qso-brp-2022', 'AGNfitter_rXv0.1/AGNfitter/') #'/Users/gcalistr/Documents/AGNfitter/'  #path to the AGNfitter code\n" +\
"    cat['filename'] = cat['path']+'data/{name}.txt'\n" +\
"    cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. \n" +\
"    cat['name'] = 0                 ## If ASCII: Column index (int) of source IDs\n" +\
"                                    ## If FITS : Column name (str). E.g. 'ID'\n" +\
"    cat['redshift'] = 1             ## If ASCII:  Column index(int) of redshift \n" +\
"                                     ## If FITS : Column name (str). E.g. z'\n" +\
"\n" +\
"   ##FREQUENCIES/WAVELENGTHS \n" +\
"    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'\n" +\
"    cat['freq/wl_list'] = np.arange(2,{length},3).tolist()  #138 147\n" +\
"                                        ## If ASCII: List of column indexes (int), \n" +\
"                                        ##           corresponding to freq/wl.                                  \n" +\
"    #cat['freq/wl_suffix'] = '_wl'      ## If FITS: common ending to wavelength column names\n" +\
"\n" +\
"    cat['use_central_wavelength'] = True # Option to use central wavelength if no wavelengths in table\n" +\
"\n" +\
"    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*\n" +\
"                                         ## 'frequency' or 'wavelength'?\n" +\
"    cat['freq/wl_unit'] = u.Angstrom #u. Angstrom u.micron     ## Astropy unit of freq or wavelength\n" +\
"\n" +\
"    ##FLUXES \n" +\
"    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'\n" +\
"    cat['flux_in_magAB'] = False # Option to calculate flux and flux_error from magnitude AB.\n" +\
"    cat['flux_unit'] = u.Jy * 1e-3            ## Astropy unit of *flux* (astropy-units)\n" +\
"    cat['flux_list'] = np.arange(3,{length + 1},3).tolist()      #139 148\n" +\
"                                        ## If ASCII: List of column indexes (int)\n" +\
"    #cat['flux_suffix'] = '_f'          ## If FITS: Common ending of all flux column names (str)    \n" +\
"    cat['fluxerr_list'] = np.arange(4,{length + 2},3).tolist() #140 149\n" +\
"                                        ## If ASCII: List of column indexes (int)\n" +\
"    #cat['fluxerr_suffix'] = '_e'       ## If FITS: common ending to fluxerr column names (str)\n" +\
"    cat['err+10%flux_moreflex'] = False\n" +\
"    ##NON-DETECTIONS                                        \n" +\
"    cat['ndflag_bool'] = False          ## Does you catalog has columns with flags 1(0) for \n" +\
"                                        ## detections (nondetections)? \n" +\
"    cat['ndflag_list'] = 'list'         ## If ASCII: List of column indexes (int)\n" +\
"                                        ## If FITS: List of column names (str)    \n" +\
"\n" +\
"    ## COSTUMIZED WORKING PATHS\n" +\
"    cat['workingpath'] = cat['path']  # Allows for a working path other than the AGNfitter code path.\n" +\
"                                      # Will include:\n" +\
"                                            # dictionary of models \n" +\
"                                            # SETTINGS_AGNFitter.py file  \n" +\
"                                            # OUTPUT\n" +\
"                                      # Specially needed in order not to alter git original repository\n" +\
"                                      # and when using an external processor.\n" +\
"                                      # Default: cat['path'] (same as AGNfitter code path) \n" +\
"                                      \n" +\
"    cat['output_folder'] =  cat['workingpath'] +'OUTPUT/' #if no special OUTPUT folder, leave default\n" +\
"\n" +\
"\n" +\
"\n" +\
"    return cat\n" +\
"\n" +\
"\n" +\
"def FILTERS_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    Set the photometric bands included in your catalog,\n" +\
"    in order to integrate the models over their response curves.\n" +\
"    =================================='''\n" +\
"\n" +\
"    filters = dict()\n" +\
"\n" +\
"    filters['dict_zarray'] = np.array([{redshift}]) # [0.058900, 0.042170] np.array([0.283, 1.58])  # Deprecated. The grid of redshifts needed to fit your catalog\n" +\
"    filters['path'] = 'models/FILTERS/' \n" +\
"    filters['filterset'] = 'example_46datapointa' ## 'filterset_default' (for the test case),\n" +\
"                                               ## for the user's case: customize, eg. filtersv1\n" +\
"\n" +\
"{filters}\n" +\
"\n" +\
"\n" +\
"    filters['add_filters']= True # If 'True' please add them below in ADD FILTERS\n" +\
"\n" +\
"    '''==================================\n" +\
"    ADD FILTERS (optional)\n" +\
"    =================================='''\n" +\
"\n" +\
"    ADDfilters=dict()\n" +\
"    ADDfilters['names'] = {filternames}  ## (string/list of strings)User especified filter names.\n" +\
"                                ## If name has been previously used, an error message will appear.\n" +\
"    ADDfilters['filenames'] = {filterfilenames}  ## (string/list of strings) File names of the new filters.\n" +\
"                                ## File format: 2 columns of 1) freq/wavelength 2) Throughput. \n" +\
"                                ## Path assumed is the cat['path'] especified above. \n" +\
"                                ## Example: 'models/FILTERS/my_new_filter.txt'\n" +\
"    ADDfilters['freq/wl_format'] = ['wavelength'] * len(ADDfilters['names']) ## Info about the column 1 of your filter file.\n" +\
"                                                                             ## Options: 'wavelength' or 'frequency'.    \n" +\
"    ADDfilters['freq/wl_unit'] =  [u.Angstrom]* len(ADDfilters['names']) ## (Astropy Unit) Info about the column 1 \n" +\
"                                                                         ## of your filter file. \n" +\
"    ADDfilters['description'] = ['description_dummy']* len(ADDfilters['names']) ## (Str) Any description the user wants to give \n" +\
"                                                                                ##  to the filter to add.\n" +\
"\n" +\
"    filters['add_filters_dict']= ADDfilters\n" +\
"\n" +\
"    return filters\n" +\
"\n" +\
"def MODELS_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    Work in progress\n" +\
"    =================================='''\n" +\
"\n" +\
"\n" +\
"    models = dict()\n" +\
"    models['path'] = 'models/' \n" +\
"    models['modelset'] = {modelset}\n" +\
"\n" +\
"\n" +\
"    models['GALAXY'] = {GALAXY}   ### Current options:\n" +\
"                                ### 'BC03' (Bruzual & Charlot 2003)\n" +\
"                                ### 'BC03_metal' (Bruzual & Charlot 2003), with metallicities\n" +\
"    models['STARBURST'] = {STARBURST} ### Current options:\n" +\
"                                ### 'DH02_CE01' (Dale & Helou 2002 + Chary & Elbaz 2001)\n" +\
"                                ### 'S17' (Schreiber et al. 2017 (submitted))  \n" +\
"                                ### 'S17_newmodel'\n" +\
"                                ### 'S17_radio' \n" +\
"\n" +\
"    models['BBB'] = {BBB} ### Current options:\n" +\
"                         ### 'R06' (Richards et al. 2006) ## Needs 2 manual changes in PARAMETERSPACE_AGNfitter.py\n" +\
"                         ### 'SN12' (Slone&Netzer 2012)\n" +\
"                         ### 'D12_S' (Done et al. 2012) for Schwarzschild BH, with x-ray predictions\n" +\
"                         ### 'D12_K' (Done et al. 2012) for Kerr BH, with x-ray predictions\n" +\
"\n" +\
"    models['TORUS'] ={TORUS} ### Current options:\n" +\
"                           ### 'S04' (Silva et al. 2004)\n" +\
"                           ### 'NK0' (Nenkova et al. 2008)\n" +\
"                           ### 'NK0_2P' (Nenkova et al. 2008) with averaged SEDs for each inclination and openning angle\n" +\
"                           ### 'NK0_3P' (Nenkova et al. 2008) with averaged SEDs for each inclination, openning angle and optical depth\n" +\
"                           ### 'NK0_4P' (Nenkova et al. 2008)\n" +\
"                           ### 'SKIRTOR' (Stalevski et al. 2016)\n" +\
"                           ### 'SKIRTORC' with parameter values used in X-CIGALE\n" +\
"                           ### 'SKIRTORM' SKIRTOR model with averaged SEDs for each inclination\n" +\
"                           ### 'SKIRTORM_2P' SKIRTOR model with averaged SEDs for each inclination and openning angle\n" +\
"\n" +\
"    models['XRAYS'] = {hasxray} ### If X-ray data is available and informative for the fit\n" +\
"\n" +\
"    models['RADIO'] = {hasradio} ### If radio data is available and informative for the fit\n" +\
"\n" +\
"    models['PRIOR_energy_balance'] = {PRIOR_energy_balance} ### Default:'Flexible'\n" +\
"                                          ### 'Flexible': Sets a lower limit to the dust emission luminosity ('starburst' model)\n" +\
"                                          ### as given by the observed attenuation in the stellar component SED.\n" +\
"                                          ### 'Restrictive': Favours a combination of parameters such that the luminosity of the cold dust \n" +\
"                                          ### and that attenuated in the stellar component are equal.\n" +\
"    models['PRIOR_AGNfraction'] = {PRIOR_AGNfraction}  ### Default: True\n" +\
"                                        ### True: - *IF* blue/UV bands (around 1500 Angstrom) are 10 times higher than expected by the galaxy luminosity function by Parsa, Dunlop et al. 2014. \n" +\
"                                        ###         this option rejects AGN-to-GAL ratios lower than 1 (log =0). It then applies a Gaussian prior probability with log ratio=2, with a sigma of 2.\n" +\
"                                        ###       - In this cases it also applies a Gaussian prior on the galaxy normalization, i.e. stellar mass (usually unconstrained in these cases) to \n" +\
"                                        ###         populate physically expected ranges for QSO hosts -> 10^9 - 10^11. \n" +\
"                                        ###       - *ELSE IF* blue/UV bands (around 1500 Angstrom) are below 10 times the expected value by Parsa, Dunlop et al. 2014. \n" +\
"                                        ###         this option gives preference to galaxy contribution in the optical UV, with Gaussian prior probability centered on AGN to GALAXY log ratios of -1. \n" +\
"                                        ###          and sigma 1, i.e. accretion disk is disfavoured at least the data strongly prefers it.\n" +\
"                                        ### False:- Non-informative prior\n" +\
"    models['PRIOR_midIR_UV'] = {PRIOR_midIR_UV}\n" +\
"    models['PRIOR_galaxy_only'] = {PRIOR_galaxy_only} ### Default:False \n" +\
"                                        ### True: sets all AGN contribution to 0.ß\n" +\
"    return models\n" +\
"\n" +\
"def MCMC_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    Set your preferences for the MCMC sampling.\n" +\
"    =================================='''\n" +\
"\n" +\
"    mc = dict()\n" +\
"\n" +\
"    mc['Nwalkers'] = {nwalkers}  ## 100 number of walkers #100\n" +\
"    mc['Nburnsets']= {nburnsets}   ## number of burn-in sets\n" +\
"    mc['Nburn'] = {nburn} ## length of each burn-in sets\n" +\
"    mc['Nmcmc'] = {nmcmc}  ## length of each burn-in sets\n" +\
"    mc['iprint'] = {iprint} ## show progress in terminal in steps of this many samples\n" +\
"\n" +\
"    return mc\n" +\
"\n" +\
"def OUTPUT_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    Set your preferences for the production of OUTPUT files. \n" +\
"    =================================='''\n" +\
"\n" +\
"    out = dict()\n" +\
"\n" +\
"    out['plot_format'] = 'pdf'\n" +\
"\n" +\
"    #CHAIN TRACES\n" +\
"    out['plot_tracesburn-in'] = {plot_tracesburn-in}    \n" +\
"    out['plot_tracesmcmc'] = {plot_tracesmcmc}\n" +\
"\n" +\
"    #BASIC OUTPUT\n" +\
"    out['Nsample'] = 1000 ## out['Nsample'] * out['Nthinning'] <= out['Nmcmc']\n" +\
"    out['Nthinning'] = 10 ## This describes thinning of the chain to sample\n" +\
"    out['writepar_meanwitherrors'] = True ##Write output values for all parameters in a file.\n" +\
"    out['plot_posteriortriangle'] = {plot_posteriortriangle} ##Plot triangle with all parameters' PDFs?\n" +\
"\n" +\
"    #INTEGRATED LUMINOSITIES\n" +\
"    out['calc_intlum'] = True  \n" +\
"    out['save_posterior_luminosities']= False\n" +\
"    out['save_posteriors']= True\n" +\
"    out['realizations2int'] = {realizations2int} #This process is very time consuming.\n" +\
"                                #Around 100-1000 is recomendend for computational reasons.\n" +\
"                                #If you want to plot posterior triangles of \n" +\
"                                #the integrated luminosities, should be > 1000.\n" +\
"    out['plot_posteriortrianglewithluminosities'] = False  # requires out['calc_intlum']=True \n" +\
"\n" +\
"    #INTEGRATION RANGES\n" +\
"    out['intlum_models'] = ['sb','bbb', 'bbbdered', 'gal', 'tor','agn_rad', 'tor+bbb','AGNfrac_IR', 'gal+bbb', 'AGNfrac_opt', 'bbb', 'tor','sb']  #leave 'sb' always as first element\n" +\
"    out['intlum_freqranges_unit'] = u.micron   #Astropy unit \n" +\
"    out['intlum_freqranges'] = np.array([[8.,1000.],[0.1,1.],[0.1,1.],[0.1,1.],[1.,30.],[0.1, 15000.], [0.1, 30], [8.,1000.],[0.4, 0.5], [0.4, 0.5], [1.2398e-4, 6.1992e-4], [6,6],   [1.,30.]])\n" +\
"    out['intlum_names'] = ['LIR(8-1000)','Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)','Lagn_rad(0.1-15000)', 'LAGN(0.1-30)', 'AGNfrac(8-1000)', 'Lgal_bbb(0.4-0.5)', 'AGNfrac(0.4-0.5)', 'Lxr(2-10keV)',  'Ltor(6)', 'Lsb(1-30)']\n" +\
"\n" +\
"    #SED PLOTTING\n" +\
"    out['realizations2plot'] = {realizations2plot}\n" +\
"    out['plot_residuals']= {plot_residuals}\n" +\
"    out['saveSEDresiduals'] = True\n" +\
"    out['plotSEDrealizations'] = {plotSEDrealizations}\n" +\
"    out['saveSEDrealizations'] = True\n" +\
"\n" +\
"    return out"

settings_template = "'''\n" +\
"AGNfitter setting file:\n" +\
"\n" +\
"required:\n" +\
"  CATALOG_settings\n" +\
"  FILTERS_settings\n" +\
"  MCMC_settings\n" +\
"  OUTPUT_settings\n" +\
"\n" +\
"For default use (test example with 2 redshifts and default filter set)\n" +\
"\n" +\
"Change only the functions which state\n" +\
"***USER INPUT NEEDED***.\n" +\
"'''\n" +\
"import os\n" +\
"def CATALOG_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    ***USER INPUT NEEDED***\n" +\
"\n" +\
"    Set the right values to be able to read your catalog's format.\n" +\
"    FITS option is not available yet.\n" +\
"    =================================='''\n" +\
"\n" +\
"\n" +\
"    cat = dict()\n" +\
"\n" +\
"\n" +\
"    ##GENERAL\n" +\
"    cat['path'] = os.getcwd().replace('lensed-qso-brp-2022', 'AGNfitter/')  #path to the AGNfitter code\n" +\
"\n" +\
"\n" +\
"    cat['filename'] = cat['path']+'data/{name}.txt'\n" +\
"    cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'.\n" +\
"                              ## FITS option not available yet.\n" +\
"\n" +\
"    cat['name'] = 0#'ID'            ## If ASCII: Column index (int) of source IDs\n" +\
"                                        ## If FITS: not yet\n" +\
"    cat['redshift'] = 1#'z'              ## If ASCII:  Column index(int) of redshift\n" +\
"                                        ## If FITS: not yet\n" +\
"\n" +\
"    ##FREQUENCIES/WAVELENGTHS\n" +\
"    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'\n" +\
"    cat['freq/wl_list'] = np.arange(2,{length},3).tolist()\n" +\
"                                        ## If ASCII: List of column indexes (int),\n" +\
"                                        ##           corresponding to freq/wl.\n" +\
"    #cat['freq/wl_suffix'] = '_wl'      ## If FITS: common ending to wavelength column names\n" +\
"    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*\n" +\
"                                         ## 'frequency' or 'wavelength'?\n" +\
"    cat['freq/wl_unit'] = u.Angstrom       ## Astropy unit of freq or wavelength\n" +\
"\n" +\
"    ##FLUXES\n" +\
"    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'\n" +\
"    cat['flux_unit'] = u.mJy             ## Astropy unit of *flux* (astropy-units)\n" +\
"    cat['flux_list'] = np.arange(3,{length + 1},3).tolist()\n" +\
"                                        ## If ASCII: List of column indexes (int)\n" +\
"    #cat['flux_suffix'] = '_f'          ## If FITS: Common ending of all flux column names (str)\n" +\
"\n" +\
"    cat['fluxerr_list'] = np.arange(4,{length + 2},3).tolist()\n" +\
"                                        ## If ASCII: List of column indexes (int)\n" +\
"    #cat['fluxerr_suffix'] = '_e'       ## If FITS: common ending to fluxerr column names (str)\n" +\
"\n" +\
"    ##NON-DETECTIONS\n" +\
"    cat['ndflag_bool'] = False          ## Does you catalog has columns with flags 1(0) for\n" +\
"                                        ## detections (nondetections)?\n" +\
"    cat['ndflag_list'] = 'list'         ## If ASCII: List of column indexes (int)\n" +\
"                                        ## If FITS: List of column names (str)\n" +\
"\n" +\
"    ## COSTUMIZED WORKING PATHS\n" +\
"    cat['workingpath'] = cat['path']  # Allows for a working path other than the AGNfitter code path.\n" +\
"                                      # Will include:\n" +\
"                                            # dictionary of models\n" +\
"                                            # SETTINGS_AGNFitter.py file\n" +\
"                                            # OUTPUT\n" +\
"                                      # Specially needed in order not to alter git original repository\n" +\
"                                      # and when using an external processor.\n" +\
"                                      # Default: cat['path'] (same as AGNfitter code path)\n" +\
"\n" +\
"    cat['output_folder'] =  cat['workingpath'] +'OUTPUT/' #if no special OUTPUT folder, leave default\n" +\
"    cat['dict_path'] = cat['workingpath'] + 'models/MODELSDICT_default'\n" +\
"\n" +\
"\n" +\
"    return cat\n" +\
"\n" +\
"\n" +\
"def FILTERS_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    Set the photometric bands included in your catalog,\n" +\
"    in order to integrate the models over their response curves.\n" +\
"    =================================='''\n" +\
"\n" +\
"    filters = dict()\n" +\
"\n" +\
"    filters['dict_zarray'] =np.array([{redshift}])  # The grid of redshifts needed to fit your catalog\n" +\
"    filters['Bandset'] = 'BANDSET_settings' # OPTIONS:\n" +\
"                                           # 'BANDSET_default' (for testing)\n" +\
"                                           # 'BANDSET_settings' (choosing relevant filters below, as given by your catalog)\n" +\
"                                           # if your filter is not included, go to DICTIONARIES_AGNfitter to add.\n" +\
"\n" +\
"{filters}\n" +\
"    return filters\n" +\
"\n" +\
"def MCMC_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    Set your preferences for the MCMC sampling.\n" +\
"    =================================='''\n" +\
"\n" +\
"    mc = dict()\n" +\
"\n" +\
"    mc['Nwalkers'] = {nwalkers}  ## number of walkers\n" +\
"    mc['Nburnsets']= {nburnsets}   ## number of burn-in sets\n" +\
"    mc['Nburn'] = {nburn} ## length of each burn-in sets\n" +\
"    mc['Nmcmc'] = {nmcmc}  ## length of each burn-in sets\n" +\
"    mc['iprint'] = {iprint} ## show progress in terminal in steps of this many samples\n" +\
"\n" +\
"    return mc\n" +\
"\n" +\
"def OUTPUT_settings():\n" +\
"\n" +\
"    '''==================================\n" +\
"    Set your preferences for the production of OUTPUT files.\n" +\
"    =================================='''\n" +\
"\n" +\
"    out = dict()\n" +\
"\n" +\
"    out['plot_format'] = 'pdf'\n" +\
"\n" +\
"    #CHAIN TRACES\n" +\
"    out['plot_tracesburn-in'] = {plot_tracesburn-in}\n" +\
"    out['plot_tracesmcmc'] = {plot_tracesmcmc}\n" +\
"\n" +\
"    #BASIC OUTPUT\n" +\
"    out['Nsample'] = 1000 ## out['Nsample'] * out['Nthinning'] <= out['Nmcmc']\n" +\
"    out['Nthinning'] = 10 ## This describes thinning of the chain to sample\n" +\
"    out['writepar_meanwitherrors'] = True ##Write output values for all parameters in a file.\n" +\
"    out['plot_posteriortriangle'] = {plot_posteriortriangle} ##Plot triangle with all parameters' PDFs?\n" +\
"\n" +\
"    #INTEGRATED LUMINOSITIES\n" +\
"    out['calc_intlum'] = True\n" +\
"    out['realizations2int'] = {realizations2int} #This process is very time consuming.\n" +\
"                                #Around 100-1000 is recomendend for computational reasons.\n" +\
"                                #If you want to plot posterior triangles of\n" +\
"                                #the integrated luminosities, should be > 1000.\n" +\
"    out['plot_posteriortrianglewithluminosities'] = False  # requires out['calc_intlum']=True\n" +\
"\n" +\
"    #INTEGRATION RANGES\n" +\
"    out['intlum_models'] = ['sb','bbb', 'bbbdered', 'gal', 'tor','sb']  #leave 'sb' always\n" +\
"                                                                        #as first element\n" +\
"    out['intlum_freqranges_unit'] = u.micron   #Astropy unit\n" +\
"    out['intlum_freqranges'] = np.array([[8.,1000.],[0.1,1.],[0.1,1.],[0.1,1.],[1.,30.],[1.,30.]])\n" +\
"    out['intlum_names'] = ['LIR(8-1000)','Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)','Lsb(1-30)']\n" +\
"\n" +\
"    #SED PLOTTING\n" +\
"    out['realizations2plot'] = {realizations2plot}\n" +\
"\n" +\
"    out['plotSEDrealizations'] = {plotSEDrealizations}\n" +\
"\n" +\
"    return out\n" +\
""
