import os
import distutils.dir_util
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.legend_handler import HandlerTuple

from src.filters import get_wavelength, FILTER_PROPERTIES

import warnings


FILTERED_SOURCES = {
    'B1152+200': ['panstarrs', 'chandra', 'luichies'],
    'B1600+434': ['panstarrs', 'chandra', 'luichies'],
    'B1608+656': ['chandra', 'luichies'],#['Koopmans+2003' ],
    'J0806+2006': ['panstarrs', 'chandra', 'luichies'],
    'J0924+0219': ['panstarrs', 'faure', 'castles', 'chandra', 'luichies'],
    'J1330+1810': ['panstarrs', 'chandra', 'luichies'],
    'J1455+1447': ['panstarrs', 'chandra', 'luichies'],
    'J1524+4409': ['panstarrs', 'chandra', 'luichies'],
    'J1633+3134': ['panstarrs', 'chandra', 'luichies'],
    'J1650+4251': ['panstarrs', 'chandra', 'luichies']
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
        self.save_location = os.path.join(save_location, name)

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

    def filter_sed(self, disallowed_sources=None, component='_total', allow_zero_error=True):
        """
        Returns a DataFrame that contains only those rows that are allowed through LensedQSO.allowed_sources
        :param disallowed_sources: passed to LensedQSO.allowed_sources
        :return:
        """
        allowed_sources = self.allowed_sources(disallowed_sources)

        compfilter = np.ones(self.sed.shape[0], dtype=bool)
        if component is not None:
            compfilter = self.sed[f'flux{component}'] > 0

        errorfilter = np.ones(self.sed.shape[0], dtype=bool)
        if not allow_zero_error:
            errorfilter= self.sed[f'flux{"" if component == "_total" else component}_err'] > 0

        return self.sed.loc[self.sed.source.isin(allowed_sources) * compfilter * errorfilter]

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
            fig.savefig(os.path.join(self.save_location, f'SED{component if component is not None else "_total"}.pdf'))
            fig.savefig(os.path.join(self.save_location, f'SED{component if component is not None else "_total"}.png'))

        return fig, ax

    def save_sed(self):
        self.sed.to_csv(os.path.join('data', self.name, self.sed_file), index=False)

    def sed_to_agn_fitter(self):
        """
        Translates the SED to a catalog that is compatible with AGNfitter.
        :return: str
        """
        # TODO: use split data, model subtraction from total fluxes
        id = self.agn_fitter_id()

        header = '# ID redshift [wavelength_angstrom flux_mJy flux_error_mJy]\n'

        l = 2

        catalog = header + f'{id} {self.props.z_qso.values[0]} '
        for i, row in self.filter_sed(component='_sub').iterrows():
            if row['flux_sub'] <= 0:
                print('Skipping SED row', i)
                continue

            if not row.upper_limit:
                catalog += f'{row.wavelength} {row.flux_sub} {row.flux_sub_err} '
            else:
                # Upper limit has error -99 as flag for AGNfitter
                catalog += f'{row.wavelength} {row.flux_sub} -99 '

            l += 3

        return catalog, l

    def agn_fitter_id(self):
        t = self.name.replace('B', '').replace('J', '').replace('+', '')
        return t if t[0] != '0' else t[1:]

    def agn_settings(self):
        return settings_template.format(**{
            'length': self.sed_to_agn_fitter()[1],
            'length + 1': self.sed_to_agn_fitter()[1] + 1,
            'length + 2': self.sed_to_agn_fitter()[1] + 2,
            'name': self.name,
            'redshift': self.props.z_qso.values[0],
            'filters': self.agn_fitter_input_filters()
        })

    def agn_fitter_input_filters(self):
        outstr = ''

        # Iterate over each row to get the dictionary settings
        for i, row in FILTER_PROPERTIES.iterrows():

            tel = FILTER_PROPERTIES.telescope.values[i]
            fil = FILTER_PROPERTIES.filtername.values[i]

            # if [tel,fil] in lqso.sed[['telescope', 'filter']].values:
            if self.filter_sed(component='_sub').loc[(self.filter_sed(component='_sub')['telescope'] == tel) & (
                    self.filter_sed(component='_sub')['filter'] == fil)].shape[0] > 0:
                f = fil.replace("'", 'prime').replace('-', 'dash').replace('/', 'slash').replace('.', 'dot')
                outstr += f"    filters['_{tel.replace(' ', 'space')}_{f}']=True\n"
            else:
                f = fil.replace("'", 'prime').replace('-', 'dash').replace('/', 'slash').replace('.', 'dot')
                outstr += f"    filters['_{tel.replace(' ', 'space')}_{f}']=False\n"

        return outstr

    def agn_fitter_output(self, rX=False, agnf_id=None, copy=False):
        if agnf_id is None:
            agnf_id = self.agn_fitter_id()
        if rX:
            print('rX path not yet')
        else:
            path = os.path.join(os.pardir, 'AGNfitter', 'OUTPUT', str(agnf_id))

        if not os.path.isdir(path):
            print('No results found.')
            print('This function only works when working from strw/vdesk.')
            return

        if copy:
            distutils.dir_util.copy_tree(path, os.path.join('data', self.name, 'agnfitter'))

        # TODO: some column names contain spaces? Like log Mstar
        result = pd.read_csv(os.path.join(path, f'parameter_outvalues_{agnf_id}.txt'),
                             delim_whitespace=True, skiprows=4, header=None, names=['tau', 'age', 'Nh', 'irlum', 'SB', 'BB', 'GA', 'TO', 'EBVbbb', 'EBVgal', 'logMstar', 'SFR_opt', 'LIR(8-1000)', 'Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)', 'Lsb(1-30)', 'SFR_IR', '-ln_like'])

        return result


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
"\n" +\
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
"    cat['path'] ='/data2/bach1/abbo/brp/AGNfitter/'  #path to the AGNfitter code\n" +\
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
"    mc['Nwalkers'] = 100  ## number of walkers\n" +\
"    mc['Nburnsets']= 2   ## number of burn-in sets\n" +\
"    mc['Nburn'] = 4000 ## length of each burn-in sets\n" +\
"    mc['Nmcmc'] = 10000  ## length of each burn-in sets\n" +\
"    mc['iprint'] = 1000 ## show progress in terminal in steps of this many samples\n" +\
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
"    out['plot_tracesburn-in'] = False\n" +\
"    out['plot_tracesmcmc'] = True\n" +\
"\n" +\
"    #BASIC OUTPUT\n" +\
"    out['Nsample'] = 1000 ## out['Nsample'] * out['Nthinning'] <= out['Nmcmc']\n" +\
"    out['Nthinning'] = 10 ## This describes thinning of the chain to sample\n" +\
"    out['writepar_meanwitherrors'] = True ##Write output values for all parameters in a file.\n" +\
"    out['plot_posteriortriangle'] = False ##Plot triangle with all parameters' PDFs?\n" +\
"\n" +\
"    #INTEGRATED LUMINOSITIES\n" +\
"    out['calc_intlum'] = True\n" +\
"    out['realizations2int'] = 100 #This process is very time consuming.\n" +\
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
"    out['realizations2plot'] = 10\n" +\
"\n" +\
"    out['plotSEDrealizations'] = True\n" +\
"\n" +\
"    return out\n" +\
""
