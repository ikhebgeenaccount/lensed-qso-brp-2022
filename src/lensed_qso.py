import os

import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.legend_handler import HandlerTuple

from src.filters import get_wavelength

import warnings


FILTERED_SOURCES = {
    'B1152+200': ['panstarrs'],
    'B1600+434': ['panstarrs'],
    'B1608+656': [],
    'J0806+2006': ['panstarrs', 'Inada'],
    'J0924+0219': ['panstarrs'],
    'J1330+1810': ['panstarrs'],
    'J1455+1447': ['panstarrs'],
    'J1524+4409': ['panstarrs'],
    'J1633+3134': ['panstarrs'],
    'J1650+4251': ['panstarrs']
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
        self.sed.fillna(0, inplace=True)

        try:
            # Read mags
            self.mags = pd.read_csv(os.path.join('data', name, mags_source))
            self.mags.fillna(0, inplace=True)
        except FileNotFoundError:
            print('No mags file found for galaxy ' + self.name)

        self.props = pd.read_csv(os.path.join('data', properties), skiprows=1)

        # filtered_sed only selects entries that have a wavelength and a flux_total
        self.filtered_sed = self.sed[(self.sed.wavelength > 0) * (self.sed.flux_total > 0)].copy()

    def filter_sed(self, disallowed_sources=None, component=None):
        """
        Returns a DataFrame that contains only those rows that are allowed through LensedQSO.allowed_sources
        :param disallowed_sources: passed to LensedQSO.allowed_sources
        :return:
        """
        allowed_sources = self.allowed_sources(disallowed_sources)

        if component is None:
            return self.sed.loc[self.sed.source.isin(allowed_sources)]
        else:
            return self.sed.loc[self.sed.source.isin(allowed_sources) * (self.sed[f'flux{component}'] > 0)]

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
                warnings.warn(self.name + ': Filtering sources ' + str(disallowed_sources))

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
            sel = data[(data.source == l) * (data.wavelength > 0) * (data[data_type] > 0)]

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
            fig.savefig(os.path.join(self.save_location, 'SED_total.pdf'))
            fig.savefig(os.path.join(self.save_location, 'SED_total.png'))

        return fig, ax

    def save_sed(self):
        self.sed.to_csv(os.path.join('data', self.name, self.sed_file), index=False)

    def sed_to_agn_fitter(self):
        """
        Translates the SED to a catalog that is compatible with AGNfitter.
        :return: str
        """
        # TODO: use split data, model subtraction from total fluxes
        id = self.name.replace('B', '')
        id = id.replace('J', '')
        id = id.replace('+', '')

        header = '# ID redshift [wavelength_angstrom flux_mJy flux_error_mJy]\n'

        catalog = header + f'{id} {self.props.loc[self.props.galaxy == self.name].z_qso.values[0]} '
        for i, row in self.filtered_sed.iterrows():
            if not row.upper_limit:
                catalog += f'{row.wavelength} {row.flux_total} {row.flux_err} '
            else:
                # Upper limit has error -99 as flag for AGNfitter
                catalog += f'{row.wavelength} {row.flux_total} -99 '

        return catalog
