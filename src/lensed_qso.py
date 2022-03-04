import os

import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.legend_handler import HandlerTuple

from src.filters import get_wavelength


class LensedQSO:

    def __init__(self, name, sed_source='sed.csv',  mags_source='mags.csv', properties='properties.txt'):
        """
        :param name: Name of the lensed qso.
        :param sed_source: source of the spectral energy density (SED). Must be located in 'data/[name]'. Default is
        'data/[name]/sed.csv'.
        """
        self.name = name

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

    def plot_spectrum(self, loglog=False, mags=False, disallowed_sources=['panstarrs'], **kwargs):
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
            data_type = 'flux_total'
            data_err = 'flux_err'
            limit = 'upper_limit'

        # For every unique source, add their data separately
        u_sources = list(data.source.unique())
        to_remove = []

        if disallowed_sources is not None:
            # Check for occurrences of disallowed sources in the unique sources
            for ds in disallowed_sources:
                for s in u_sources:
                    if ds.lower() in s.lower():
                        to_remove.append(s)

            # Remove all found disallowed sources
            for r in to_remove:
                u_sources.remove(r)

        for l in u_sources:
            # Filter based on source
            # Only take those that have both a wavelength and a total flux
            sel = data[(data.source == l) * (data.wavelength > 0) * (data[data_type] > 0)]

            # If there are no entries after filtering, continue to next source
            if len(sel) == 0:
                continue

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
        ax.legend(legend_list, u_sources, loc='upper left', handler_map={tuple: HandlerTuple(ndivide=None)})
        ax.set_xscale('log')

        if loglog:
            ax.set_yscale('log')

        ax.set_title(f'{self.name} SED')
        ax.set_xlabel('$\mathit{Wavelength}\ (\mathrm{\AA})$')
        if mags:
            ax.set_ylabel('$\mathit{mag}$')
        else:
            ax.set_ylabel('$\mathit{Flux\ density}\ (\mathrm{mJy})$')

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
