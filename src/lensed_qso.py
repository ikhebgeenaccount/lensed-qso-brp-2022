import os

import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.legend_handler import HandlerTuple


class LensedQSO:

    def __init__(self, name, sed_source='sed.csv'):
        """
        :param name: Name of the lensed qso.
        :param sed_source: source of the spectral energy density (SED). Must be located in 'data/[name]'. Default is
        'data/[name]/sed.csv'.
        """
        self.name = name
        self.sed = pd.read_csv(os.path.join('data', name, sed_source))

    def plot_spectrum(self, loglog=False, **kwargs):
        fig, ax = plt.subplots(figsize=(10, 8))

        legend_list = []

        # upper_limits = ()

        # For every unique source, add their data separately
        for l in self.sed.source.unique():
            # Filter based on source
            # Only take those that have both a wavelength and a total flux
            sel = self.sed[(self.sed.source == l) * (self.sed.wavelength > 0) * (self.sed.flux_total > 0)]

            # If there are no entries after filtering, continue to next source
            if len(sel) == 0:
                continue

            # Separate upper limits from regular data points
            sel_upper_limit = sel[sel.upper_limit == 1]
            sel_reg = sel[sel.upper_limit == 0]

            # Plot regular data points and upper limits separately, upper limits with special marker
            le_1, _, _ = ax.errorbar(sel_reg.wavelength, sel_reg.flux_total, sel_reg.flux_err, fmt='o', label=l, **kwargs)

            if len(sel_upper_limit) > 0:
                le_2, _, _ = ax.errorbar(sel_upper_limit.wavelength, sel_upper_limit.flux_total, sel_upper_limit.flux_err,
                                         fmt='o', label=l, marker='v', color=le_1.get_color(), **kwargs)
                legend_list.append((le_1, le_2))
            else:
                legend_list.append(le_1)

            # upper_limits += (le_2, )

        # ax.legend(legend_list + [upper_limits], list(self.sed.source.unique()) + ['upper limit'], loc='upper left', handler_map={tuple: HandlerTuple(ndivide=None)})
        ax.legend(legend_list, self.sed[(self.sed.wavelength > 0) * (self.sed.flux_total > 0)].source.unique(), loc='upper left', handler_map={tuple: HandlerTuple(ndivide=None)})
        ax.set_xscale('log')

        if loglog:
            ax.set_yscale('log')

        ax.set_title(f'{self.name} SED')
        ax.set_xlabel('$\mathit{Wavelength}\ (\mathrm{\AA})$')
        ax.set_ylabel('$\mathit{Flux\ density}\ (\mathrm{mJy})$')

        return fig, ax
