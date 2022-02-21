import os

import matplotlib.pyplot as plt
import pandas as pd


class LensedQSO:

	def __init__(self, name, sed_source='sed.csv'):
		"""
		:param name: Name of the lensed qso.
		:param sed_source: source of the spectral energy density (SED). Must be located in 'data/[name]'. Default is
		'data/[name]/sed.csv'.
		"""
		self.name = name
		self.sed = pd.read_csv(os.path.join('data', name, sed_source))

	def plot_spectrum(self, **kwargs):
		fig, ax = plt.subplots()

		# For every unique source, add their data separately
		for l in self.sed.source.unique():
			sel = self.sed[self.sed.source == l]
			ax.errorbar(sel.wavelength, sel.flux, sel.fluxerr, fmt='o', label=l, **kwargs)

		ax.legend()
		ax.set_xscale('log')

		return fig, ax
