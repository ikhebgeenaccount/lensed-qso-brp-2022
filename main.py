import matplotlib.pyplot as plt

from src.fits_to_csv import sdss_fits_to_csv

gal_df = sdss_fits_to_csv('J0806+2006', 'spec-1922-53315-0582.fits')

fig, ax = plt.subplots()

ax.errorbar(gal_df['wavelength'], gal_df['flux'], yerr=gal_df['fluxerr'], fmt='o', capsize=2, capthick=1, elinewidth=1, markersize=2)

plt.show()