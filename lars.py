import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes

GALAXIES = ['B1152+200', 'B1600+434', 'B1608+656', 'J0806+2006', 'J0924+0219', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'J1633+3134', 'J1650+4251']

if __name__ == '__main__':
    for g in GALAXIES:
        try:
            mags_to_fluxes(g)
        except FileNotFoundError:
            continue

        print(g)
        # TODO: flux ratios and difference between same SDSS and PanSTARRS bands

        lqso = LensedQSO(g)
        # lqso.plot_spectrum()
        lqso.plot_spectrum(loglog=True)
        lqso.plot_spectrum(mags=True)

    # lqso.filtered_sed['error %'] = lqso.filtered_sed.flux_err / lqso.filtered_sed.flux_total
    #
    # print(lqso.filtered_sed[['source', 'error %']])

    # print(lqso.sed_to_agn_fitter())

    plt.show()
