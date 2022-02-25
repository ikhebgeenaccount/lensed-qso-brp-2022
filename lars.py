import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO

if __name__ == '__main__':
    lqso = LensedQSO('J0806+2006')
    lqso.plot_spectrum()
    lqso.plot_spectrum(loglog=True)

    lqso.filtered_sed['error %'] = lqso.filtered_sed.flux_err / lqso.filtered_sed.flux_total

    print(lqso.filtered_sed[['source', 'error %']])

    print(lqso.sed_to_agn_fitter())

    plt.show()
