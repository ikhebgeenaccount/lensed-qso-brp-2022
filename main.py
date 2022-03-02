import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes

if __name__ == '__main__':
    gal_name = 'J1524+4409'
    mags_to_fluxes(gal_name)
    lqso = LensedQSO(gal_name)
    lqso.plot_spectrum()
    lqso.plot_spectrum(loglog=True)
    
    plt.show()
