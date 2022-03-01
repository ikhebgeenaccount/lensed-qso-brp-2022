import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes

if __name__ == '__main__':
    mags_to_fluxes('J1524+4409')
    lqso = LensedQSO('J1524+4409')
    lqso.plot_spectrum()
    lqso.plot_spectrum(loglog=True)
    
    plt.show()
