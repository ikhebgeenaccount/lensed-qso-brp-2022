import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes

if __name__ == '__main__':
    mags_to_fluxes('B1152+200')
    lqso = LensedQSO('B1152+200')
    lqso.plot_spectrum()
    lqso.plot_spectrum(loglog=True)
    
    plt.show()
