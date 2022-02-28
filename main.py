import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes

if __name__ == '__main__':
    lqso = LensedQSO('J1633+3134')
    lqso.plot_spectrum()
    lqso.plot_spectrum(loglog=True)
    
    mags_to_fluxes('J1633+3134')
    plt.show()
