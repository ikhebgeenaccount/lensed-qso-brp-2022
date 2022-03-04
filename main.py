import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes
from src.ned_to_sed import ned_table_to_sed

if __name__ == '__main__':
    galaxy = 'J0924+0219'
    lqso = LensedQSO(galaxy)
    
    mags_to_fluxes(lqso)
    ned_table_to_sed(lqso,'NED_11', allowed_sources=['Chandra', 'WISE', '2MASS'])
    lqso.plot_spectrum(loglog=True)
    
    plt.show()
