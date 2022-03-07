import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes
from src.ned_to_sed import ned_table_to_sed

if __name__ == '__main__':
    galaxy = 'J1650+4251'
    lqso = LensedQSO(galaxy)
    
    mags_to_fluxes(lqso)
    ned_table_to_sed(lqso,'ned_galex', allowed_sources=['Chandra', 'WISE', '2MASS', 'Galex'])
    lqso.plot_spectrum(loglog=True)
    
    plt.show()
