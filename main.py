import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO

if __name__ == '__main__':
    lqso = LensedQSO('B1152+200')
    lqso.plot_spectrum()
    lqso.plot_spectrum_loglog()
    plt.show()
