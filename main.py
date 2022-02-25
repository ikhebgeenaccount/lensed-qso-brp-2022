import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO

if __name__ == '__main__':
    lqso = LensedQSO('J1633+3134')
    lqso.plot_spectrum()
    lqso.plot_spectrum(loglog=True)
    plt.show()
