import matplotlib.pyplot as plt

from src.lensed_qso import LensedQSO

if __name__ == '__main__':
    lqso = LensedQSO('J0806+2006')
    lqso.plot_spectrum()

    plt.show()
