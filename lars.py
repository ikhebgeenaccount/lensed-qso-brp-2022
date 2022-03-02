import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.lensed_qso import LensedQSO
from src.mags_to_fluxes import mags_to_fluxes

# ratio avg: 1.4543629655671957, ratio std: 0.37071640199561534

GALAXIES = ['B1152+200', 'B1600+434', 'B1608+656', 'J0806+2006', 'J0924+0219', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'J1633+3134', 'J1650+4251']
SDSS_V_PANSTARRS_GS = ['B1152+200', 'B1600+434', 'J0806+2006', 'J1633+3134', 'J1455+1447', 'J1524+4409']

if __name__ == '__main__':    
    flux_discrepancy = pd.DataFrame(columns=['galaxy', 'filter', 'sdss', 'panstarrs', 'ratio', 'diff'])
    
    ratios = []
    diffs = []
    
    panstarrs_mag = 'ap'
    
    for g in SDSS_V_PANSTARRS_GS:
        try:
            mags_to_fluxes(g)
        except FileNotFoundError:
            continue

        print(f'{g}, PanSTARRS mag choice: {panstarrs_mag}')
        # TODO: flux ratios and difference between same SDSS and PanSTARRS bands
        
        lqso = LensedQSO(g)
        # lqso.plot_spectrum()
        lqso.plot_spectrum(loglog=True)
        lqso.plot_spectrum(mags=True)
        
        for i, r in lqso.sed.iterrows():
            if 'SDSS' not in r['source'] and panstarrs_mag not in r['source']:
                continue
            if r['filter'] in ['g', 'r', 'i', 'z']:
                source = 'sdss' if 'SDSS' in r.source else 'panstarrs'
                flux_discrepancy.loc[r['filter'], source] = r.flux_total
        
        flux_discrepancy['ratio'] = flux_discrepancy.sdss.div(flux_discrepancy.panstarrs.where(flux_discrepancy.panstarrs != 0, np.nan))
        # Make all ratios > 1 so we can compare
        flux_discrepancy['ratio'] = flux_discrepancy['ratio'].apply(lambda r: 1. / r if r < 1  else r)
        flux_discrepancy['diff'] = flux_discrepancy['sdss'] - flux_discrepancy['panstarrs']
        print(flux_discrepancy)
        print()
        
        ratios.append(list(flux_discrepancy['ratio']))
        diffs.append(list(flux_discrepancy['diff']))
    
    ratios = np.array(ratios).flatten()
    diffs= np.array(diffs).flatten()
    
    print(f'ratio avg: {np.nanmean(ratios)}, ratio std: {np.nanstd(ratios)}')
    print(f'diff avg: {np.nanmean(diffs)}, diff std: {np.nanstd(diffs)}')
            

    # lqso.filtered_sed['error %'] = lqso.filtered_sed.flux_err / lqso.filtered_sed.flux_total
    #
    # print(lqso.filtered_sed[['source', 'error %']])

    # print(lqso.sed_to_agn_fitter())

    plt.show()
