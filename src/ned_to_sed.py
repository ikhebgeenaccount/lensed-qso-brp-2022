import pandas as pd
import numpy as np
import os.path
import warnings

from src.lensed_qso import LensedQSO

def ned_table_to_sed(galaxy, ned_file='ned.txt', wavelength_conversion=1e4, flux_conversion=1e3):
    ned_df = pd.read_csv(os.path.join('data', galaxy, ned_file), delimiter='|')
    
    # Filter out non-measurement
    fs = np.sum(ned_df['Flux Density'].isna())
    if fs > 0:
        warnings.warn(f'Filtered out {fs} measurements, check if this is correct!')
        ned_df = ned_df[ned_df['Flux Density'].notna()]    
    
    if not np.all(np.where(ned_df['Upper limit of uncertainty'] == ned_df['Lower limit of uncertainty'], 1, 0)):
        # print(ned_df[['Upper limit of uncertainty', 'Lower limit of uncertainty']])
        raise NotImplementedError('Upper and lower limit of uncertainty are different, this cannot be handled.')
        
    if np.all(np.where(ned_df['Upper limit of Flux Density'] == np.nan, 1, 0)) or np.all(np.where(ned_df['Lower limit of Flux Density'] == np.nan, 1, 0)):
        raise NotImplementedError('Table has upper or lower limit for flux density, this cannot be handled.')

    ned_df['wavelength'] = ned_df['Wavelength']
    ned_df['flux'] = ned_df['Flux Density']
    ned_df['flux_err'] = ned_df['Upper limit of uncertainty']
    ned_df['observed_passband'] = ned_df['Observed Passband']
    
    # TODO: convert wavelength, flux, flux_err
    # TODO: strip strings? check if spaces in observed_passband
    
    lqso = LensedQSO(galaxy)
    
    print(ned_df[['wavelength', 'flux', 'flux_err', 'observed_passband']])
    
    print(lqso.sed.size)
    
    # TODO: append doesn' t work?
    lqso.sed.append(ned_df[['wavelength', 'flux', 'flux_err', 'observed_passband']])
    
    print(lqso.sed.size)
