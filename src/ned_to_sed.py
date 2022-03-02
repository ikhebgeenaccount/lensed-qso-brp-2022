import pandas as pd
import numpy as np
import os.path
import warnings

from src.lensed_qso import LensedQSO


def ned_table_to_sed(galaxy, ned_file='ned.txt', wavelength_conversion=1e4, flux_conversion=1e3):
    """
    Reads a ned.txt file of a NED bar-separated table and enters it into the SED of the galaxy.
    :param galaxy:
    :param ned_file:
    :param wavelength_conversion: default: um to Angstrom
    :param flux_conversion: default: Jy to mJy
    :return:
    """
    ned_df = pd.read_csv(os.path.join('data', galaxy, ned_file), delimiter='|')
    
    # Filter out non-measurement
    fs = np.sum(ned_df['Flux Density'].isna())
    if fs > 0:
        warnings.warn('Filtered out ' + str(fs) + ' measurements, check if this is correct!')
        ned_df = ned_df[ned_df['Flux Density'].notna()]    
    
    if not np.all(np.where(ned_df['Upper limit of uncertainty'] == ned_df['Lower limit of uncertainty'], 1, 0)):
        # print(ned_df[['Upper limit of uncertainty', 'Lower limit of uncertainty']])
        raise NotImplementedError('Upper and lower limit of uncertainty are different, this cannot be handled.')
        
    if np.all(np.where(ned_df['Upper limit of Flux Density'] == np.nan, 1, 0)) or np.all(np.where(ned_df['Lower limit of Flux Density'] == np.nan, 1, 0)):
        raise NotImplementedError('Table has upper or lower limit for flux density, this cannot be handled.')

    # Give fields we want proper name, convert
    ned_df['wavelength'] = ned_df['Wavelength'] * wavelength_conversion
    ned_df['flux_total'] = ned_df['Flux Density'] * flux_conversion
    ned_df['flux_err'] = ned_df['Upper limit of uncertainty']
    ned_df['observed_passband'] = ned_df['Observed Passband'].apply(lambda v: v.strip())

    wls = ned_df['wavelength'].unique()

    lqso = LensedQSO(galaxy)

    # Check if observed_passband column exists, if not add
    if not 'observed_passband' in lqso.sed.columns:
        lqso.sed['observed_passband'] = ''

    # Combine measurements for each unique wavelength into single value and error
    for wl in wls:
        sel = ned_df.loc[ned_df['wavelength'] == wl]
        weights = 1. / np.power(sel['flux_err'], 2)

        flux_total = np.average(sel['flux_total'], weights=weights)
        flux_err = np.sqrt(1. / np.sum(weights))  # TODO: sanity check error calculation
        observed_passband = sel['observed_passband'].values[0]

        # Check if entry is not in lqso.sed yet
        if lqso.sed[(lqso.sed['wavelength'] == wl) & (lqso.sed['observed_passband'] == observed_passband)].empty:
            # Add it
            lqso.sed = lqso.sed.append({
                'wavelength': wl,
                'flux_total': flux_total,  # weighted average
                'flux_err': flux_err,
                'observed_passband': observed_passband
            }, ignore_index=True)
        else:
            # Overwrite it
            index = lqso.sed.index[(lqso.sed['wavelength'] == wl) & (lqso.sed['observed_passband'] == observed_passband)]
            lqso.sed.loc[index, ['wavelength', 'flux_total', 'flux_err', 'observed_passband']] = \
                [wl, flux_total, flux_err, observed_passband]

    # Save to SED
    lqso.save_sed()
