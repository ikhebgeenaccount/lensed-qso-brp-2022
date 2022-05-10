from astropy.io import fits
from astropy.table import Table

import pandas as pd
import numpy as np
import os


"""
The _reader functions return a pandas.Dataframe with columns
['name', 'z', 'SFR', 'SFR_pe', 'SFR_me', 'Mstar', 'Mstar_pe', 'Mstar_me']
in units SFR: [M_sun/yr]; Mstar: [M_sun].
(pe = plus error, me = min error, i.e. SFR^{+SFR_pe}_{-SFR_me})

file must be the complete path to the file to be read.

col_names must be an array of length 8, which is equal to the amount of columns in the resulting DataFrame.
col_names is used to translate from the file that is being read to create the DataFrame, which means that the
order of names in col_names must obey the order of names in the list above. col_names can contain strings for
columns that only need a namechange, or callables (functions/lamba expressions) for columns that require some
arithmetic or other change.
"""
def _fits_reader(files, col_names, join_col=None, hdu_index=1):
    df = None
    for file in files:
        with fits.open(file) as hdul:

            if df is None:
                df = Table(hdul[hdu_index].data).to_pandas()
            else:
                odf = Table(hdul[hdu_index].data).to_pandas()
                df = df.merge(odf, on=join_col, how='left')

    return _extract_cols(df, col_names)


def _csv_reader(files, col_names, join_col=None, read_csv_kwargs=None):
    if read_csv_kwargs is None:
        read_csv_kwargs = {}

    df = None
    for file in files:
        tdf= pd.read_csv(file, **read_csv_kwargs)

        if df is None:
            df = tdf
        else:
            df = df.merge(tdf, on=join_col, how='left')

    return _extract_cols(df, col_names)


def _extract_cols(df, col_names):
    """
    Extract the columns that are defined in COLUMNS, using the col_names translation from the given df.
    """
    rdf = pd.DataFrame()
    for i, c in enumerate(COLUMNS):
        if callable(col_names[i]):
            rdf[c] = col_names[i](df)
        else:
            rdf[c] = df[col_names[i]]

    return rdf


COLUMNS = ['name', 'redshift', 'logSFR', 'logSFR_pe', 'logSFR_me', 'logMstar', 'logMstar_pe', 'logMstar_me']

FILES = {
    'SMGs, Birkin+2020': _fits_reader([os.path.join('data', 'context_main_seq', 'table_a2.fits'),
                                      os.path.join('data', 'context_main_seq', 'table_a1.fits')],
                                     join_col='ID', col_names=['ID', 'z_phot', 'SFR', lambda df: df['SFR_u'] - df['SFR'],
                                                               lambda df: df['SFR'] - df['SFR_l'], 'Mstar',
                                                               lambda df: df['Mstar_u'] - df['Mstar'], lambda df: df['Mstar'] - df['Mstar_l']]),
    'SFGs, Tacconi+2013': _csv_reader([os.path.join('data', 'context_main_seq', 'PHIBBS_table1.csv'),
                                       os.path.join('data', 'context_main_seq', 'PHIBBS_table2.csv')], join_col='Source',
                                      col_names=['Source', 'z_CO', lambda df: np.log10(df['SFR^d']), lambda df: np.zeros(len(df['Source'])),
                                                 lambda df: np.zeros(len(df['Source'])), lambda df: np.log10(df['M_*^g']),
                                                 lambda df: np.zeros(len(df['Source'])), lambda df: np.zeros(len(df['Source']))],
                                      read_csv_kwargs={'delimiter': '\t', 'skipinitialspace': True, 'na_values': ['...', '... ']}),
    'SMGs, Cunha+2015': _csv_reader([os.path.join('data', 'context_main_seq', 'Cunha.csv')], col_names=['ID', 'z_phot'] + COLUMNS[2:],
                                    read_csv_kwargs={'delimiter': '\t'}),
}
