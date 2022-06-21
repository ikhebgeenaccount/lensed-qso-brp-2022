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
def _file_reader(files, cols_orig, cols_res=None, join_col=None, open_kwargs=None, post_filter=None, **kwargs):
    if open_kwargs is None:
        open_kwargs = {}
    if cols_res is None:
        cols_res = COLUMNS

    df = None
    for file in files:
        ext = file.split('.')[-1]
        if ext == 'fits':
            tdf = _fits_reader(file, fits_open_kwargs=open_kwargs, **kwargs)
        elif ext == 'csv':
            tdf = _csv_reader(file, read_csv_kwargs=open_kwargs)
        else:
            raise ValueError('Extension ' + ext + ' does not have an associated file reader.')

        if df is None:
            df = tdf
        else:
            df = df.merge(tdf, on=join_col, how='left')

    if callable(post_filter):
        df = post_filter(df)

    return _extract_cols(df, cols_orig, cols_res)


def _fits_reader(file, fits_open_kwargs=None, **kwargs):
    if 'hdu_index' in kwargs:
        hdu_index = kwargs['hdu_index']
    else:
        hdu_index = 1

    with fits.open(file, **fits_open_kwargs) as hdul:
        df = Table(hdul[hdu_index].data).to_pandas()

    return df


def _csv_reader(file, read_csv_kwargs=None):
    return pd.read_csv(file, **read_csv_kwargs)


def _extract_cols(df, col_names, cols_res):
    """
    Extract the columns that are defined in COLUMNS, using the col_names translation from the given df.
    """
    rdf = pd.DataFrame()
    for i, c in enumerate(cols_res):
        if callable(col_names[i]):
            rdf[c] = col_names[i](df)
        else:
            rdf[c] = df[col_names[i]]

    return rdf


COLUMNS = ['name', 'redshift', 'logSFR', 'logSFR_pe', 'logSFR_me', 'logMstar', 'logMstar_pe', 'logMstar_me']

FILES = {
    'SMGs, Birkin+2020, z=[1.2-6.1]' : _file_reader([os.path.join('data', 'context_main_seq', 'table_a2.fits'),
                                      os.path.join('data', 'context_main_seq', 'table_a1.fits')],
                                     join_col='ID', cols_orig=['ID', 'z_phot', 'SFR', lambda df: df['SFR_u'] - df['SFR'],
                                                               lambda df: df['SFR'] - df['SFR_l'], 'Mstar',
                                                               lambda df: df['Mstar_u'] - df['Mstar'], lambda df: df['Mstar'] - df['Mstar_l']]),
    'SMGs, Da Cunha+2015, z=[1.3-6.1]': _file_reader([os.path.join('data', 'context_main_seq', 'Cunha.csv')], cols_orig=['ID', 'z_phot'] + COLUMNS[2:],
                                                    open_kwargs={'delimiter': '\t'}),
    # 'ULIRGs, Da Cunha+2010, z=[0.03-0.5]': _file_reader([os.path.join('data', 'context_main_seq', 'ulirgs_cunha1.csv'),
    #                                     os.path.join('data', 'context_main_seq', 'ulirgs_cunha2.csv')], join_col='Galaxy',
    #                                                    cols_orig=['Galaxy', 'z', lambda df: df['logsSFR'] + df['logMstar'],
    #                                                               lambda df: np.zeros(len(df['Galaxy'])), lambda df: np.zeros(len(df['Galaxy']))] +
    #                                   COLUMNS[-3:], open_kwargs={'delimiter': '\t'}),
    'Local SFGs, Sun+2020, z=0': _file_reader([os.path.join('data', 'context_main_seq', 'sun.csv')],
                                             cols_orig=['galaxy', lambda df: df['d'] * 70 / 3e5, lambda df: np.log10(df['SFR']),
                                                        lambda df: np.zeros(len(df['galaxy'])), lambda df: np.zeros(len(df['galaxy'])),
                                                        lambda df: np.log10(df['Mstar'] * 1e9),
                                                        lambda df: np.zeros(len(df['galaxy'])), lambda df: np.zeros(len(df['galaxy']))],
                                             open_kwargs={'delim_whitespace': True}),
    'SFGs, Tacconi+2013, z=[1-2.4]': _file_reader([os.path.join('data', 'context_main_seq', 'PHIBBS_table1.csv'),
                                       os.path.join('data', 'context_main_seq', 'PHIBBS_table2.csv')], join_col='Source',
                                                 cols_orig=['Source', 'z_CO', lambda df: np.log10(df['SFR^d']), lambda df: .35 * df['SFR^d'] / (df['SFR^d'] * np.log(10.)),
                                                            lambda df: .35 * df['SFR^d'] / (df['SFR^d'] * np.log(10.)), lambda df: np.log10(df['M_*^g']),
                                                            lambda df: .3 * df['M_*^g'] / (df['M_*^g'] * np.log(10.)), lambda df: .3 * df['M_*^g'] / (df['M_*^g'] * np.log(10.))],
                                                 open_kwargs={'delimiter': '\t', 'skipinitialspace': True, 'na_values': ['...', '... ']}),


    # ULIRGs from Cunha+2010 don't work: they list specific star formation rate instead of star formation rate.

    # COSMOS is a lot
    # 'COSMOS, Laigle+2016': _file_reader([os.path.join('data', 'context_main_seq', 'COSMOS2015_Laigle+_v1.1.fits')],
    #                                     col_names=['NUMBER', 'ZPDF', 'SFR_MED', 'SFR_MED_MAX68', 'SFR_MED_MIN68',
    #                                                'MASS_MED', 'MASS_MED_MAX68', 'MASS_MED_MIN68'],
    #                                     fits_open_kwargs={'mmap': False},
    #                                     post_filter=lambda df: df[(df['TYPE'] == 0) & (df['ZPDF'] > .9) & (df['ZPDF'] < 1.6)
    #                                                               & (df['SFR_MED'] > 0.) & (df['MASS_MED'] > 6.)]),
    'Quasars, Jarvis+2020, z=[0.1-0.2]': _file_reader([os.path.join('data', 'context_main_seq', 'qsos_jarvis.csv')],
                                                  cols_orig=['Name', 'z', lambda df: np.log10(df['SFR']), lambda df: df['SFR_pe'] / (df['SFR'] * np.log(10.)),
                                                             lambda df: df['SFR_me'] / (df['SFR'] * np.log(10.))] + COLUMNS[-3:],
                                                  open_kwargs={'delimiter': '\t'}),
    # 'Local catalog, Salim+2018' : _file_reader([os.path.join('data', 'context_main_seq', 'GSWLC-D2.dat')],
    #                                           col_names=['objid', 'z', 'logSFR', 'logSFR_err', 'logSFR_err', 'logMstar', 'logMstar_err', 'logMstar_err'],
    #                                           open_kwargs={'delim_whitespace': True, 'index_col': False},
    #                                           post_filter=lambda df: df[np.random.choice(a=[True, False], size=len(df), p=[1000/50000, 1-1000/50000]) &
    #                                                                     (df['logMstar'] > -99)]),
    'AGN, Kakkad+2017, z=1.5': _file_reader([os.path.join('data', 'context_main_seq', 'agns_kakkad.csv')],
                                            cols_orig=['ID', 'za', lambda df: np.log10(df['SFRe']), lambda df: np.zeros(len(df)),
                                                       lambda df: np.zeros(len(df)), 'Log Mc', lambda df: np.zeros(len(df)), lambda df: np.zeros(len(df))],
                                            open_kwargs={'delimiter': '\t'})
}
