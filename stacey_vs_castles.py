import pandas as pd

import os

from matplotlib import pyplot as plt

from src.lensed_qso import LensedQSO
from src.sed_compilation import ned_api, castles, ned

# Print all rows and columns
pd.set_option('display.max_rows', None)
pd.set_option('max_columns', None)

# Read Stacey and CASTLES tables
sta = pd.read_csv(os.path.join('data', 'stacey_table.csv'), delimiter='\t', names=['lens', 'type', 'dtheta', 'zs', '250um', '350um', '500um', 'comments', 'references'], header=0)
cas = pd.read_csv(os.path.join('data', 'castles_table.csv'), delimiter='\t', names=['no', 'lens', 'G', 'zs', 'zl', 'RA', 'DEC', 'EBV', 'ms', 'ml', 'FGhz', 'Nim', 'size', 'dt', 'sigma'], header=0)

# Remove leading and trailing whitespace from names
sta['lens'] = sta['lens'].str.strip()
cas['lens'] = cas['lens'].str.strip()

# Make sure we have something that is equal for both, to merge
sta['identifier'] = sta['lens'].str.extract('[A-Z ]*([0-9]+[+−.][0-9]+)')
cas['identifier'] = cas['lens'].str.extract('[A-Z ]*([0-9]+[+-.][0-9]+)')

# Merge the DataFrames
df_comb = sta.merge(cas, on='identifier', suffixes=['_sta', '_cas'])

print(df_comb[['lens_sta', '250um', '350um', '500um']])

print(df_comb.columns)

# Find NED data point counts for each lqso
counts = []
for cas_name, sta_name in zip(df_comb['lens_cas'], df_comb['lens_sta']):
    lqso = LensedQSO(cas_name, aliases=[sta_name])
    print(sta_name)
    # Fetch CASTLES data
    castles.update_mags(lqso)

    # NED data
    ned.update_sed(lqso)

df_comb['ned_count'] = counts

print(df_comb[['lens_sta', 'lens_cas', 'ned_count']].sort_values(by='ned_count'))
print(sum(df_comb['ned_count'] > 0))

# Histogram of NED data point counts
fig, ax = plt.subplots()

ax.hist(df_comb['ned_count'], bins=20)

# IDEA: foreground subtraction: for systems that have no foreground data points, take average of all other lensing galaxy, then scale to luminosity distance

plt.show()
