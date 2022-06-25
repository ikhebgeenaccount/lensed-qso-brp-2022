import pandas as pd

import os

from matplotlib import pyplot as plt

from src.sed_compilation import ned_api

sta = pd.read_csv(os.path.join('data', 'stacey_table.csv'), delimiter='\t', names=['lens', 'type', 'dtheta', 'zs', '250um', '350um', '500um', 'comments', 'references'], header=0)
cas = pd.read_csv(os.path.join('data', 'castles_table.csv'), delimiter='\t', names=['no', 'lens', 'G', 'zs', 'zl', 'RA', 'DEC', 'EBV', 'ms', 'ml', 'FGhz', 'Nim', 'size', 'dt', 'sigma'], header=0)

sta['identifier'] = sta['lens'].str.extract('[A-Z ]*([0-9]+[+−.][0-9]+)')
cas['identifier'] = cas['lens'].str.extract('[A-Z ]*([0-9]+[+-.][0-9]+)')

print(sta := sta.merge(cas, on='identifier')[['lens_x', '250um', '350um', '500um']])

counts = []
for name in sta['lens_x']:
    nedsed = ned_api.access_sed(name)

    print(nedsed.shape[0])

    counts.append(nedsed.shape[0])

sta['ned_count'] = counts

print(sta[['lens_x', 'ned_count']])
print(sum(sta['ned_count'] > 0))

fig, ax = plt.subplots()

ax.hist(sta['ned_count'], bins=20)

plt.show()
