import pandas as pd

import os

sta = pd.read_csv(os.path.join('data', 'stacey_table.csv'), delimiter='\t', names=['lens', 'type', 'dtheta', 'zs', '250um', '350um', '500um', 'comments', 'references'], header=0)
cas = pd.read_csv(os.path.join('data', 'castles_table.csv'), delimiter='\t', names=['no', 'lens', 'G', 'zs', 'zl', 'RA', 'DEC', 'EBV', 'ms', 'ml', 'FGhz', 'Nim', 'size', 'dt', 'sigma'], header=0)

sta['identifier'] = sta['lens'].str.extract('[A-Z ]*([0-9]+[+−.][0-9]+)')
cas['identifier'] = cas['lens'].str.extract('[A-Z ]*([0-9]+[+-.][0-9]+)')

print(sta.merge(cas, on='identifier')[['lens_x', '250um', '350um', '500um']])
