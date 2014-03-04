import pandas as pd
import numpy as np

a = pd.read_csv('instr_config_basic.tpl', delimiter='\t')
b = pd.read_csv('lincoln_remap.txt', delimiter='\t')

for ii in range(0, b['#remap'].shape[0], 1):
    ant = a[a['antenna'] == b['#remap'][ii]-1]
    print ant[ant['pol'] == b['pol'][ii]].values[0,3]

z = (b['#remap']-1).values
for item in z:
    print item