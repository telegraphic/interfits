#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import ujson as json
import numpy as np
from lsl.common.stations import lwa1
SOL = 299792458

toECI = lwa1.getECITransform()


a   = pd.read_csv('leda64nm-mapping.csv')
ant = pd.read_csv('leda64nm-antennas.csv')

a_sorted = a.sort('LEDA_input', ascending=1)
ant_sorted = ant.sort('STAND', ascending=1)

stands = a_sorted['stand']
stand_names = ["stand_%i"%s for s in stands.values]

stabxyz = []
delays = []

for s in stands:
    srow = ant[ant['STAND'] == s].values[0, 1:4]
    
    # Convert from topocentric to ECI coordinates for the stands
    topo = np.round(srow, 3)
    topo.shape += (1,)
    eci = toECI*np.matrix(topo)
    eci = np.array(eci).flatten()
    stabxyz.append(list(np.round(eci, 3)))
    
    # Adding in delays
    sd = ant[ant['STAND'] == s][' DELAY'].values * SOL
    sd = sd.tolist()
    delays.append(sd)
    
print stabxyz
print stand_names
print delays

json.dump(stand_names, open('stands-temp.json', 'w'))
json.dump(stabxyz, open('stabxyz-temp.json', 'w'))
json.dump(delays, open('delays-temp.json', 'w'))