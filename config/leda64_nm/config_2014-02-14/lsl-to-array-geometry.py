import pandas as pd
import ujson as json
SOL = 299792458

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
    stabxyz.append(list(np.round(srow, 3)))
    
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