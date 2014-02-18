#!/usr/bin/env python
# -*- coding: utf-8 -*-

from lsl.common.stations import lwa1
ants = lwa1.getAntennas()

f = open('lwa1-antennas.csv', 'w')
f.write("STAND, POSX, POSY, POSZ, POL, DELAY\n")
for a in ants:
    s = a.stand
    id, p, x, y, z = s.id, a.pol, s.x, s.y, s.z
    d = a.cable.delay(50e6)
    f.write("%s, %s, %s, %s, %s, %s\n"%(id, x, y, z, p, d))
f.close()
