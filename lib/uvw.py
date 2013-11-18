#! /usr/bin/env python
# encoding: utf-8
"""
uvw.py
======

Helper utilities from reading and writing JSON data with numpy arrays. The native JSON
encoder cannot handle numpy arrays and will chuck a hissy-fit, hence these helper functions.
"""

import ephem
import numpy as np

LIGHT_SPEED = 299792458 # From Google

class AntArray(ephem.Observer):
    """ An antenna array class.
    Based on pyEphem's Observer class. Probably very similar to the one in AIPY.

    Parameters
    ----------
    lat: dd:mm:ss
    latitude of array centre, e.g. 44:31:24.88
    long: dd:mm:ss
    longitude of array centre, e.g. 11:38:45.56
    elev: float
    elevation in metres of array centrem e.g. 28.0
    date: datetime object
    date and time of observation, e.g. datetime.now()
    antennas: np.array([x,y,z])
    numpy array of antenna positions, in xyz coordinates in meters,
    relative to the array centre.
    """
    def __init__(self, lat, long, elev, date, antennas):
        super(Array, self).__init__()
        self.lat = lat
        self.long = long
        self.elev = elev
        self.date = date
        self.antennas = antennas
    def update(self, date):
        """Update antenna with a new datetime"""
        self.date = date

def makeSource(name, ra, dec, flux=0, epoch=2000):
    """ Create a pyEphem FixedBody

    Parameters
    ----------
    name: string
    Name of source, e.g. CasA
    ra: hh:mm:ss
    right ascension, e.g. 23:23:26
    dec: dd:mm:ss
    declination e.g. 58:48:22.21
    flux: float
    flux brightness in Jy (not actually used here)
    epoch: J2000
    Defaults to J2000, i.e. 2000"
    """
    line = "%s,f,%s,%s,%s,%d"%(name,ra,dec,flux,epoch)
    body = ephem.readdb(line)
    return body

def computeUVW(xyz, H, d, conjugate=False, in_microseconds=True):
    """ Converts X-Y-Z coordinates into U-V-W
    Uses the transform from Thompson Moran Swenson (4.1, pg86)

    Parameters
    ----------
    xyz: should be a numpy array [x,y,z]
    H: float (degrees) is the hour angle of the phase reference position
    d: float (degrees) is the declination
    in_seconds (bool): Return in microseconds (True) or meters (False)
    Returns
    """
    sin = np.sin
    cos = np.cos

    H, d = np.deg2rad(H), np.deg2rad(d)

    xyz = np.matrix(xyz) # Cast into a matrix

    trans = np.matrix([
          [sin(H), cos(H), 0],
          [-sin(d)*cos(H), sin(d)*sin(H), cos(d)],
          [cos(d)*cos(H), -cos(d)*sin(H), sin(H)]
        ])

    uvw = trans * xyz.T
    uvw = np.array(uvw)
    if conjugate:
        uvw = -1 * uvw
    if in_microseconds:
        return uvw[:,0] / LIGHT_SPEED
    else:
        return uvw[:,0]