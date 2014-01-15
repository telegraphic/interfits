#! /usr/bin/env python
# encoding: utf-8
"""
uvw.py
======

Helper utilities from reading and writing JSON data with numpy arrays. The native JSON
encoder cannot handle numpy arrays and will chuck a hissy-fit, hence these helper functions.
"""

import time
import ephem
import numpy as np
import time, calendar
from datetime import datetime

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

def computeUVW(xyz, H, d, conjugate=False, in_microseconds=True, t_matrix=None):
    """ Converts X-Y-Z coordinates into U-V-W
    Uses the transform from Thompson Moran Swenson (4.1, pg86)

    Parameters
    ----------
    xyz: should be a numpy array [x,y,z]
    H: float (float, radians) is the hour angle of the phase reference position
    d: float (float, radians) is the declination
    in_seconds (bool): Return in microseconds (True) or meters (False)
    t_matrix (np.matrix): Default none. Can alternatively supply the transform
                          matrix to be used (useful for speed optimization).
    Returns uvw vector (in microseconds)
    """
    xyz = np.matrix(xyz) # Cast into a matrix

    if H > 2*np.pi or d > 2*np.pi:
        raise TypeError("HA and DEC are not in radians!")

    if t_matrix is None:
        sin = np.sin
        cos = np.cos

        #H, d = np.deg2rad(H), np.deg2rad(d)

        trans = np.matrix([
              [sin(H), cos(H), 0],
              [-sin(d)*cos(H), sin(d)*sin(H), cos(d)],
              [cos(d)*cos(H), -cos(d)*sin(H), sin(d)]
            ])
    else:
        trans = t_matrix

    uvw = trans * xyz.T
    uvw = np.array(uvw)
    if conjugate:
        uvw = -1 * uvw
    if in_microseconds:
        return uvw[:,0] / LIGHT_SPEED
    else:
        return uvw[:,0]


def convertToJulianTuple(timestamp):
    """ Convert a list of timestamps into DATE and TIME since julian midnight

    FITS-IDI requires DATE parameter is set to Julian date @ midnight. TIME is
    then indexed against this DATE.
    
    timestamp (float): timestamp of type 'float' is preferred, but this should
                       handle time tuples, strings and datetime objects too.
    """
    
    # Figure out exactly what type of timestamp we've got.
    if type(timestamp) in (str, unicode):
        tt = time.strptime(timestamp, "%Y-%m-%dT%H:%M:%S")
        ts = calendar.timegm(tt)
    elif type(timestamp) is float:
        ts = timestamp
    elif type(timestamp) is tuple:
        ts = calendar.timegm(timestamp)
    elif type(timestamp) is type(datetime.now):
        date_str = timestamp.strftime("%Y-%m-%dT%H:%M:%S")
        tt = time.strptime(date_str, "%Y-%m-%dT%H:%M:%S")
        ts = calendar.timegm(tt)
    else:
        print type(timestamp)
        raise TypeError
        
    # DATE is julian date at midnight that day
    # TIME is in DAYS since midnight
    julian = ephem.julian_date(time.gmtime(ts)[:6])
    # Ephem returns julian date at NOON, we need at MIDNIGHT
    julian_midnight = int(julian) - 0.5

    time_elapsed = ephem.julian_date(time.gmtime(ts)[:6]) - julian_midnight

    return julian_midnight, time_elapsed
