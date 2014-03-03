#! /usr/bin/env python
# -*- coding: utf-8 -*-
# encoding: utf-8
"""
coords.py
=========

Coordinate and time transforms required for generating correct timestamps and UVW
coordinates.

"""

import time
import re
import ephem
import numpy as np
import time, calendar
from datetime import datetime

from numpy import sin, cos

LIGHT_SPEED = 299792458.0 # From Google

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

def generateBaselineIds(n_ant=32, autocorrs=True):
    """ Generate a list of unique baseline IDs and antenna pairs

    This uses the MIRIAD definition for >256 antennas:
    bl_id = 2048*ant1 + ant2 + 65536

    n_ant: number of antennas in the array
    autocorrs: include antenna autocorrelations?

    TODO: Rewrite for speed (itertools?)
    """

    bls, ant_arr = [], []
    for ii in range(1, n_ant + 1):
        for jj in range(1, n_ant + 1):
            if jj >= ii:
                if autocorrs is False and ii == jj:
                    pass
                else:
                    ant_arr.append((ii, jj))
                    if ii > 255 or jj > 255:
                        bl_id = ii * 2048 + jj + 65536
                    else:
                        bl_id = 256 * ii + jj
                    bls.append(bl_id)
    return bls, ant_arr

def computeUVW(xyz, H, d, in_seconds=True, conjugate=False):
    """ Converts X-Y-Z baselines into U-V-W

    Parameters
    ----------
    xyz: should be a numpy array [x,y,z] of baselines (NOT ANTENNAS!)
    H: float (float, radians) is the hour angle of the phase reference position
    d: float (float, radians) is the declination
    conjugate: (bool): Conjugate UVW coordinates?
    in_seconds (bool): Return in seconds (True) or meters (False)
    Returns uvw vector (in microseconds)

    Notes
    -----
    The transformation is precisely that from the matrix relation 4.1 in
    Thompson Moran Swenson:
              [sin(H), cos(H), 0],
              [-sin(d)*cos(H), sin(d)*sin(H), cos(d)],
              [cos(d)*cos(H), -cos(d)*sin(H), sin(d)]

    A right-handed Cartesian coordinate system is used where X
    and Y are measured in a plane parallel to the earth's equator, X in the meridian
    plane (defined as the plane through the poles of the earth and the reference point
    in the array), Y is measured toward the east, and Z toward the north pole. In terms
    of hour angle H and declination d, the coordinates (X, Y, Z) are measured toward
    (H = 0, d = 0), (H = -6h, d = O), and (d = 90"), respectively.

    Here (H, d) are usually the hour angle and declination of the phase reference
    position.

    """
    is_list = True
    if type(xyz) in (list, tuple):
        x, y, z = xyz
    if type(xyz) == type(np.array([1])):
        is_list = False
        try:
            #print xyz.shape
            x, y, z = np.split(xyz, 3, axis=1)
        except:
            print xyz.shape
            raise

    sh, sd = sin(H), sin(d)
    ch, cd = cos(H), cos(d)
    u  = sh * x + ch * y
    v  = -sd * ch * x + sd * sh * y + cd * z
    w  = cd * ch * x - cd * sh * y + sd * z

    if is_list:
        uvw = np.array((u, v, w))
    else:
        uvw = np.column_stack((u, v, w))

    if conjugate:
        uvw *= -1
    if in_seconds:
        return uvw / LIGHT_SPEED
    else:
        return uvw

def computeBaselineVectors(xyz, autocorrs=True):
    """ Compute all the possible baseline combos for antenna array

    Parameters
    ----------
    xyz: np.array of antenna positions.
    autocorrs: (bool) baselines should contain autocorrs (zero-length baselines)

    Returns
    -------
    XYZ array of all antenna combinations, in ascending antenna IDs.
    """
    try:
        assert type(xyz) == type(np.array([1]))
    except AssertionError:
        xyz = np.array(xyz)
        if xyz.shape[1] != 3:
            raise ValueError("Cannot understand input XYZ array")

    n_ant = xyz.shape[0]
    if autocorrs:
        bls = np.zeros([n_ant * (n_ant - 1) / 2 + n_ant, 3])
    else:
        bls = np.zeros([n_ant * (n_ant - 1) / 2, 3])

    bls_idx = 0
    for ii in range(n_ant):
        for jj in range(n_ant):
            if jj >= ii:
                if autocorrs is False and ii == jj:
                    pass
                else:
                    bls[bls_idx] = xyz[ii] - xyz[jj]
                    bls_idx += 1
    return bls

def coordTransform(xyz, input='ENU', output='NED'):
    """ Coordinate frame transformations.

    Understood frames are right-hand-rule based:
    ENU - East North Up
    NED - North East Down
    """

    is_list = True
    if type(xyz) in (list, tuple):
        a, b, c = xyz
    if type(xyz) == type(np.array([1])):
        is_list = False
        try:
            a, b, c = np.split(xyz, 3, axis=1)
        except:
            print xyz.shape
            raise

    if input == 'ENU':
        x, y, z = a, b, c
    elif input == 'NED':
        x, y, z = b, a, -c
    else:
        raise CoordinateError("Unknown coordinate system: %s"%input)

    if output == 'ENU':
        if is_list:
            return np.array((x, y, z))
        else:
            return np.column_stack((x, y, z))
    elif output == 'NED':
        if is_list:
            return np.array((y, x, -z))
        else:
            return np.column_stack((y, x, -z))
    else:
        raise CoordinateError("Unknown coordinate system: %s"%output)


def geo2ecef(lat, lon, elev):
    """
    Convert latitude (rad), longitude (rad), elevation (m) to earth-
    centered, earth-fixed coordinates.

    Note: This routine from LSL 1.0.0
    http://fornax.phys.unm.edu/lwa/trac/wiki
    """

    WGS84_a = 6378137.000000
    WGS84_b = 6356752.314245
    N = WGS84_a**2 / np.sqrt(WGS84_a**2*np.cos(lat)**2 + WGS84_b**2*np.sin(lat)**2)

    x = (N+elev)*np.cos(lat)*np.cos(lon)
    y = (N+elev)*np.cos(lat)*np.sin(lon)
    z = ((WGS84_b**2/WGS84_a**2)*N+elev)*np.sin(lat)

    return (x, y, z)

def ecef2geo(x, y, z):
    """
    Convert earth-centered, earth-fixed coordinates to (rad), longitude
    (rad), elevation (m) using Bowring's method.

    Note: This routine from LSL 1.0.0
    http://fornax.phys.unm.edu/lwa/trac/wiki
    """

    WGS84_a = 6378137.000000
    WGS84_b = 6356752.314245
    e2 = (WGS84_a**2 - WGS84_b**2) / WGS84_a**2
    ep2 = (WGS84_a**2 - WGS84_b**2) / WGS84_b**2

    # Distance from rotation axis
    p = np.sqrt(x**2 + y**2)

    # Longitude
    lon = np.arctan2(y, x)
    p = np.sqrt(x**2 + y**2)

    # Latitude (first approximation)
    lat = np.arctan2(z, p)

    # Latitude (refined using Bowring's method)
    psi = np.arctan2(WGS84_a*z, WGS84_b*p)
    num = z + WGS84_b*ep2*np.sin(psi)**3
    den = p - WGS84_a*e2*np.cos(psi)**3
    lat = np.arctan2(num, den)

    # Elevation
    N = WGS84_a**2 / np.sqrt(WGS84_a**2*np.cos(lat)**2 + WGS84_b**2*np.sin(lat)**2)
    elev = p / np.cos(lat) - N

    return lat, lon, elev

def convertToJulianTuple(timestamp):
    """ Convert a list of timestamps into DATE and TIME since julian midnight

    FITS-IDI requires DATE parameter is set to Julian date @ midnight. TIME is
    then indexed against this DATE.
    
    timestamp (float): timestamp of type 'float' is preferred, but this should
                       handle time tuples, strings and datetime objects too.
    """
    
    # Figure out exactly what type of timestamp we've got.
    if type(timestamp) in (str, unicode):
        tt = parse_timestring(timestamp)
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
        raise TypeError("Unknown timestamp type '%s'" % str(type(timestamp)))
        
    # DATE is julian date at midnight that day
    # TIME is in DAYS since midnight
    julian = ephem.julian_date(time.gmtime(ts)[:6])
    # Ephem returns julian date at NOON, we need at MIDNIGHT
    julian_midnight = int(julian) - 0.5

    time_elapsed = ephem.julian_date(time.gmtime(ts)[:6]) - julian_midnight

    return julian_midnight, time_elapsed

def parse_timestring(tstring):
    """ Convert timestring into timestamp """
    try:
        if re.match("\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d.\d+$", tstring):
            return time.strptime(tstring, "%Y-%m-%dT%H:%M:%S.%f")
        elif re.match("\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d.\d$", tstring):
            return time.strptime(tstring, "%Y-%m-%dT%H:%M:%S.%f")
        elif re.match("\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d$", tstring):
            return time.strptime(tstring, "%Y-%m-%dT%H:%M:%S")
        else:
            raise ValueError("Cannot parse %s"%tstring)
    except ValueError:
        raise ValueError("Unable to parse timestring '%s'" % tstring)
