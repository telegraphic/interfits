# -*- coding: utf-8 -*-

"""
interfits.py
============

Python class which reads a variety of different visibility formats.
Currently reads UV-FITS and FITS-IDI.
"""


import re
import os
import sys
from datetime import datetime
import h5py
import ephem
from colorama import Fore
import pyfits as pf
import numpy as np
from lxml import etree

from lib.pyFitsidi import *
from lib.json_numpy import *
from lib import coords


__version__ = '0.0'
__all__ = ['LinePrint', 'PrintLog', 'InterFits', '__version__', '__all__']


class LinePrint():
    """
    Print things to stdout on one line dynamically
    """

    def __init__(self, data):
        sys.stdout.write("\r\x1b[K" + data.__str__())
        sys.stdout.flush()


class PrintLog(object):
    """
    Print / log based on verbosity level

    Verbosity level: 0 - 4
    0: Print nothing
    1: Print H1
    2: Print H1 + H2
    3: Print H1 + H2 + H3
    4: Print all
    5: Print all + debugging info
    """

    def __init__(self, verbosity=4):
        self.vlevel = verbosity

    def h1(self, headstr):
        """ Print a string as a header """
        if self.vlevel >= 1:
            print(Fore.GREEN + '\n%s' % headstr)
            line = ""
            for ii in range(len(headstr)):
                line += "-"
            print(Fore.GREEN + line)
            print(Fore.WHITE)

    def h2(self, headstr):
        """ Print a string as a header """
        if self.vlevel >= 2:
            print '\n###  ', headstr

    def h3(self, headstr):
        """ Print a string as a 3rd level header """
        if self.vlevel >= 3:
            print "\t", headstr

    def pp(self, text):
        """ Print a text string """
        if self.vlevel >= 4:
            print text

    def pa(self, text):
        """ Always print """
        print text

    def debug(self, text):
        """ Print debug information to screen """
        if self.vlevel >= 5:
            print text

    def err(self, text):
        """ Print out an error message / warning string """
        print(Fore.RED + str(text) + Fore.WHITE)

    def warn(self, text):
        """ Print out a warning message """
        print text


class VerificationError(Exception):
    """ Custom data verification exception """
    pass


class InterFits(object):
    """ InterFits: UV-data interchange class
    """

    def __init__(self, filename=None, verbosity=4):
        self.filename = filename

        self.pp = PrintLog(verbosity=verbosity)

        # Set up some basic details
        self.telescope = ""
        self.instrument = ""
        self.correlator = ""
        self.source = ""
        self.date_obs = ""
        self.obs_code = ""
        self.t_int = 0.0

        # Set up dictionaries to store data
        self.h_antenna = {}
        self.h_array_geometry = {}
        self.h_source = {}
        self.h_uv_data = {}
        self.h_frequency = {}
        self.h_common = {}
        self.h_params = {}
        self.h_flag = {}

        self.d_antenna = {}
        self.d_array_geometry = {}
        self.d_source = {}
        self.d_uv_data = {}
        self.d_frequency = {}
        self.d_flag = {}

        self.write_flags = False

        self.stokes_codes = {
            1: 'Stokes I',
            2: 'Stokes Q',
            3: 'Stokes U',
            4: 'Stokes V',
            -1: 'RR',
            -2: 'LL',
            -3: 'RL',
            -4: 'LR',
            -5: 'XX',
            -6: 'YY',
            -7: 'XY',
            -8: 'YX'
        }

        if filename:
            self.readFile(filename)

    def __repr__(self):
        to_print = ""
        to_print += "Telescope:  %s\n" % self.telescope
        to_print += "Instrument: %s\n" % self.instrument
        to_print += "Object:     %s\n" % self.source
        to_print += "Date obs:   %s\n" % self.date_obs
        return to_print

    def readFile(self, filename=None, filetype=None):
        """ Check file type, and load corresponding

        filename (str): name of file. Alternatively, if a psrdada header dictionary
                        is passed, data will be loaded from shared memory. File type
                        is inferred from extension (unless filetype arg is also passed).
        filetype (str): Defaults to none. If passed, treat file as having an explicit
                        type. Useful for when extension does not match data.
        """
        # Check what kind of file to load

        if filetype is not None:
            self._readFile(filetype)

        else:
            if filename is None:
                pass
            elif type(filename) is tuple:
                # Tuple is header_dict and numpy data array
                matched = True
                head, data = filename[0], filename[1]
                self.readDada(header_dict=head, data_arr=data)
            else:
                file_ext = os.path.splitext(filename)[1][1:]
                self._readFile(file_ext)

        self._initialize_site()

    def _readFile(self, filetype):
        """ Lookup dictionary (case statement) for file types """
        return {
            'uvfits': self.readUvfits,
            'fitsidi': self.readFitsidi,
            'FITS_1': self.readFitsidi,
            'fidi': self.readFitsidi,
            'idifits': self.readFitsidi,
            'hdf5': self.readHdf5,
            'hdf': self.readHdf5,
            'h5': self.readHdf5,
            'json': self.readJson
        }.get(filetype, self.readError)()

    def _initialize_site(self):
        """ Setup site (ephem observer)

        Uses ecef2geo function (Bowring's method), to convert
        ECEF to Lat-Long-Elev, then creates an ephem observer.
        """
        x = self.h_array_geometry["ARRAYX"]
        y = self.h_array_geometry["ARRAYY"]
        z = self.h_array_geometry["ARRAYZ"]
        lat_rad, long_rad, elev = coords.ecef2geo(x, y, z)

        self.site      = ephem.Observer()
        self.site.lon  = long_rad
        self.site.lat  = lat_rad
        self.site.elev = elev

        self.site_info(follow_verbosity=True)

    def site_info(self, follow_verbosity=False):
        """ Print information about the site """
        if not follow_verbosity:
            self.pp.pa("Telescope: %s" % self.telescope)
            self.pp.pa("Latitude:  %s" % self.site.lat)
            self.pp.pa("Longitude: %s" % self.site.long)
            self.pp.pa("Elevation: %s" % self.site.elev)

    def readError(self):
        """ Raise an error if file cannot be read """
        raise IOError("Cannot read %s" % self.filename)

    def searchKeys(self, pattern, header):
        """ Search through a header, returning a list of matching keys 
        
        pattern: str
            regular expression pattern to match
        header: pyfits header
            header of a pyfits HDU to search within
        """

        matches = []
        keys = header.keys()

        for k in keys:
            m = re.search(pattern, k)
            if m:
                matches.append(m.string)
        matches.sort()
        return matches

    def readUvfits(self):
        """ Open and read the contents of a uv-fits file 
        
        Notes
        -----
        
        Regular axes for Data matrix
        (from FITS-IDI document, assuming uv-fits adheres)
        
        Name        Mandatory   Description
        -------     ---------   -----------
        COMPLEX     yes         Real, imaginary, weight
        STOKES      yes         Stokes parameter
        FREQ        yes         Frequency (spectral channel)
        BAND        no          Band number
        RA          yes         Right ascension of phase center
        DEC         yes         Declination of phase center
        -------     ---------   -----------

        """

        self.pp.h1("Opening uvfits data")

        self.fits = pf.open(self.filename)
        #print self.fits

        # Set if source / freq tables are present
        is_source = False
        is_freq   = False

        for tbl in self.fits:
            try:
                tbl_name = tbl.header["EXTNAME"].strip()
                if tbl_name == 'AIPS AN':
                    self.antdata = tbl.data
                    self.anthead = tbl.header
                if tbl_name == "AIPS FQ":
                    is_freq = True
                    self.freqdata = tbl.data
                    self.freqhead = tbl.header
                if tbl_name == "AIPS SU":
                    is_source = True
                    self.sourcedata = tbl.data
                    self.sourcehead = tbl.header
            except KeyError:
                # UV_DATA table doesn't have EXTNAME!
                self.uvdata = tbl.data
                self.uvhead = tbl.header

        # Load basic metadata
        self.telescope = self.uvhead['TELESCOP'].strip()
        self.instrument = self.anthead['ARRNAM'].strip()
        self.source = self.uvhead['OBJECT'].strip()
        self.date_obs = self.uvhead['DATE-OBS'].strip()
        self.n_ant = self.antdata.shape[0]

        if self.pp.vlevel >= 3:
            self.fits.info()

        # Load array geometry data
        self.pp.h2("Loading array geometry")
        ag_keywords = ['ARRAYX', 'ARRAYY', 'ARRAYZ', 'ARRNAM', 'FREQ']
        ag_data = ['ANNAME', 'STABXYZ', 'NOSTA', 'MNTSTA', 'STAXOF']
        for k in ag_keywords:
            self.h_array_geometry[k] = self.anthead[k]
        for k in ag_data:
            self.d_array_geometry[k] = self.antdata[k]

        self.h_common['RDATE'] = self.anthead['RDATE']

        # Load antenna table data
        self.pp.h2("Loading antenna table")
        an_keywords = ['NOPCAL']
        an_data = ['POLTYA', 'POLAA', 'POLCALA', 'POLTYB', 'POLAB', 'POLCALB']
        for k in an_keywords:
            self.h_antenna[k] = self.anthead[k]
        for k in an_data:
            self.d_antenna[k] = self.antdata[k]

        # Load frequency table data
        # This is the first of the non-straightforward conversions
        # Need to find frequency index. It's probably CTYPE4 but safer to look it up
        self.pp.h2("Loading frequency table")
        ctypes = self.searchKeys('CTYPE\d', self.uvhead)
        c = [(ct, self.uvhead[ct], self.uvhead['NAXIS%s' % ct.lstrip('CTYPE')]) for ct in ctypes]
        freq_cid = c[[x[1] == 'FREQ' for x in c].index(True)][0].lstrip('CTYPE')

        ch_width = float(self.uvhead['CDELT%s' % freq_cid])
        ref_freq = float(self.uvhead['CRVAL%s' % freq_cid])
        ref_pixl = float(self.uvhead['CRPIX%s' % freq_cid])
        bw = ch_width * int(self.uvhead['NAXIS%s' % freq_cid])
        self.h_common['REF_FREQ'] = ref_freq
        self.h_common['CHAN_BW'] = ch_width
        self.h_common['REF_PIXL'] = ref_pixl
        self.h_common['NO_CHAN'] = int(self.uvhead['NAXIS%s' % freq_cid])

        self.d_frequency['TOTAL_BANDWIDTH'] = bw
        self.d_frequency['CH_WIDTH'] = ch_width

        # Load source table data
        self.pp.h2("Loading source table")
        if not is_source:
            self.d_source['SOURCE'] = self.source
            self.d_source['RAEPO']  = self.uvhead['OBSRA']
            self.d_source['DECEPO'] = self.uvhead['OBSDEC']
        else:
            self.d_source['SOURCE'] = self.sourcedata["SOURCE"]
            self.d_source['RAEPO']  = self.sourcedata["RAEPO"]
            self.d_source['DECEPO'] = self.sourcedata["DECEPO"]

        # Load UV-DATA
        self.pp.h2("Loading UV-data")
        self.h_uv_data['TELESCOP'] = self.telescope
        self.h_uv_data['DATE-OBS'] = self.date_obs
        self.h_uv_data['INSTRUME'] = self.instrument
        uv_datacols = ['UU', 'VV', 'WW', 'BASELINE', 'DATE']
        for k in uv_datacols:
            self.d_uv_data[k] = self.uvdata[k]

        s = self.uvdata['DATA'].shape
        # Find stokes axis type and values
        stokes_axid = 0
        for ct in ctypes:
            if self.uvhead[ct].strip() == 'STOKES':
                stokes_axid = int(ct.lstrip('CTYPE'))
                break

        stokes_axis_len = int(self.uvhead['NAXIS%i' % stokes_axid])
        stokes_code = int(self.uvhead['CRVAL%i' % stokes_axid])
        stokes_delt = int(self.uvhead['CDELT%i' % stokes_axid])
        stokes_vals = range(stokes_code, stokes_code + stokes_delt * stokes_axis_len, stokes_delt)
        self.stokes_vals = stokes_vals
        self.stokes_axis = [self.stokes_codes[i] for i in stokes_vals]

        # Should have 7 axes, likely
        if len(s) != 7:
            if len(s) == 6:
                self.pp.pp("Reshaping uv-data...")
                new_shape = (s[0], s[1], s[2], 1, s[3], s[4], s[5])
                s = new_shape
                self.d_uv_data['DATA'] = self.uvdata['DATA'].reshape(new_shape)
            else:
                self.pp.err("ERROR: Data axis not understood, incorrect length")
                self.pp.err(c)
                raise
        else:
            self.d_uv_data['DATA'] = self.uvdata['DATA']

        # Best line in the history of indexing below
        # Note the 0:2 and *2 at the end is to not include weights
        self.pp.pp("Converting DATA column to FLUX convention...")
        self.d_uv_data['FLUX'] = self.d_uv_data['DATA'][..., 0:2].astype('float32').reshape(s[0], s[4] * s[5] * 2)

        self.h_params["NSTOKES"] = len(self.stokes_vals)
        self.h_params["NBAND"] = self.d_uv_data['DATA'].shape[-4]
        self.h_params["NCHAN"] = self.d_uv_data['DATA'].shape[-3]

        self.d_uv_data.pop('DATA')

        num_rows = self.d_uv_data['FLUX'].shape[0]

        self.pp.pp("NOTE: Setting INTTIM to 1.0 (not supplied by uvfits).")
        self.pp.pp("NOTE: Setting FREQID to 1 (for FITS-IDI tables)")
        self.pp.pp("NOTE: Setting SOURCE to 1 (for FITS-IDI tables)")
        self.d_uv_data['INTTIM'] = np.ones(num_rows)
        self.d_uv_data['FREQID'] = np.ones(num_rows)
        self.d_uv_data['SOURCE'] = np.ones(num_rows)

    def readFitsidi(self, from_file=True, load_uv_data=True):
        """ Open and read the contents of a FITS-IDI file

        Notes
        -----

        Regular axes for Data matrix (from FITS-IDI document)

        Name        Mandatory   Description
        -------     ---------   -----------
        COMPLEX     yes         Real, imaginary, weight
        STOKES      yes         Stokes parameter
        FREQ        yes         Frequency (spectral channel)
        BAND        no          Band number
        RA          yes         Right ascension of phase center
        DEC         yes         Declination of phase center
        -------     ---------   -----------

        """

        if from_file:
            self.pp.h1("Opening FITS-IDI data")
            self.fits = pf.open(self.filename)
        else:
            pass  # Assumes that self.fits is already populated

        # Match tables
        opt_tbl_flag = False        # FLAG table optional
        for tbl in self.fits:
            try:
                if tbl.header['EXTNAME'] == 'ARRAY_GEOMETRY':
                    self.tbl_array_geometry = tbl
                elif tbl.header['EXTNAME'] == 'FREQUENCY':
                    self.tbl_frequency = tbl
                elif tbl.header['EXTNAME'] == 'ANTENNA':
                    self.tbl_antenna = tbl
                elif tbl.header['EXTNAME'] == 'SOURCE':
                    self.tbl_source = tbl
                elif tbl.header['EXTNAME'] == 'UV_DATA':
                    self.tbl_uv_data = tbl
                elif tbl.header['EXTNAME'] == 'FLAG':
                    self.tbl_flag = tbl
                    opt_tbl_flag = True
                else:
                    self.pp.warn("Warning: %s not recognized" % tbl.header["EXTNAME"])
            except KeyError:
                pass

        # Load basic metadata
        if load_uv_data:
            self.telescope = self.tbl_uv_data.header['TELESCOP'].strip()
            try:
                self.date_obs = self.tbl_uv_data.header['DATE-OBS'].strip()
                self.h_uv_data['DATE-OBS'] = self.date_obs
            except AttributeError:
                self.date_obs = 0.0
            except KeyError:
                self.pp.warn("Warning: DATE-OBS keyword not found in UV_DATA header")
                try:
                    self.date_obs = self.tbl_uv_data.header['RDATE'].strip()
                    self.pp.warn("Warning: Using RDATE instead of DATE-OBS (not found)")
                except KeyError:
                    self.date_obs = 0.0

        else:
            self.telescope, self.date_obs = '', 0.0
        self.instrument = self.tbl_array_geometry.header['ARRNAM'].strip()
        self.source = str(self.tbl_source.data['SOURCE'][0]).lstrip("\"[\'").rstrip("\']\"")
        self.n_ant = self.tbl_antenna.data.shape[0]

        if from_file:
            if self.pp.vlevel >= 3:
                self.fits.info()
            #self.site_info(follow_verbosity=True)

        # Load array geometry data
        self.pp.h2("Loading array geometry")
        ag_keywords = ['ARRAYX', 'ARRAYY', 'ARRAYZ', 'ARRNAM', 'FREQ']
        ag_data = ['ANNAME', 'STABXYZ', 'NOSTA', 'MNTSTA', 'STAXOF']
        for k in ag_keywords:
            self.h_array_geometry[k] = self.tbl_array_geometry.header[k]
        for k in ag_data:
            self.d_array_geometry[k] = self.tbl_array_geometry.data[k]

        # Load antenna table data
        self.pp.h2("Loading antenna table")
        an_keywords = ['NOPCAL']
        an_data = ['POLTYA', 'POLAA', 'POLCALA', 'POLTYB', 'POLAB', ' POLCALB']
        for k in an_keywords:
            self.h_antenna[k] = self.tbl_antenna.header[k]
        for k in an_data:
            try:
                self.d_antenna[k] = self.tbl_antenna.data[k]
            except KeyError:
                self.pp.warn("\tWARNING: %s key error raised." % k)

        # Load frequency table data
        # This is the first of the non-straightforward conversions
        self.pp.h2("Loading frequency table")
        frq_keywords = ['FREQID', 'BANDFREQ', 'CH_WIDTH', 'TOTAL_BANDWIDTH', 'SIDEBAND']
        for k in frq_keywords:
            try:
                self.d_frequency[k] = self.tbl_frequency.data[k]
            except KeyError:
                self.pp.warn("\tWARNING: %s key error raised." % k)

        # Load source table data
        self.pp.h2("Loading source table")
        src_keywords = ['SOURCE_ID', 'SOURCE', 'QUAL', 'CALCODE', 'FREQID', 'IFLUX',
                        'QFLUX', 'UFLUX', 'VFLUX', 'ALPHA', 'FREQOFF', 'RAEPO', 'DECEPO',
                        'EQUINOX', 'RAAPP', 'DECAPP', 'SYSVEL', 'VELTYP', 'VELDEF', 'RESTFREQ',
                        'PMRA', 'PMDEC', 'PARALLAX']
        for k in src_keywords:
            try:
                self.d_source[k] = self.tbl_source.data[k]
                if k == 'SOURCE':
                    if type(self.d_source[k]) is not str:
                        self.d_source[k] = self.d_source[k][0]
            except KeyError:
                self.pp.warn("\tWARNING: %s key error raised." % k)

        # Load common (mandatory) keywords
        self.pp.h2("Loading common keywords")
        com_keywords = ['STK_1', 'NO_BAND', 'NO_STKD', 'REF_PIXL', 'REF_FREQ', 'CHAN_BW', 'NO_CHAN', 'RDATE']
        for k in com_keywords:
            try:
                self.h_common[k] = self.tbl_frequency.header[k]
            except KeyError:
                self.pp.warn("\tWARNING: %s key error raised." % k)

        # Also fill in the parameter header dictionary (needed for XML generation).
        self.h_params["NSTOKES"] = self.h_common["NO_STKD"]
        self.h_params["NBAND"] = self.h_common["NO_BAND"]
        self.h_params["NCHAN"] = self.h_common["NO_CHAN"]

        # Load UV-DATA
        if load_uv_data:
            self.pp.h2("Loading UV-data")
            uv_keywords = ['TELESCOP']
            for k in uv_keywords:
                try:
                    self.h_uv_data[k] = self.tbl_uv_data.header[k]
                except KeyError:
                    self.pp.warn("\tWARNING: %s key error raised." % k)

            uv_datacols = ['UU', 'VV', 'WW', 'BASELINE', 'DATE', 'FLUX', 'INTTIM', 'FREQID', 'SOURCE']
            for k in uv_datacols:
                self.d_uv_data[k] = self.tbl_uv_data.data[k]

            try:
                self.d_uv_data["TIME"] = self.tbl_uv_data.data["TIME"]
            except KeyError:
                self.pp.warn("\tWARNING: TIME column does not exist.")
                raise

            self.d_uv_data["FLUX"] = self.d_uv_data["FLUX"].astype('float32')

            # Find stokes axis type and values
            stokes_axid = 0
            ctypes = self.searchKeys('CTYPE\d', self.tbl_uv_data.header)
            for ct in ctypes:
                if self.tbl_uv_data.header[ct].strip() == 'STOKES':
                    stokes_axid = int(ct.lstrip('CTYPE'))
                    break

            stokes_axis_len = int(self.tbl_uv_data.header['MAXIS%i' % stokes_axid])
            stokes_code = int(self.tbl_uv_data.header['CRVAL%i' % stokes_axid])
            stokes_delt = int(self.tbl_uv_data.header['CDELT%i' % stokes_axid])
            stokes_vals = range(stokes_code, stokes_code + stokes_delt * stokes_axis_len, stokes_delt)
            self.stokes_vals = stokes_vals
            self.stokes_axis = [self.stokes_codes[i] for i in stokes_vals]

        if opt_tbl_flag:
            self.write_flags = True
            self.pp.h2("Loading FLAG table")
            flag_keywords = ['SOURCE_ID', 'ARRAY', 'ANTS', 'FREQID', 'BANDS', 'CHANS', 'PFLAGS', 'REASON', 'SEVERITY']
            for k in flag_keywords:
                try:
                    self.d_frequency[k] = self.tbl_frequency.data[k]
                except KeyError:
                    self.pp.warn("\tWARNING: %s key error raised." % k)

    def readJson(self, filename=None):
        """ Read JSON data into InterFits dictionaries.

        Notes
        -----
        This expects the filename to be any *.json file in a directory of .json files
        """

        if filename:
            filepath = os.path.split(filename)[0]
        else:
            filepath = os.path.split(self.filename)[0]

        self.pp.h1("Opening JSON data")

        if not os.path.exists(filepath):
            raise IOError("Cannot read directory %s" % filepath)

        try:
            self.pp.h2("Loading common keywords")
            self.h_common = load_json(os.path.join(filepath, 'h_common.json'))
            self.h_params = load_json(os.path.join(filepath, 'h_params.json'))
        except IOError:
            self.pp.warn("Warning: Could not load common keywords")
        try:
            self.pp.h2("Loading antenna table")
            self.h_antenna = load_json(os.path.join(filepath, 'h_antenna.json'))
            self.d_antenna = load_json(os.path.join(filepath, 'd_antenna.json'))
        except IOError:
            self.pp.warn("Warning: Could not load antenna table")
        try:
            self.pp.h2("Loading frequency table")
            self.h_frequency = load_json(os.path.join(filepath, 'h_frequency.json'))
            self.d_frequency = load_json(os.path.join(filepath, 'd_frequency.json'))
        except IOError:
            self.pp.warn("Warning: Could not load frequency table")
        try:
            self.pp.h2("Loading array geometry table")
            self.h_array_geometry = load_json(os.path.join(filepath, 'h_array_geometry.json'))
            self.d_array_geometry = load_json(os.path.join(filepath, 'd_array_geometry.json'))
        except IOError:
            self.pp.warn("Warning: Could not load array geometry table")
        try:
            self.pp.h2("Loading source table")
            self.h_source = load_json(os.path.join(filepath, 'h_source.json'))
            self.d_source = load_json(os.path.join(filepath, 'd_source.json'))
        except IOError:
            self.pp.warn("Warning: Could not load frequency table")
        try:
            self.pp.h2("Loading UV_DATA (header)")
            self.h_uv_data = load_json(os.path.join(filepath, 'h_uv_data.json'))
        except IOError:
            self.pp.warn("Warning: Could not load UV_DATA table header")
        try:
            self.pp.h2("Loading FLAG")
            self.h_flag = load_json(os.path.join(filepath, 'h_flag.json'))
            self.d_flag = load_json(os.path.join(filepath, 'd_flag.json'))
        except IOError:
            self.pp.warn("Warning: Could not load UV_DATA table header")

        try:
            self.telescope = self.h_uv_data['TELESCOP']
        except KeyError:
            self.pp.warn("Warning: Could not load TELESCOP from UV_DATA header")
        try:
            self.date_obs = self.h_uv_data['DATE-OBS']
        except KeyError:
            self.pp.warn("Warning: Could not load DATE-OBS from UV_DATA header")
        try:
            s = self.d_source['SOURCE']
            if type(s) is list:
                s = str(s[0])
            self.source = s
        except KeyError:
            self.pp.warn("Warning: Could not load SOURCE from UV_DATA header")
        try:
            self.instrument = self.h_array_geometry['ARRNAM']
        except KeyError:
            self.pp.warn("Warning: Could not load ARRNAM from UV_DATA header")

    def readHdf5(self):
        """ Read data from HDF5 file. """
        self.pp.h1("Reading HDF5 file %s" % self.filename)
        self.hdf = h5py.File(self.filename, "r")
        try:
            ifds = [self.h_antenna, self.h_source, self.h_array_geometry, self.h_frequency, self.h_uv_data,
                    self.d_antenna, self.d_source, self.d_array_geometry, self.d_frequency, self.d_uv_data,
                    self.h_common, self.h_params, self.h_flag, self.d_flag]

            ifd_names = ["h_antenna", "h_source", "h_array_geometry", "h_frequency", "h_uv_data",
                         "d_antenna", "d_source", "d_array_geometry", "d_frequency", "d_uv_data",
                         "h_common", "h_params", "h_flag", "d_flag"]

            for ii in range(len(ifds)):
                ifd = ifds[ii]
                ifd_name = ifd_names[ii]

                if ifd_name in self.hdf.keys():
                    self.pp.h3("Reading %s" % ifd_name)
                    h5d = self.hdf[ifd_name]
                    self.pp.pp(h5d)
                    for key in h5d.keys():
                        if ifd_name.startswith("h_"):
                            # Convert back to native types from numpy types
                            int_types = type(np.int64(1)), type(np.int32(1))
                            float_types = type(np.float64(1)), type(np.float32(1))
                            str_types = type(np.string_("hi")), type(np.unicode_("hi"))

                            if type(h5d[key][0]) in int_types:
                                #print "INT"
                                ifd[key] = int(h5d[key][0])
                                #print type(ifd[key])
                            if type(h5d[key][0]) in float_types:
                                #print "FLOAT"
                                ifd[key] = float(h5d[key][0])
                                #print type(ifd[key])
                            if type(h5d[key][0]) in str_types:
                                #print "STRING"
                                ifd[key] = str(h5d[key][0])
                                #print type(ifd[key])
                        else:
                            ifd[key] = h5d[key][:]

            # TODO: Unhardcode this
            self.stokes_axis = ['XX', 'YY', 'XY', 'YX']
            self.stokes_vals = [-5, -6, -7, -8]

            self.n_ant = len(self.d_array_geometry['NOSTA'])

        except ValueError:
            self.hdf.close()
            self.pp.warn(key)
            self.pp.warn(h5d[key])
            raise

        self.date_obs = self.h_uv_data["DATE-OBS"]
        try:
            self.telescope = self.h_uv_data["TELESCOP"]
            self.instrument = self.h_array_geometry["ARRNAM"]
        except KeyError:
            self.telescope = 'UNKNOWN'
            self.instrument = 'UNKNOWN'
        self.source = self.d_source["SOURCE"][0]

    def setXml(self, table, keyword, value):
        """ Find a header parameter and replace it """
        try:
            self.xmlroot = self.xmlData.getroot()
            if type(value) == str or type(value) == type(u'Unicode'):
                # Make sure strings are in single quotes
                value = "'" + value.strip("'") + "'"
            self.xmlroot.find(table).find(keyword).text = str(value)
        except:
            self.pp.err("Error: Something went wrong with XML parsing")
            self.pp.err("%s, %s, %s" % (table, keyword, value))

    def s2arr(self, val):
        """ Put a single value into a numpy array """
        return np.array([val])

    def setDefaults(self, n_uv_rows):
        """ FIll headers and data with default data """

        zero_vec = np.zeros(n_uv_rows).astype('float32')
        ones_vec = np.ones(n_uv_rows).astype('float32')
        self.d_uv_data["DATE"] = zero_vec
        self.d_uv_data["UU"] = zero_vec
        self.d_uv_data["VV"] = zero_vec
        self.d_uv_data["WW"] = zero_vec
        self.d_uv_data["FREQID"] = ones_vec
        self.d_uv_data["INTTIM"] = ones_vec
        self.d_uv_data["SOURCE"] = ones_vec

        self.stokes_axis = ['XX', 'YY', 'XY', 'YX']
        self.stokes_vals = [-5, -6, -7, -8]

        self.d_array_geometry["ANNAME"] = \
            np.array(["Stand%d" % (i + 1) for i in range(len(self.d_array_geometry["ANNAME"]))])
        self.d_array_geometry["NOSTA"] = \
            np.array([i + 1 for i in range(len(self.d_array_geometry["NOSTA"]))])

        self.d_frequency["FREQID"] = self.s2arr(1)
        self.d_frequency["BANDFREQ"] = self.s2arr(0)

        #self.d_source["EQUINOX"]   = self.s2arr('J2000')
        self.d_source["SOURCE"] = self.s2arr('ZENITH')
        self.d_source["SOURCE_ID"] = self.s2arr(1)
        #self.d_source["VELDEF"]    = self.s2arr("RADIO")
        #self.d_source["VELTYP"]    = self.s2arr("GEOCENTR")

        self.h_uv_data["DATE-OBS"] = '2013-01-01T00:00:00.0'

    def generateFitsidiXml(self, xmlbase=None, filename_out=None):
        """ Generate XML file that encodes fitsidi structure 
        
        xmlbase: str
            name of basic input xml file
        filename_out: str
            name of output file
        """

        if xmlbase is None:
            dirname, filename = os.path.split(os.path.abspath(__file__))
            xmlbase = os.path.join(dirname, 'config/config.xml')

        self.xmlData = etree.parse(xmlbase)

        # Look in the config file or fits-idi convention for
        # info about what these should be and refer to
        self.setXml("PARAMETERS", "NSTOKES", self.h_params["NSTOKES"])
        self.setXml("PARAMETERS", "NBAND", self.h_params["NBAND"])
        self.setXml("PARAMETERS", "NCHAN", self.h_params["NCHAN"])

        # Common headers - required for each table
        self.setXml("COMMON", "OBSCODE", self.obs_code)
        self.setXml("COMMON", "STK_1", self.stokes_vals[0])
        self.setXml("COMMON", "REF_FREQ", self.h_common['REF_FREQ'])
        self.setXml("COMMON", "CHAN_BW", self.h_common['CHAN_BW'])
        self.setXml("COMMON", "REF_PIXL", self.h_common["REF_PIXL"])
        self.setXml("COMMON", "RDATE", self.h_common["RDATE"])

        self.setXml("ANTENNA", "NOPCAL", self.h_antenna["NOPCAL"])

        # Support tables
        self.setXml("PRIMARY", "CORRELAT", self.correlator)
        self.setXml("ARRAY_GEOMETRY", "ARRNAM", self.h_array_geometry["ARRNAM"])

        self.setXml("ARRAY_GEOMETRY", "FREQ", self.h_array_geometry['FREQ'])
        self.setXml("ARRAY_GEOMETRY", "ARRAYX", self.h_array_geometry["ARRAYX"])
        self.setXml("ARRAY_GEOMETRY", "ARRAYY", self.h_array_geometry["ARRAYY"])
        self.setXml("ARRAY_GEOMETRY", "ARRAYZ", self.h_array_geometry["ARRAYZ"])

        # UV-DATA
        self.setXml("UV_DATA", "DATE-OBS", self.date_obs)
        self.setXml("UV_DATA", "TELESCOP", self.telescope)
        stokes_delt = self.stokes_vals[1] - self.stokes_vals[0]
        self.setXml("UV_DATA", "CDELT2", stokes_delt)
        self.setXml("UV_DATA", "CRVAL2", self.stokes_vals[0])
        if type(self.d_frequency['CH_WIDTH']) is np.ndarray:
            self.setXml("UV_DATA", "CDELT3", self.d_frequency['CH_WIDTH'][0])
        else:
            self.setXml("UV_DATA", "CDELT3", self.d_frequency['CH_WIDTH'])
        self.setXml("UV_DATA", "CRVAL3", self.h_common['REF_FREQ'])

        if filename_out:
            if os.path.isfile(filename_out):
                os.remove(filename_out)
            self.pp.pp("Writing to %s" % filename_out)
            with open(filename_out, 'w') as f:
                f.write(etree.tostring(self.xmlData))

    def exportJson(self, dirname_out, dump_uv_data=False):
        """ Export data as a directory of JSON files.

        dirname_out: str
            name of directory to output files into.
        dump_uv_data: bool
            Dump UV DATA into a dictionary? Defaults to false as this is very large.
        """

        if not os.path.exists(dirname_out):
            os.mkdir(dirname_out)

        self.pp.h1("Creating JSON-Numpy dictionaries in %s" % dirname_out)
        dump_json(self.h_antenna, os.path.join(dirname_out, 'h_antenna.json'))
        dump_json(self.h_array_geometry, os.path.join(dirname_out, 'h_array_geometry.json'))
        #dump_json(self.h_common, os.path.join(dirname_out, 'h_common.json'))
        dump_json(self.h_frequency, os.path.join(dirname_out, 'h_frequency.json'))
        #dump_json(self.h_params, os.path.join(dirname_out, 'h_params.json'))
        dump_json(self.h_uv_data, os.path.join(dirname_out, 'h_uv_data.json'))
        dump_json(self.h_source, os.path.join(dirname_out, 'h_source.json'))
        dump_json(self.d_antenna, os.path.join(dirname_out, 'd_antenna.json'))
        dump_json(self.d_array_geometry, os.path.join(dirname_out, 'd_array_geometry.json'))
        dump_json(self.d_frequency, os.path.join(dirname_out, 'd_frequency.json'))
        dump_json(self.d_source, os.path.join(dirname_out, 'd_source.json'))

        if self.write_flags:
            dump_json(self.d_flag, os.path.join(dirname_out, 'd_flag.json'))
            dump_json(self.h_flag, os.path.join(dirname_out, 'h_flag.json'))

        if dump_uv_data:
            dump_json(self.d_uv_data, os.path.join(dirname_out, 'd_uv_data.json'))

    def exportHdf5(self, filename_out, compression=None, scaleoffset=None):
        """ Export data as HDF5 file

        filename_out: str
            name of output files into.
        compression: str
            Apply compression to data? Defaults to none, can be one of:
            None, gzip (slow), lzf, szip (if supported), scaleoffset
        scaleoffset: int
            Scale offset, used only if using scale offset filter-based compression.
            This is a lossy compression scheme, high integers = higher loss

        """
        self.pp.h1("Exporting to %s" % filename_out)
        self.hdf = h5py.File(filename_out, "w")

        self.h_uv_data["TELESCOP"]      = self.telescope
        self.h_array_geometry["ARRNAM"] = self.instrument

        ifds = [self.h_antenna, self.h_source, self.h_array_geometry, self.h_frequency, self.h_uv_data,
                self.d_antenna, self.d_source, self.d_array_geometry, self.d_frequency, self.d_uv_data,
                self.h_common, self.h_params]

        ifd_names = ["h_antenna", "h_source", "h_array_geometry", "h_frequency", "h_uv_data",
                     "d_antenna", "d_source", "d_array_geometry", "d_frequency", "d_uv_data",
                     "h_common", "h_params"]

        if self.write_flags:
            ifds.append(self.h_flag)
            ifds.append(self.d_flag)
            ifd_names.append("h_flag")
            ifd_names.append("d_flag")

        for ii in range(len(ifds)):
            ifd = ifds[ii]
            ifd_name = ifd_names[ii]

            self.pp.h2("Creating %s" % ifd_name)
            hgroup = self.hdf.create_group(ifd_name)
            for key in ifd:

                if type(ifd[key]) is list:
                    if type(ifd[key][0]) is unicode:
                        ifd[key] = np.array(ifd[key], dtype='str')

                if type(ifd[key]) in (str, int, float, unicode):
                    # Key must be unicode!
                    key = str(key)

                    if type(ifd[key]) is unicode:
                        ifd[key] = str(ifd[key])
                    if compression is not None:
                        try:
                            hgroup.create_dataset(key, data=[ifd[key]],
                                                  compression=compression, scaleoffset=scaleoffset)
                        except TypeError:
                            hgroup.create_dataset(key, data=[ifd[key]])
                    else:
                        hgroup.create_dataset(key, data=[ifd[key]])
                else:
                    if compression is not None:
                        try:
                            hgroup.create_dataset(key, data=[ifd[key]],
                                                  compression=compression, scaleoffset=scaleoffset)
                        except TypeError:
                            hgroup.create_dataset(key, data=ifd[key])
                    else:
                        hgroup.create_dataset(key, data=ifd[key])
        self.hdf.close()

    def exportFitsidi(self, filename_out, config_xml=None, verbose=False):
        """ Export data as FITS IDI 
        
        filename_out (str): output filename
        config_xml (str): path to config file
        
        """
        self.pp.h1("Exporting to FITS-IDI")
        self.pp.h2('Generating FITS-IDI XML schema')
        if config_xml is None:
            dirname, this_file = os.path.split(os.path.abspath(__file__))
            config_xml = os.path.join(dirname, 'config/config.xml')
        xmlfile = filename_out.replace(".fitsidi", "").replace(".fits", "") + ".xml"
        self.generateFitsidiXml(config_xml, xmlfile)
        config_xml = xmlfile

        self.pp.h2('Creating Primary HDU')
        hdu_primary = make_primary(config=config_xml)
        self.pp.debug(hdu_primary.header.ascard)

        self.pp.h2('Creating ARRAY_GEOMETRY')
        tbl_array_geometry = make_array_geometry(config=config_xml, num_rows=self.n_ant)
        self.pp.debug(tbl_array_geometry.header.ascard)

        self.pp.h2('Creating ANTENNA')
        tbl_antenna = make_antenna(config=config_xml, num_rows=self.n_ant)
        self.pp.debug(tbl_antenna.header.ascard)

        self.pp.h2('Creating FREQUENCY')
        tbl_frequency = make_frequency(config=config_xml, num_rows=1)
        self.pp.debug(tbl_frequency.header.ascard)

        self.pp.h2('Creating SOURCE')
        tbl_source = make_source(config=config_xml, num_rows=1)
        self.pp.debug(tbl_source.header.ascard)

        self.pp.h2('Creating UV_DATA')
        num_rows = self.d_uv_data['FLUX'].shape[0]
        uvd = self.d_uv_data
        if type(uvd['BASELINE']) is list:
            uvd['BASELINE'] = np.array(uvd['BASELINE'])

        try:
            jtime = uvd['TIME']
        except KeyError:
            self.pp.warn("\t Warning: TIME column does not exist.")
            jtime = None

        # TODO: Fix time and date to julian date
        tbl_uv_data = make_uv_data(config=config_xml, num_rows=num_rows,
                                   uu_data=uvd['UU'], vv_data=uvd['VV'], ww_data=uvd['WW'],
                                   date_data=uvd['DATE'], time_data=jtime,
                                   baseline_data=uvd['BASELINE'].astype('int32'),
                                   source_data=uvd['SOURCE'].astype('int32'), freqid_data=uvd['FREQID'].astype('int32'),
                                   inttim_data=uvd['INTTIM'],
                                   weights_data=None, flux_data=uvd['FLUX'], weights_col=False)

        self.pp.debug(tbl_uv_data.header.ascard)

        if self.write_flags:
            n_rows_flag = 0
            if self.d_flag.get("SOURCE_ID"):
                n_rows_flag = len(self.d_flag["SOURCE_ID"])
            if n_rows_flag > 0:
                self.pp.h2('Creating FLAG')
                tbl_flag = make_flag(config=config_xml, num_rows=n_rows_flag)
                self.pp.debug(tbl_flag.header.ascardlist())

        self.pp.h2('Filling in data')
        self.pp.h3("ARRAY_GEOMETRY")
        for i in range(self.n_ant):
            for k in ['ANNAME', 'STABXYZ', 'NOSTA', 'MNTSTA', 'STAXOF']:
                tbl_array_geometry.data[k][i] = self.d_array_geometry[k][i]

        self.pp.h3("ANTENNA")
        for i in range(self.n_ant):
            # TODO: 'POLCALB' and POLCALA               
            for k in ['POLTYA', 'POLAA', 'POLTYB', 'POLAB']:
                try:
                    tbl_antenna.data['ANNAME'][i] = self.d_array_geometry['ANNAME'][i]
                    tbl_antenna.data['ANTENNA_NO'][i] = i + 1
                    tbl_antenna.data['ARRAY'][i] = 1
                    tbl_antenna.data['FREQID'][i] = 1
                    tbl_antenna.data['NO_LEVELS'][i] = 255
                    tbl_antenna.data[k][i] = self.d_antenna[k][i]
                except:
                    self.pp.warn("Warning: keyword error: %s" % k)

        self.pp.h3("FREQUENCY")
        tbl_frequency.data["FREQID"][0] = 1
        tbl_frequency.data['BANDFREQ'][0] = 0
        tbl_frequency.data['CH_WIDTH'][0] = self.d_frequency['CH_WIDTH']
        tbl_frequency.data['TOTAL_BANDWIDTH'][0] = self.d_frequency['TOTAL_BANDWIDTH']

        self.pp.h3("SOURCE")
        if type(self.d_source['SOURCE']) == str:
            for k in ['SOURCE', 'RAEPO', 'DECEPO']:
                try:
                    tbl_source.data['SOURCE_ID'][0] = 1
                    tbl_source.data['EQUINOX'][0] = 'J2000'
                    tbl_source.data[k][0] = self.d_source[k]
                except:
                    self.pp.warn("Warning: keyword error: %s" % k)
                    raise
        else:
            n_rows = len(self.d_source['SOURCE'])
            #print self.d_source
            for i in range(n_rows):
                tbl_source.data['SOURCE_ID'][i] = i + 1
                tbl_source.data['EQUINOX'][i] = 'J2000'
                for k in ['SOURCE', 'RAEPO', 'DECEPO']:
                    try:
                        tbl_source.data[k][i] = self.d_source[k][i]
                    except:
                        self.pp.warn("Warning: keyword error: %s" % k)

        self.pp.h3("UV_DATA")
        self.pp.h3("(Pre-filled)")
        # NOTE: This is now superfluous, thanks to the make_uv_data call above
        #for i in range(self.d_uv_data['DATA'].shape[0]):
        #    LinePrint("Row %i of %i"%(i+1, self.d_uv_data['DATA'].shape[0]))
        #    for k in ['UU','VV','WW','BASELINE','DATE']:
        #        try:
        #            tbl_uv_data.data[k][i] = self.d_uv_data[k][i]
        #            #tbl_uv_data.data['FLUX'][i] = self.d_uv_data['FLUX'][i]
        #        except:
        #            raise

        if self.write_flags:
            if n_rows_flag > 0:
                self.pp.h3("FLAG")
                for i in range(n_rows_flag):
                    flag_keywords = ['SOURCE_ID', 'ARRAY', 'ANTS', 'FREQID', 'BANDS', 'CHANS', 'PFLAGS', 'REASON',
                                     'SEVERITY']
                    for k in flag_keywords:
                        try:
                            tbl_flag.data[k][i] = self.d_flag[k][i]
                        except:
                            self.pp.warn("Warning: keyword error: %s" % k)

        # Add history and comments to header
        hdu_primary.header.add_comment("FITS-IDI: FITS Interferometry Data Interchange Convention")
        hdu_primary.header.add_comment("defined at http://fits.gsfc.nasa.gov/registry/fitsidi.html")
        now = datetime.now()
        datestr = now.strftime("Interfits: File created %Y-%m-%dT%H:%M:%S")
        hdu_primary.header.add_history(datestr)
        self.pp.h1('Creating HDU list')
        hdus = [hdu_primary,
                tbl_array_geometry,
                tbl_frequency,
                tbl_antenna,
                tbl_source,
                tbl_uv_data]

        if self.write_flags and n_rows_flag > 0:
            hdus.append(tbl_flag)
        hdulist = pf.HDUList(hdus)
        if self.pp.vlevel >= 3:
            hdulist.info()

        self.pp.pp('\nVerifying integrity...')
        hdulist.verify()

        if os.path.isfile(filename_out):
            self.pp.pp('Removing existing file %s...' % filename_out)
            os.remove(filename_out)
        self.pp.pp('Writing to file %s...' % filename_out)
        hdulist.writeto(filename_out)

    def verify_baseline_order(self):
        """ Check baseline IDs are in order """

        self.pp.pp("Verification: Checking uv_data baseline order...")

        bls = [int(b) for b in self.d_uv_data['BASELINE']]

        # Generate lower an upper triangular matrices
        bl_lower, bl_upper = [], []
        for i in range(1, self.n_ant + 1):
            for j in range(1, self.n_ant + 1):
                if j >= i:
                    if j >= 256:
                        bl_lower.append(2048 * i + j + 65536)
                    else:
                        bl_lower.append(256 * i + j)
                elif j <= i:
                    if j >= 256:
                        bl_upper.append(2048 * j + i + 65536)
                    else:
                        bl_upper.append(256 * j + i)

        # Check every baseline is right, over all dumps
        n_bls = self.n_ant * (self.n_ant - 1) / 2 + self.n_ant
        n_dumps = len(bls) / n_bls

        lower_t, upper_t = False, False
        for i in range(n_dumps):
            if bls[i * n_bls:(i + 1) * n_bls] == bl_lower:
                #print "LOWER"
                upper_t = False
                lower_t = True
            if bls[i * n_bls:(i + 1) * n_bls] == bl_upper:
                #print "UPPER"
                lower_t = False
                upper_t = True

            if not lower_t and not upper_t:
                raise VerificationError("Baseline order neither upper or lower triangular.")

        if lower_t:
            self.pp.pp("Verification: OK. Baselines in lower triangular order.")
        if upper_t:
            self.pp.pp("Verification: OK. Baselines in upper triangular order.")

        return True

    def verify_uv_table(self):
        """ Basic diagnostics on UV_DATA table """
        self.pp.pp("Verification: checking UV_DATA for null entries")

        freq_ids   = self.d_uv_data["FREQID"]
        source_ids = self.d_uv_data["SOURCE"]
        baselines  = self.d_uv_data["BASELINE"]

        if 0 in freq_ids:
            raise VerificationError("FREQID in UV_DATA references non-existent FREQ with ID 0")
        if 0 in source_ids:
            raise VerificationError("SOURCE in UV_DATA references non-existent SOURCE with ID 0")
        if 0 in baselines:
            raise VerificationError("BASELINE in UV_DATA references non-existent BASELINE with ID 0")

        self.pp.pp("Verification: OK. UV_DATA does not contain null (zero) entries in required fields")
        return True

    def verify_frequency_axis(self):
        """ Verify frequency values are sensical """
        try:
            f = self.formatFreqs()
            self.pp.pp("Verification: OK. Frequency axis spans valid range.")
        except ValueError:
            raise VerificationError("Frequency values are fubarred")

    def verify(self):
        """ Run a series of diagnostics to test data validity """
        self.pp.h1("Data verification")
        self.verify_uv_table()
        self.verify_baseline_order()
        self.verify_frequency_axis()

    def formatStokes(self):
        """ Return data as complex stokes vector """
        data = self.d_uv_data["FLUX"]

        #xx_data = data[:, 0::8] + 1j * data[:, 1::8]
        #yy_data = data[:, 2::8] + 1j * data[:, 3::8]
        #xy_data = data[:, 4::8] + 1j * data[:, 5::8]
        #yx_data = data[:, 6::8] + 1j * data[:, 7::8]

        if self.d_uv_data["FLUX"].dtype == 'float32':
            data = data.view('complex64')
        elif self.d_uv_data["FLUX"].dtype == 'float64':
            data = data.view('complex128')

        if self.h_params["NSTOKES"] == 1:
            return data
        elif self.h_params["NSTOKES"] == 2:
            return np.array((data[:, ::2], data[:, 1::2]))
        elif self.h_params["NSTOKES"] == 4:
            data2 = np.array((data[:, ::4], data[:, 1::4], data[:, 2::4], data[:, 3::4]))
            return data2
        else:
            raise ValueError("NSTOKES in h_params is not valid!")

    def formatFreqs(self):
        """ Convert FITS keywords to frequency array """
        ref_delt = self.h_common["CHAN_BW"]
        ref_pix  = self.h_common["REF_PIXL"]
        ref_val  = self.h_common["REF_FREQ"]
        num_pix  = self.h_common["NO_CHAN"]
        freqs    = np.arange(0, num_pix, 1) * ref_delt + (ref_val - ref_pix * ref_delt)
        try:
            assert np.min(freqs) >= 0
            assert ref_delt > 0
        except Exception, e:
            self.pp.err("CHAN_BW: %s\n REF_PIXL: %s\n REF_FREQ: %s\n NO_CHAN %s" \
                        % (ref_delt, ref_pix, ref_val, num_pix))
            self.pp.pp("Error: %s" % str(e))
            raise ValueError("Frequency values are fubarred.")

        return freqs

    def get_antenna_id(self, bl_id):
        """ Convert baseline ID into an antenna pair.

        Uses MIRIAD convention for antennas > 256
        Returns a tuple of antenna IDs.
        """
        if bl_id > 65536:
            ant1 = (bl_id - 65536) / 2048
            ant2 = (bl_id - 65536) % 2048
        else:
            ant1 = bl_id / 256
            ant2 = bl_id % 256
        return ant1, ant2

    def get_baseline_id(self, ant1, ant2):
        """ Convert antenna pair into baseline ID """
        if ant1 > 255 or ant2 > 255:
            bl_id = ant1 * 2048 + ant2 + 65536
        else:
            bl_id = ant1 * 256 + ant2
        return bl_id

    def search_baselines(self, ref_ant, autocorrs=True):
        """ Retrieve baseline ids that contain a given antenna

        ref_ant (int):  Antenna of interest
        triangle (str): Defaults to upper matrix order (lower not supported yet)

        returns: a list of all matching baseline IDs

        TODO: Add different triangle support
        """
        try:
            assert ref_ant <= self.n_ant
        except AssertionError:
            raise RuntimeError("Ref ant is larger than number of antennas: %i > %i" % (ref_ant, self.n_ant))

        bls, ant_arr = [], []
        for ii in range(1, self.n_ant + 1):
            for jj in range(1, self.n_ant + 1):
                if jj >= ii:
                    if autocorrs is False and ii == jj:
                        pass
                    else:
                        ant_arr.append((ii, jj))
                        if ii > 255 or jj > 255:
                            bl_id = ii * 2048 + jj + 65536
                        else:
                            bl_id = 256 * ii + jj
                        if ii == ref_ant or jj == ref_ant:
                            bls.append(bl_id)
        return bls

    def extract_integrations(self, start=None, stop=None):
        """ Extract a subset of integrations from the file

        Parameters
        ----------
        start: int
            Start of slice point. Integrations are numbered from zero (0).
        stop: int
            End of slice. For example, start=0, stop=1 will return the first integration.
        """

        bls = self.d_uv_data["BASELINE"]
        n_bls = self.n_ant * (self.n_ant - 1) / 2 + self.n_ant
        n_dumps = len(bls) / n_bls

        if start is None:
            start = 0
        if stop is None:
            stop = n_dumps

        try:
            assert stop <= n_dumps
            assert start <= stop
        except AssertionError:
            raise ValueError("Integration start and stop points invalid. (%s, %s)" % (start, stop))

        for k in self.d_uv_data.keys():
            self.d_uv_data[k] = self.d_uv_data[k][start * n_bls:stop * n_bls]

    def extract_antenna(self, antenna_id, timestamps=False):
        """ Extract autocorrelation of a give antenna

        Parameters
        ----------
        antenna_id: int
            ID of antenna to extract
        timestamps: bool
            Default False. Returns (timestamps, data) tuple if true,
            else returns only data

        Returns
        -------
        Returns array with dimensions (n_stokes, n_int, n_channel)
        if timestamps arg is set to True, returns (timestamps, data)
        """

        bls   = self.d_uv_data["BASELINE"]
        bl_id = antenna_id * 256 + antenna_id

        try:
            stokes = self.stokes
        except AttributeError:
            self.stokes = self.formatStokes()
            stokes = self.stokes

        data = stokes[:, bls == bl_id]
        if timestamps is False:
            return data
        else:
            ts = self.d_uv_data["TIME"]
            ts = ts[bls == bl_id]
            return ts, data

