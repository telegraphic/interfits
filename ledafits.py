#! /usr/bin/env python
# encoding: utf-8
"""
ledafits.py
============

Extension of interfits.py with LEDA-OVRO specific methods, such as ability to read
L-files and DADA files.
"""

import time
import calendar
from datetime import datetime, timedelta
import ephem

from interfits import *
from lib import dada, uvw
import leda_config

# Load globals from config file
OFFSET_DELTA, INT_TIME, N_INT = leda_config.OFFSET_DELTA, leda_config.INT_TIME, leda_config.N_INT_PER_FILE 

class HeaderDataUnit(object):
    """ Very basic object with header and data units """
    def __init__(self, header, data):
        self.header = header
        self.data   = data

class LedaFits(InterFits):
    """ LEDA extension of InterFits class 
    
    This adds ability to read LA, LC and dada files, and adds helper functions
    for computing UVW coordinates, generating timestamps, and computing zenith RA/DEC.
    """

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

    def _readFile(self, filetype):
        """ Lookup dictionary (case statement) for file types """
        return {
                'uvfits': self.readUvfits,
                'fitsidi': self.readFitsidi,
                'fidi': self.readFitsidi,
                'idifits': self.readFitsidi,
                'hdf5': self.readHdf5,
                'hdf': self.readHdf5,
                'h5': self.readHdf5,
                'lfile': self.readLfile,
                'LA': self.readLfile,
                'LC': self.readLfile,
                'json': self.readJson,
                'dada': self.readDada
        }.get(filetype, self.readError)()

    def _readLfile(self, n_ant=256, n_pol=2, n_chans=109, n_stk=4):
        """ Main L-File reading subroutine.
        Opens L-files and forms a visibility matrix.
        See readLfile for main routine """

        n_antpol = n_ant * n_pol
        n_blcc = n_antpol * (n_antpol - 1) / 2

        filename = self.filename.rstrip('.LA').rstrip('.LC')

        # Autocorrs
        #h2("Opening autocorrs (.LA)")
        lfa = np.fromfile(filename + '.LA', dtype='float32')
        lfa = lfa.reshape([len(lfa) / n_chans / n_antpol, n_antpol, n_chans, 1])
        lfa = np.concatenate((lfa, np.zeros_like(lfa)), axis=3)
        print lfa.shape

        # Cross-corrs
        #h2("Opening cross-corrs (.LC)")
        lfc = np.fromfile(filename + '.LC', dtype='float32')
        lfc = lfc.reshape([len(lfc) / n_chans / n_blcc / 2, n_blcc, n_chans, 2])
        print lfc.shape

        #h2("Forming visibility matrix")
        # Create a visibility matrix, and use indexing to populate upper triangle
        n_dumps = lfa.shape[0]
        vis = np.zeros([n_dumps, n_antpol, n_antpol, n_chans, 2], dtype='float32')
        iup = np.triu_indices(n_antpol, 1)
        idiag = (np.arange(0, n_antpol), np.arange(0, n_antpol))

        for ii in range(0, vis.shape[0]):
            LinePrint("%i of %i" % (ii, vis.shape[0]))
            vis[ii][iup] = lfc[ii]
            vis[ii][idiag] = lfa[ii]
        print "\n", vis.shape

        return vis

    def readLfile(self, n_ant=256, n_pol=2, n_chans=109, n_stk=4, config_xml=None):
        """ Read a LEDA L-file 
        
        filename: str
            name of L-file
        n_ant: int
            Number of antennas. Defaults to 256
        n_pol: int
            Number of polarizations. Defaults to 2 (dual-pol)
        n_chans: int
            Number of channels in file. Defaults to 109 (LEDA-512 default)
        n_stk: int
            Number of stokes parameters in file. Defaults to 4
        config_xml: str
            Filename of XML schema file. If None, will default to [filename].xml

        Notes
        -----

        .LA and .LC are binary data streams.

        .LA files store autocorrelations in the following way:
        t0 |Ant1 109 chans XX | Ant1 109 chans YY| Ant2 ... | AntN ...
        t1 |Ant1 109 chans XX | Ant1 109 chans YY| Ant2 ... | AntN ...
        These are REAL VALUED (1x float)

        .LC files store ant1 XY and are upper triangular, so
        1x1y| 1x2x | 1x2y | ... | 1x32y | 
              2x2y | 2x3x | ... |  ...  |
                     3x3y | ... |  ...  |
        These are COMPLEX VALUED (2xfloats)

        """

        h1("Opening L-file")
        filename = self.filename.rstrip('.LA').rstrip('.LC')
        if config_xml is None:
            config_xml = filename + '.xml'
        try:
            self.xmlData = etree.parse(config_xml)
        except IOError:
            print "ERROR: Cannot open %s" % config_xml
            exit()

        # Load visibility data
        h2("Loading visibility data")
        vis = self._readLfile()

        h2("Generating baseline IDs")
        # Create baseline IDs using MIRIAD >255 antenna format (which sucks)
        bls, ant_arr = self.generateBaselineIds(n_ant)

        bl_lower = []
        for dd in range(vis.shape[0]):
            bl_lower += bls

        h2("Converting visibilities to FLUX columns")
        flux = np.zeros([len(bl_lower), n_chans * n_stk * 2], dtype='float32')
        for ii in range(len(bl_lower)):
            ant1, ant2 = ant_arr[ii]
            idx1, idx2 = 2 * (ant1 - 1), 2 * (ant2 - 1)
            xx = vis[0, idx1, idx2]
            yy = vis[0, idx1 + 1, idx2 + 1]
            xy = vis[0, idx1, idx2 + 1]
            yx = vis[0, idx1 + 1, idx2]
            flux[ii] = np.column_stack((xx, yy, xy, yx)).flatten()

        self.d_uv_data["BASELINE"] = bl_lower
        self.d_uv_data["FLUX"] = flux

        h1("Generating FITS-IDI schema from XML")
        hdu_primary = make_primary(config=config_xml)
        tbl_array_geometry = make_array_geometry(config=config_xml, num_rows=n_ant)
        tbl_antenna = make_antenna(config=config_xml, num_rows=n_ant)
        tbl_frequency = make_frequency(config=config_xml, num_rows=1)
        tbl_source = make_source(config=config_xml, num_rows=1)

        #h1('Creating HDU list')
        hdulist = pf.HDUList(
            [hdu_primary,
             tbl_array_geometry,
             tbl_frequency,
             tbl_antenna,
             tbl_source
            ])
        #print hdulist.info()
        #hdulist.verify()

        # We are now ready to back-fill Interfits dictionaries using readfitsidi
        self.fits = hdulist
        self.readFitsidi(from_file=False, load_uv_data=False)

        h2("Populating interfits dictionaries")
        # Create interfits dictionary entries
        self.setDefaultsLeda(n_uv_rows=len(bl_lower))
        self.telescope = self.h_uv_data["TELESCOP"]
        self.source    = self.d_source["SOURCE"][0]


    def readDada(self, n_int=10,  n_stk=4, xmlbase=None, header_dict=None, data_arr=None):
        """ Read a LEDA DADA file.

        header_dict (dict): psrdada header. Defaults to None. If a dict is passed, then instead of
                            loading data from file, data will be loaded from data_arr
        data_arr (np.ndarray): data array. This should be a preformatted FLUX data array.
        """
        h1("Loading DADA data")
        if type(header_dict) is dict:
            h2("Loading from shared memory")
            d = HeaderDataUnit(header_dict, data_arr)
            flux = data_arr
            h2("Generating baseline IDs")
            bls, ant_arr = self.generateBaselineIds(n_ant)
            bl_lower = []
            while len(bl_lower) < len(flux):
                bl_lower += bls
        else:
            h2("Loading visibility data")
            d   = dada.DadaSubBand(self.filename, n_int=n_int)
            vis = d.data
            try:
                n_chans = int(d.header["NCHAN"])
                n_pol   = int(d.header["NPOL"])
                n_ant   = int(d.header["NSTATION"])
            except ValueError:
                print "WARNING: Cannot load NCHAN / NPOL / NSTATION from dada file"
                raise

            h2("Generating baseline IDs")
            bls, ant_arr = self.generateBaselineIds(n_ant)
            bl_lower = []
            for dd in range(vis.shape[0] / n_int):
                bl_lower += bls

        if not header_dict:
            h2("Converting visibilities to FLUX columns")
            # Need to convert into real & imag
            # Not a major performance hit as these are super fast.
            re_vis = np.real(vis)
            im_vis = np.imag(vis)
            # assert np.allclose(vis, re_vis + np.complex(1j)*im_vis)
            flux = np.zeros([len(bl_lower) * n_int, n_chans * n_stk * 2], dtype='float32')
            for int_num in range(n_int):
                idx = int_num * len(bl_lower)
                for ii in range(len(bl_lower)):
                    ant1, ant2 = ant_arr[ii]
                    ant1, ant2 = ant1 - 1, ant2 - 1
                    re_xx = re_vis[int_num, ant1, ant2, :, 0, 0]
                    re_yy = re_vis[int_num, ant1, ant2, :, 0, 1]
                    re_xy = re_vis[int_num, ant1, ant2, :, 1, 0]
                    re_yx = re_vis[int_num, ant1, ant2, :, 1, 1]
                    im_xx = im_vis[int_num, ant1, ant2, :, 0, 0]
                    im_yy = im_vis[int_num, ant1, ant2, :, 0, 1]
                    im_xy = im_vis[int_num, ant1, ant2, :, 1, 0]
                    im_yx = im_vis[int_num, ant1, ant2, :, 1, 1]

                    flux[idx + ii] = np.column_stack((re_xx, im_xx, re_yy, im_yy, re_xy, im_xy, re_yx, im_yx)).flatten()
            #print flux.shape

        self.d_uv_data["BASELINE"] = np.array([bl_lower for ii in range(n_int)]).flatten()
        self.d_uv_data["FLUX"] = flux


        h1("Generating FITS-IDI schema from XML")
        if xmlbase is None:
            dirname, filename = os.path.split(os.path.abspath(__file__))
            xmlbase = os.path.join(dirname, 'config/config.xml')
        self.xmlData = etree.parse(xmlbase)

        hdu_primary        = make_primary(config=self.xmlData)
        tbl_array_geometry = make_array_geometry(config=self.xmlData, num_rows=n_ant)
        tbl_antenna        = make_antenna(config=self.xmlData, num_rows=n_ant)
        tbl_frequency      = make_frequency(config=self.xmlData, num_rows=1)
        tbl_source         = make_source(config=self.xmlData, num_rows=1)

        #h1('Creating HDU list')
        hdulist = pf.HDUList(
            [hdu_primary,
             tbl_array_geometry,
             tbl_frequency,
             tbl_antenna,
             tbl_source
            ])
        #print hdulist.info()
        #hdulist.verify()

        # We are now ready to back-fill Interfits dictionaries using readfitsidi
        self.fits = hdulist
        self.stokes_axis = ['XX', 'YY', 'XY', 'YX']
        self.stokes_vals = [-5, -6, -7, -8]
        self.readFitsidi(from_file=False, load_uv_data=False)

        h2("Populating interfits dictionaries")
        self.setDefaults(n_uv_rows=len(bl_lower))
        self.obs_code = ''
        self.correlator = d.header["INSTRUMENT"]
        self.instrument = d.header["INSTRUMENT"]
        self.telescope  = d.header["TELESCOPE"]
        
        # Compute time offset
        h2("Computing UTC offsets")
        dt_obj = datetime.strptime(d.header["UTC_START"], "%Y-%m-%d-%H:%M:%S")
        byte_offset = int(d.header["OBS_OFFSET"])
        bytes_per_avg = int(d.header["BYTES_PER_AVG"])
        num_int = byte_offset / bytes_per_avg
        time_offset = num_int * leda_config.INT_TIME
        dt_obj = dt_obj + timedelta(seconds=time_offset)
        date_obs = dt_obj.strftime("%Y-%m-%dT%H:%M:%S")
        dd_obs   = dt_obj.strftime("%Y-%m-%d")
        
        print "UTC START:   %s"%d.header["UTC_START"]
        print "TIME OFFSET: %s"%timedelta(seconds=time_offset)
        print "NEW START:   %s"%date_obs
        
        self.date_obs   = date_obs
        self.h_params["NSTOKES"] = 4
        self.h_params["NBAND"]   = 1
        self.h_params["NCHAN"]   = int(d.header["NCHAN"])
        self.h_common["NO_CHAN"] = int(d.header["NCHAN"])
        self.h_common["REF_FREQ"] = float(d.header["CFREQ"]) * 1e6
        self.h_common["CHAN_BW"]  = float(d.header["BW"]) * 1e6 / self.h_params["NCHAN"]
        self.h_common["REF_PIXL"] = self.h_params["NCHAN"] / 2
        self.h_common["RDATE"]    = dd_obs  # Ignore time

        self.d_frequency["CH_WIDTH"]  = self.h_common["CHAN_BW"]
        self.d_frequency["TOTAL_BANDWIDTH"] = float(d.header["BW"]) * 1e6
        self.stokes_axis = ['XX', 'YY', 'XY', 'YX']
        self.stokes_vals = [-5, -6, -7, -8]

        self.d_array_geometry["ANNAME"] = ["Stand%03d"%i for i in range(len(self.d_array_geometry["ANNAME"]))]
        self.d_array_geometry["NOSTA"]  = [i for i in range(len(self.d_array_geometry["NOSTA"]))]

        print d.header
        #print self.d_frequency
        #print self.h_common
        #print self.h_params

    def setDefaultsLeda(self, n_uv_rows=None):
        """ set LEDA specific default values """
        if n_uv_rows is None:
            n_uv_rows = len(self.d_uv_data["BASELINE"])
            
        self.setDefaults(n_uv_rows)

        self.d_frequency["CH_WIDTH"]        = self.s2arr(leda_config.CH_WIDTH)
        self.d_frequency["TOTAL_BANDWIDTH"] = self.s2arr(leda_config.SUB_BW)
        self.h_uv_data["TELESCOP"]          = leda_config.TELESCOP
        self.h_array_geometry["ARRNAM"]     = leda_config.ARRNAM
    
    def loadAntArr(self):
        """ Loads ANTENNA and ARRAY_GEOMETRY tables as set in leda_config """
        h1("Loading ANTENNA and ARRAY_GEOMETRY from JSON")
        self.h_array_geometry = load_json(leda_config.json_h_array_geometry)
        self.d_array_geometry = load_json(leda_config.json_d_array_geometry)
        self.h_antenna        = load_json(leda_config.json_h_antenna)     
        self.d_antenna        = load_json(leda_config.json_d_antenna)      

    def computeSiderealTime(self, ts=None):
        """ Computes the LST for a given timestamp.

        ts (float): Timestamp to use. If none is given, DATE-OBS value
                    is used in lieu.

        returns LST in degrees
        """

        h2("Computing LST from UTC")
        if ts is not None:
            dt_utc = datetime.utcfromtimestamp(ts)
        else:
            dt_utc = datetime.strptime(self.date_obs, "%Y-%m-%dT%H:%M:%S")
        ov = leda_config.ovro
        ov.date = dt_utc
        lst, lst_deg = ov.sidereal_time(), ov.sidereal_time() / 2 / np.pi * 360
        print "UTC: %s"%dt_utc
        print "LST: %s (%s)"%(lst, lst_deg)
        return lst_deg

    def generateBaselineIds(self, n_ant=None):
        """ Generate a list of unique baseline IDs and antenna pairs
        
        This uses the MIRIAD definition for >256 antennas:
        bl_id = 2048*ant1 + ant2 + 65536
        
        n_ant: number of antennas in the array
        """
        if n_ant is None:
            n_ant = self.n_ant

        bls, ant_arr = [], []
        for ii in range(1, n_ant + 1):
            for jj in range(1, n_ant + 1):
                if jj >= ii:
                    ant_arr.append((ii, jj))
                    if ii > 255 or jj > 255:
                        bl_id = ii * 2048 + jj + 65536
                    else:
                        bl_id = 256 * ii + jj
                    bls.append(bl_id)
        return bls, ant_arr

    def generateUVW(self, src='ZEN', update_src=True, conjugate=False, use_stored=False):
        """ Generate UVW coordinates based on timestamps and array geometry

        Updates UVW coordinates to phase to a given source. Uses pyEphem observer
        along with methods is lib.uvw for computations

        src (str): Source to phase to. Sources are three capital letters:
            ZEN: Zenith (RA will be computed from timestamps)
            CYG: Cygnus A
            CAS: Cassiopeia A
            TAU: Taurus A
            VIR: Virgo A

        use_stored (bool): If True, uses stored UVW coordinates (does not recompute).
                           this is faster than recomputing.
        update_src (bool): Default True, update the SOURCE table.
        conjugate (bool): Conjuagte UVW coordinates? Do this if things are flipped in map.
        """

        h1("Generating UVW coordinates")
        # First, compute LST
        tt_source = time.strptime(self.date_obs, "%Y-%m-%dT%H:%M:%S")
        ts_source = calendar.timegm(tt_source)
        lst_deg = self.computeSiderealTime(ts_source)

        # Find HA and DEC of source
        if src.upper() == 'ZEN':
            H, d = 0, np.deg2rad(float(leda_config.latitude))
            dec_deg  = float(leda_config.latitude)
            ra_deg   = lst_deg
        else:
            try:
                src_names = leda_config.src_names
                src_ras   = leda_config.src_ras
                src_decs  = leda_config.src_decs
                idx = src_names.index(src.upper())
                h2("Phasing to %s"%src_names[idx])
                ra_deg, dec_deg = src_ras[idx], src_decs[idx]

                # Now we have the RA and DEC, need to convert into hour angle
                H = np.deg2rad(lst_deg - ra_deg)
                d = np.deg2rad(dec_deg)

            except ValueError:
                print "Error: Cannot phase to %s"%src
                print "Choose one of: ZEN, ",
                for src in src_names:
                    print "%s, "%src,
                print "\n"
                raise

        print "LST:        %2.3f deg"%lst_deg
        print "Source RA:  %2.3f deg"%ra_deg
        print "Source DEC: %2.3f deg"%dec_deg
        print "HA:         %2.3f deg"%np.rad2deg(H)

        # Recreate list of baselines
        h2("Computing UVW coordinates for %s"%src)
        bl_ids, ant_arr = self.generateBaselineIds(self.n_ant)
        n_iters = int(len(self.d_uv_data["BASELINE"]) / len(bl_ids))
        
        if use_stored:
            h2("Loading stored values")
            self.loadUVW()
        else:
            n_ant = len(self.d_array_geometry['ANNAME'])
            xyz   = self.d_array_geometry['STABXYZ']

            # Pre-compute UVW tranformation matrix
            sin, cos = np.sin, np.cos
            try:
                assert H < 2 * np.pi and d < 2 * np.pi
            except AssertionError:
                raise ValueError("HA and DEC are too large (may not be in radians).")
            t_matrix = np.matrix([
              [sin(H), cos(H), 0],
              [-sin(d)*cos(H), sin(d)*sin(H), cos(d)],
              [cos(d)*cos(H), -cos(d)*sin(H), sin(H)]
            ])

            # Compute baseline vectors
            bl_veclist, uvw_list = [], []
            for ant_pair in ant_arr:
                ii, jj = ant_pair[0] - 1, ant_pair[1] - 1
                bl_vec = xyz[ii] - xyz[jj]
                bl_veclist.append(bl_vec)
                uvw_list.append(uvw.computeUVW(bl_vec, H, d, conjugate=conjugate, t_matrix=t_matrix))
            uvw_arr = np.array(uvw_list)
            
            # Fill with data
            uu, vv, ww = [], [], []
            for ii in range(n_iters):
                uu.append(uvw_arr[:, 0])
                vv.append(uvw_arr[:, 1])
                ww.append(uvw_arr[:, 2])
            self.d_uv_data["UU"]   = np.array(uu).ravel()
            self.d_uv_data["VV"]   = np.array(vv).ravel()
            self.d_uv_data["WW"]   = np.array(ww).ravel()
            
        h2("Generating timestamps")
        dd, tt = [], []        
        for ii in range(n_iters):
            jd, jt = uvw.convertToJulianTuple(self.date_obs)
            tdelta = leda_config.INT_TIME * ii
            jds = [jd for jj in range(len(ant_arr))]
            jts = [jt + tdelta for jj in range(len(ant_arr))]
            dd.append(jds)
            tt.append(jts)

        self.d_uv_data["DATE"] = np.array(dd, dtype='float64').ravel()
        self.d_uv_data["TIME"] = np.array(tt, dtype='float64').ravel()

        if update_src:
            h2("Updating SOURCE table")
            self.d_source["SOURCE"] = self.s2arr(src)
            self.d_source["RAEPO"]  = self.s2arr(ra_deg)
            self.d_source["DECEPO"] = self.s2arr(dec_deg)

        #print self.d_uv_data["DATE"]
        #print self.d_uv_data["TIME"]
        #print np.array(uu).shape, np.array(vv).shape, np.array(ww).shape

    def dumpUVW(self, filename):
        """ Dump precomputed UVW coordinates to file 
        
        filename (str): name of output file (.json format)
        """
        d ={}
        d["UU"]       = self.d_uv_data["UU"]
        d["VV"]       = self.d_uv_data["VV"]
        d["WW"]       = self.d_uv_data["WW"]
        d["BASELINE"] = self.d_uv_data["BASELINE"]
        
        h2("Dumping UVW coords to %s"%filename)
        dump_json(d, filename)
    
    def loadUVW(self, filename=None):
        """ Load precomputed UVW coordinates from file
        
        filename (str): name of input file. If not set, uses default
                        from leda_config file.
        """
        
        h2("Loading UVW coordinates from file")
        if not filename:
            filename = leda_config.json_uvw_coordinates
        d = load_json(filename)
        self.d_uv_data["UU"]       = d["UU"]      
        self.d_uv_data["VV"]       = d["VV"]      
        self.d_uv_data["WW"]       = d["WW"]      
        self.d_uv_data["BASELINE"] = d["BASELINE"]
    
    def leda_set_value(self, key, value):
        """ Set values which are commonly incorrect from uvfits writer """

        if key == 'ARRNAM':
            self.h_array_geometry['ARRNAM'] = value
        if key == 'INTTIM':
            self.d_uv_data['INTTIM'][:] = value
        if key == 'TELESCOP':
            self.h_uv_data['TELESCOP'] = value

    def readAntennaLocations(self, filename):
        """ Read antenna locations file """
        atab = np.genfromtxt(filename, comments='#', dtype='str')
        stand_ids  = atab[:, 0]
        stand_east = atab[:, 1].astype('float')
        stand_west = atab[:, 2].astype('float')
        stand_elev = atab[:, 3].astype('float')

    def remove_miriad_baselines(self):
        """ Remove baseline data for all antennas with IDs > 255

        Miriad-type UVFITS files use the convention
            ant1*256+ant2 if ants < 255
            ant1*2048+ant2+65536 if ants >255
        The miriad convention screws up import into many reduction packages.
        """

        h1("Removing MIRIAD baselines")
        bls = np.array(self.d_uv_data["BASELINE"])

        if self.n_ant > 255:
            self.n_ant = 255

        max_bl = 255 * 256 + 255
        ok_bls = bls < max_bl
        #print ok_bls
        for k in self.d_uv_data.keys():
            try:
                self.d_uv_data[k] = np.array(self.d_uv_data[k])
                self.d_uv_data[k] = self.d_uv_data[k][ok_bls]
                #print len(self.d_uv_data[k])
            except TypeError:
                print k
                print self.d_uv_data[k]
                raise
            except ValueError:
                print k
                print len(self.d_uv_data[k]),
                print len(ok_bls)
                raise

        #for k in self.d_antenna.keys():
        #    self.d_antenna[k] = self.d_antenna[k][0:self.n_ant]
        #    #    print len(self.d_antenna[k])

        #for k in self.d_array_geometry.keys():
        #    self.d_array_geometry[k] = self.d_array_geometry[k][0::self.n_ant]
        #    #    print len(self.d_array_geometry[k])

        h2("Fixing NOPCAL (setting to zero)")
        self.h_antenna["NOPCAL"] = 0

        h2("Setting INTTIME to %s"%INT_TIME)
        self.d_uv_data["INTTIM"] = np.ones_like(self.d_uv_data["INTTIM"]) * INT_TIME

    def flag_antenna(self, antenna_id, reason=None, severity=0):
        """ Flag antenna as bad

        antenna_id (int): ID of antenna to flag. Starts at 1.
        reason (str): Defaults to None; short (<24 char) reason for flagging
        severity (int): -1 Not assigned, 0 Known bad, 1 Probably bad, 2 Maybe bad
        """
        h2("Flagging antenna %i"%antenna_id)

        flag_keywords = ['SOURCE_ID', 'ARRAY', 'ANTS', 'FREQID', 'BANDS', 'CHANS', 'PFLAGS', 'REASON', 'SEVERITY']
        try:
            self.d_flag["SOURCE_ID"]
        except KeyError:
            for k in flag_keywords:
                self.d_flag[k] = []

        if reason is None:
            reason = "Known bad antenna."

        flag_k_zeros  = ['SOURCE_ID', 'ARRAY', 'FREQID']
        flag_k_one    = ['BANDS']
        for k in flag_k_one:
            self.d_flag[k].append(1)
        for k in flag_k_zeros:
            self.d_flag[k].append(0)

        self.d_flag["ANTS"].append((antenna_id, 0))
        self.d_flag["REASON"].append(reason)
        self.d_flag["PFLAGS"].append((1,1,1,1))
        self.d_flag["SEVERITY"].append(severity)
        self.d_flag["CHANS"].append((0, 4096))

    def apply_cable_delays(self):
        """ Apply antenna cable delays

        Each cable introduces a phase shift of
            phi = 2 pi f t
        Visibility is VpVq*, so we need to apply
            exp(-i  (phip - phiq))
        to compensate for cable delay

        TODO: Phase values can be precomputed; implement this for speed boost.
        """

        h1("Applying cable delays")
        #t0 = time.time()
        # Load antenna Electrical Lengths
        sol   = leda_config.SPEED_OF_LIGHT
        els   = load_json(leda_config.json_antenna_el_lens)["EL"]
        tdelts = els / sol

        # Generate frequency array from metadata
        ref_delt = self.h_common["CHAN_BW"]
        ref_pix  = self.h_common["REF_PIXL"]
        ref_val  = self.h_common["REF_FREQ"]
        num_pix  = self.h_common["NO_CHAN"]
        freqs    = np.arange(0,num_pix,1) * ref_delt + (ref_val - ref_pix * ref_delt)
        print freqs.shape
        # Compute phase delay for each antenna pair
        try:
            assert self.d_uv_data["FLUX"].dtype == 'float32'
        except AssertionError:
            print self.d_uv_data["FLUX"].dtype
            raise

        flux  = self.d_uv_data["FLUX"].view('complex64')

        bls, ant_arr = self.generateBaselineIds()
        w = 2 * np.pi * freqs # Angular freq

        #t1 = time.time()
        for ii in range(len(bls)):
            ant1, ant2 = ant_arr[ii]
            bl         = bls[ii]
            td1, td2   = tdelts[ant1-1], tdelts[ant2-1]

            # Compute phases for X and Y pol on antennas A and B
            pxa, pya, pxb, pyb = w * td1[0], w * td1[1], w * td2[0], w * td2[1]

            # Corrections require negative sign (otherwise reapplying delays)
            e_xx = np.exp(-1j * (pxa - pxb))
            e_yy = np.exp(-1j * (pya - pyb))
            e_xy = np.exp(-1j * (pxa - pyb))
            e_yx = np.exp(-1j * (pya - pxb))

            phase_corrs = np.column_stack((e_xx, e_yy, e_xy, e_yx)).flatten()

            #if not ii%5000:
            #    plt.subplot(211)
            #    plt.plot(np.angle(phase_corrs[::4]))
            #    plt.subplot(212)
            #    plt.plot(np.angle(flux[ii][::4]))
            flux[ii] = flux[ii] * phase_corrs
            #if not ii%5000:
            #    plt.plot(np.angle(flux[ii][::4]))
            #    plt.xlabel("Ant %s, Ant %s"%(ant1, ant2))
            #    plt.show()
            #    print td1, td2
            #    #exit()
            #print flux[ii].shape
            #print phases.shape
#
            #time.sleep(2)
        #t2 = time.time()

        #print t2 - t0
        #print t2 - t1







