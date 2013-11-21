#! /usr/bin/env python
# encoding: utf-8
"""
ledafits.py
============

Extension of interfits.py with LEDA-OVRO specific methods, such as ability to read
L-files and DADA files.
"""

import time
from datetime import datetime
from interfits import *
from lib import dada, uvw
import leda_config

# Load globals from config file
OFFSET_DELTA, INT_TIME, N_INT = leda_config.OFFSET_DELTA, leda_config.INT_TIME, leda_config.N_INT

class LedaFits(InterFits):
    """ LEDA extension of InterFits class 
    
    This adds ability to read LA, LC and dada files, and adds helper functions
    for computing UVW coordinates, generating timestamps, and computing zenith RA/DEC.
    """

    def readFile(self, filename):
        # Check what kind of file to load
        if filename:
            matched = False
            regex = '([0-9A-Za-z-_]+).uvfits'
            match = re.search(regex, filename)
            if match:
                self._readFile('uvfits')
                matched = True

            regex = '([0-9A-Za-z-_]+).fitsidi'
            match = re.search(regex, filename)
            if match:
                self._readFile('fitsidi')
                matched = True

            regex = '([0-9A-Za-z-_]+).hdf'
            match = re.search(regex, filename)
            if match:
                self._readFile('hdf5')
                matched = True

            regex = '([0-9A-Za-z-_]+).json'
            match = re.search(regex, filename)
            if match:
                self._readFile('json')
                matched = True

            regex = '([0-9A-Za-z-_]+).LA'
            match = re.search(regex, filename)
            if match:
                self._readFile('lfile')
                matched = True

            regex = '([0-9A-Za-z-_]+).LC'
            match = re.search(regex, filename)
            if match:
                self._readFile('lfile')
                matched = True

            regex = '([0-9A-Za-z-_]+).dada'
            match = re.search(regex, filename)
            if match:
                self._readFile('dada')
                matched = True

            if matched == False:
                raise IOError("Unknown file format: %s"%filename)

    def _readFile(self, filetype):
        """ Lookup dictionary (case statement) for file types """
        return {
            'uvfits': self.readUvfits,
            'fitsidi': self.readFitsidi,
            'hdf5': self.readHdf5,
            'lfile': self.readLfile,
            'json': self.readJson,
            'dada': self.readDada
        }.get(filetype)()

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


    def readDada(self, n_ant=256, n_pol=2, n_chans=109, n_stk=4, xmlbase=None):
        """ Read a LEDA DADA file
        """

        h2("Loading visibility data")
        d = dada.DadaSubBand(self.filename)
        vis = d.data
        #print d.header
        print vis.shape

        # Need to convert into real & imag
        # Not a major performance hit as these are super fast.
        re_vis = np.real(vis)
        im_vis = np.imag(vis)
        # assert np.allclose(vis, re_vis + np.complex(1j)*im_vis)

        h2("Generating baseline IDs")
        bls, ant_arr = self.generateBaselineIds(n_ant)
        bl_lower = []
        for dd in range(vis.shape[0]):
            bl_lower += bls

        h2("Converting visibilities to FLUX columns")
        flux = np.zeros([len(bl_lower), n_chans * n_stk * 2])
        for ii in range(len(bl_lower)):
            ant1, ant2 = ant_arr[ii]
            ant1, ant2 = ant1 - 1, ant2 - 1
            re_xx = re_vis[0, ant1, ant2, :, 0, 0]
            re_yy = re_vis[0, ant1, ant2, :, 0, 1]
            re_xy = re_vis[0, ant1, ant2, :, 1, 0]
            re_yx = re_vis[0, ant1, ant2, :, 1, 1]
            im_xx = im_vis[0, ant1, ant2, :, 0, 0]
            im_yy = im_vis[0, ant1, ant2, :, 0, 1]
            im_xy = im_vis[0, ant1, ant2, :, 1, 0]
            im_yx = im_vis[0, ant1, ant2, :, 1, 1]
            flux[ii] = np.column_stack((re_xx, im_xx, re_yy, im_yy, re_xy, im_xy, re_yx, im_yx)).flatten()
        #print flux.shape
        self.d_uv_data["BASELINE"] = bl_lower
        self.d_uv_data["FLUX"] = flux

        h1("Generating FITS-IDI schema from XML")
        if xmlbase is None:
            dirname, filename = os.path.split(os.path.abspath(__file__))
            xmlbase = os.path.join(dirname, 'config/config.xml')
        self.xmlData = etree.parse(xmlbase)

        hdu_primary = make_primary(config=self.xmlData)
        tbl_array_geometry = make_array_geometry(config=self.xmlData, num_rows=n_ant)
        tbl_antenna = make_antenna(config=self.xmlData, num_rows=n_ant)
        tbl_frequency = make_frequency(config=self.xmlData, num_rows=1)
        tbl_source = make_source(config=self.xmlData, num_rows=1)

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
        
        # Convert UTC_START into dateTtime format
        dd_obs = d.header["UTC_START"][0:10]
        tt_obs = d.header["UTC_START"][11:]
        date_obs = dd_obs + 'T' + tt_obs
        
        self.date_obs   = date_obs
        self.h_params["NSTOKES"] = 4
        self.h_params["NBAND"]   = 1
        self.h_params["NCHAN"]   = int(d.header["NCHAN"])
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
        print self.d_frequency
        print self.h_common
        print self.h_params

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
        

    def generateBaselineIds(self, n_ant):
        """ Generate a list of unique baseline IDs and antenna pairs
        
        This uses the MIRIAD definition for >256 antennas:
        bl_id = 2048*ant1 + ant2 + 65536
        
        n_ant: number of antennas in the array
        """
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

    def generateUVW(self, conjugate=False, H=0, d=float(leda_config.latitude), use_stored=False):
        """ Generate UVW coordinates based on timestamps and array geometry

        Uses pyEphem observer + lib.uvw for computations
        
        H (float): Hour angle, in degrees. Defaults to 0
        d (float): Declination, in degrees. Defaults to LWA-OVRO declination.
        use_stored (bool): If True, uses stored UVW coordinates (does not recompute).
                           this is faster than recomputing. 
        """
        

        # Recreate list of baselines
        h1("Generating UVW coordinates")
        bl_ids, ant_arr = self.generateBaselineIds(self.n_ant)
        n_iters = int(len(self.d_uv_data["BASELINE"]) / len(bl_ids))
        
        if use_stored:
            h2("Loading stored values")
            self.loadUVW()
        else:
            n_ant = len(self.d_array_geometry['ANNAME'])
            xyz   = self.d_array_geometry['STABXYZ']
            
            # Compute baseline vectors
            bl_veclist, uvw_list = [], []
            for ant_pair in ant_arr:
                ii, jj = ant_pair[0] - 1, ant_pair[1] - 1
                bl_vec = xyz[ii] - xyz[jj]
                bl_veclist.append(bl_vec)
                uvw_list.append(uvw.computeUVW(bl_vec, H, d, conjugate=conjugate))
            uvw_arr = np.array(uvw_list)
            
            # Fill with data
            h2("Computing UV_DATA coordinates")
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
            jds = [jd for ii in range(len(ant_arr))]
            jts = [jt for ii in range(len(ant_arr))]
            dd.append(jds)
            tt.append(jts)
            
        self.d_uv_data["DATE"] = np.array(dd).ravel().astype('float64')
        self.d_uv_data["TIME"] = np.array(tt).ravel().astype('float64')
        
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
        #self.n_ant = self.n_ant - 1

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

        #for k in self.d_antenna.keys():
        #    self.d_antenna[k] = self.d_antenna[k][0:256]
        #    #    print len(self.d_antenna[k])

        #for k in self.d_array_geometry.keys():
        #    self.d_array_geometry[k] = self.d_array_geometry[k][0:256]
        #    #    print len(self.d_array_geometry[k])

        h2("Fixing NOPCAL (setting to zero)")
        self.h_antenna["NOPCAL"] = 0

        h2("Setting INTTIME to %s"%INT_TIME)
        self.d_uv_data["INTTIM"] = np.ones_like(self.d_uv_data["INTTIM"]) * INT_TIME