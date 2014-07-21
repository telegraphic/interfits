# -*- coding: utf-8 -*-

"""
ledafits.py
============

Extension of interfits.py with LEDA-OVRO specific methods, such as ability to read
L-files and DADA files.
"""

import calendar
from datetime import datetime, timedelta
import ephem
import re
import os
from lxml import etree
import numpy as np
import pyfits as pf
from lib.json_numpy import *

from interfits import *
from lib import dada, coords
from lib.pyFitsidi import *
import ledafits_config

__version__ = '0.0'
__all__ = ['HeaderDataUnit', 'LedaFits', '__version__', '__all__']


class HeaderDataUnit(object):
    """ Very basic object with header and data units """

    def __init__(self, header, data):
        self.header = header
        self.data = data


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

    def _readLfile(self, n_ant=32, n_pol=2, n_chans=600):
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
        self.pp.debug(lfa.shape)

        # Cross-corrs
        #h2("Opening cross-corrs (.LC)")
        lfc = np.fromfile(filename + '.LC', dtype='float32')
        lfc = lfc.reshape([len(lfc) / n_chans / n_blcc / 2, n_blcc, n_chans, 2])
        self.pp.debug(lfc.shape)

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
        self.pp.debug("\n %s" % vis.shape[0])

        return vis

    def readLfile(self, n_ant=32, n_chans=600, n_stk=4, config_xml=None):
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

        self.pp.h1("Opening L-file")
        filename = self.filename.rstrip('.LA').rstrip('.LC')
        if config_xml is None:
            config_xml = filename + '.xml'
        try:
            self.xmlData = etree.parse(config_xml)
        except IOError:
            self.pp.err("ERROR: Cannot open %s" % config_xml)
            exit()

        # Load visibility data
        self.pp.h2("Loading visibility data")
        vis = self._readLfile(n_ant=n_ant, n_pol=2, n_chans=n_chans)

        self.pp.h2("Generating baseline IDs")
        # Create baseline IDs using MIRIAD >255 antenna format (which sucks)
        bls, ant_arr = coords.generateBaselineIds(n_ant)

        bl_lower = []
        for dd in range(vis.shape[0]):
            bl_lower += bls

        self.pp.h2("Converting visibilities to FLUX columns")
        flux = np.zeros([len(bl_lower), n_chans * n_stk * 2], dtype='float32')
        for ii in range(len(bl_lower)):
            ant1, ant2 = ant_arr[ii % len(ant_arr)]
            idx1, idx2 = 2 * (ant1 - 1), 2 * (ant2 - 1)
            xx = vis[0, idx1, idx2]
            yy = vis[0, idx1 + 1, idx2 + 1]
            xy = vis[0, idx1, idx2 + 1]
            yx = vis[0, idx1 + 1, idx2]
            flux[ii] = np.column_stack((xx, yy, xy, yx)).flatten()

        self.d_uv_data["BASELINE"] = bl_lower
        self.d_uv_data["FLUX"] = flux

        self.pp.h1("Generating FITS-IDI schema from XML")
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
             tbl_source])

        # We are now ready to back-fill Interfits dictionaries using readfitsidi
        self.fits = hdulist
        self.readFitsidi(from_file=False, load_uv_data=False)

        self.pp.h2("Populating interfits dictionaries")
        # Create interfits dictionary entries
        self.setDefaultsLeda(n_uv_rows=len(bl_lower))
        self.telescope = self.h_uv_data["TELESCOP"]
        self.source = self.d_source["SOURCE"][0]

    def readDada(self, n_int=None, xmlbase=None, header_dict=None, data_arr=None, inspectOnly=False):
        """ Read a LEDA DADA file.

            header_dict (dict): psrdada header. Defaults to None. If a dict is passed, then instead of
                                loading data from file, data will be loaded from data_arr
            data_arr (np.ndarray): data array. This should be a preformatted FLUX data array.
            """

        self.pp.h1("Loading DADA data")
        if type(header_dict) is dict:
            self.pp.h2("Loading from shared memory")
            d = HeaderDataUnit(header_dict, data_arr)
            flux = data_arr
            self.pp.h2("Generating baseline IDs")
            bls, ant_arr = coords.generateBaselineIds(d.n_ant)
            bl_lower = []
            while len(bl_lower) < len(flux):
                bl_lower += bls
        else:
            self.pp.h2("Loading visibility data")
            d = dada.DadaReader(self.filename, n_int, inspectOnly=inspectOnly)
            vis = d.data
            self.dada_header = d.header
            try:
                n_ant = d.n_ant
                n_int = d.n_int
                self.n_ant = n_ant
            except ValueError:
                raise RuntimeError("Cannot load NCHAN / NPOL / NSTATION from dada file")

        if not header_dict:
            self.pp.h2("Converting visibilities to FLUX columns")
            do_remap = False
            if d.header["TELESCOPE"] in ('LEDA', 'LWAOVRO', 'LWA-OVRO', 'LEDAOVRO', 'LEDA512', 'LEDA-OVRO'):
                do_remap = False
            flux = self._vis_matrix_to_flux(vis, remap=do_remap)
            bls, ant_arr = coords.generateBaselineIds(n_ant)
            bl_lower = []
            for dd in range(vis.shape[0] / n_int):
                bl_lower += bls

        self.d_uv_data["BASELINE"] = np.array([bl_lower for ii in range(n_int)]).flatten()
        self.d_uv_data["FLUX"] = flux

        self.pp.h1("Generating FITS-IDI schema from XML")
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
             tbl_source])

        self.fits = hdulist

        self.stokes_vals = [-5, -6, -7, -8]
        self.readFitsidi(from_file=False, load_uv_data=False)
        #hdulist.verify()

        self.pp.h2("Populating interfits dictionaries")
        self.setDefaults(n_uv_rows=len(bl_lower * n_int))
        self.obs_code = ''
        self.correlator = d.header["INSTRUMENT"]
        self.instrument = d.header["INSTRUMENT"]
        self.telescope = d.header["TELESCOPE"]

        # Compute the integration time
        tsamp = float(d.header["TSAMP"]) * 1e-6  # Sampling time per channel, in microseconds
        navg = int(d.header["NAVG"])  # Number of averages per integration
        int_tim = tsamp * navg  # Integration time is tsamp * navg
        self.t_int = d.t_int

        # Compute time offset
        self.pp.h2("Computing UTC offsets")
        dt_obj = datetime.strptime(d.header["UTC_START"], "%Y-%m-%d-%H:%M:%S")
        time_offset = d.t_offset  # Time offset since observation began
        dt_obj = dt_obj + timedelta(seconds=time_offset)
        date_obs = dt_obj.strftime("%Y-%m-%dT%H:%M:%S")
        dd_obs = dt_obj.strftime("%Y-%m-%d")

        self.pp.pp("UTC START:   %s" % d.header["UTC_START"])
        self.pp.pp("TIME OFFSET: %s" % timedelta(seconds=time_offset))
        self.pp.pp("NEW START:   %s" % date_obs)

        self.date_obs = date_obs
        self.h_params["NSTOKES"] = 4
        self.h_params["NBAND"] = 1
        self.h_params["NCHAN"] = d.n_chans
        self.h_common["NO_CHAN"] = d.n_chans
        self.h_common["REF_FREQ"] = d.c_freq_mhz * 1e6
        self.h_common["CHAN_BW"] = d.chan_bw_mhz * 1e6
        self.h_common["REF_PIXL"] = d.n_chans / 2 + 1
        self.h_common["RDATE"] = dd_obs  # Ignore time component
        self.h_common["STK_1"] = -5

        self.d_frequency["CH_WIDTH"] = d.chan_bw_mhz * 1e6
        self.d_frequency["TOTAL_BANDWIDTH"] = d.bandwidth_mhz * 1e6
        self.stokes_axis = ['XX', 'YY', 'XY', 'YX']
        self.stokes_vals = [-5, -6, -7, -8]

        self.d_array_geometry["ANNAME"] = ["Stand%03d" % i for i in range(len(self.d_array_geometry["ANNAME"]))]
        self.d_array_geometry["NOSTA"] = [i for i in range(len(self.d_array_geometry["NOSTA"]))]

        self.d_uv_data["INTTIM"] = np.ones_like(self.d_uv_data["INTTIM"]) * d.t_int

        # Recreate list of baselines
        bl_ids, ant_arr = coords.generateBaselineIds(self.n_ant, autocorrs=True)
        n_iters = int(len(self.d_uv_data["BASELINE"]) / len(bl_ids))

        self.pp.h2("Generating timestamps")
        dd, tt = [], []
        for ii in range(n_iters):
            jd, jt = coords.convertToJulianTuple(self.date_obs)
            tdelta = int_tim * ii / 86400.0  # In days
            jds = [jd for jj in range(len(ant_arr))]
            jts = [jt + tdelta for jj in range(len(ant_arr))]
            dd.append(jds)
            tt.append(jts)

        self.d_uv_data["DATE"] = np.array(dd, dtype='float64').ravel()
        self.d_uv_data["TIME"] = np.array(tt, dtype='float64').ravel()

        # Load array geometry from file, based on TELESCOP name
        self.loadAntArr()
        #self.phase_to_src('ZEN')

    def _initialize_site(self):
        """ Setup site (ephem observer)

        Uses ecef2geo function (Bowring's method), to convert
        ECEF to Lat-Long-Elev, then creates an ephem observer.
        Note: this overrides the Interfits method, adding checks for LEDA telescopes,
        and preferentially loads those data.
        """
        if self.telescope in ('LEDA', 'LWAOVRO', 'LWA-OVRO', 'LEDAOVRO', 'LEDA512', 'LEDA-OVRO'):
            self.pp.h3("Data appears to be from LWAOVRO")
            self.site = ledafits_config.ovro
            self.z_elength = load_json(ledafits_config.json_antenna_el_lens)

        elif self.telescope in ('LWA1', 'LWA-1', 'LWA-NM', 'LWANM', 'LEDA64', 'LEDA64-NM'):
            self.pp.h3("Data appears to be from LWA1")
            self.site = ledafits_config.lwa1
            self.z_elength = load_json(ledafits_config.json_antenna_el_lens_nm)

        else:
            x = self.h_array_geometry["ARRAYX"]
            y = self.h_array_geometry["ARRAYY"]
            z = self.h_array_geometry["ARRAYZ"]
            lat_rad, long_rad, elev = coords.ecef2geo(x, y, z)

            self.site      = ephem.Observer()
            self.site.lon  = long_rad
            self.site.lat  = lat_rad
            self.site.elev = elev

        self.site_info(follow_verbosity=True)

    def _vis_matrix_to_flux(self, vis, remap=False):
        """Convert a visibility matrix to FITS-IDI flux standard

        Notes
        -----
        Visibility matrix should have shape:
            (n_int, ant1, ant2, chans, pola, polb)
        FITS-IDI is a flattened row on a per-baseline basis:
            (xx, yy, xy, yx)
        where each xx is (re_chan0, im_chan0, ... re_chanX, im_chanX).

        We only read one triangle of the visibility matrix, with ant1 >= ant2
        """

        # h2("Generating baseline IDs")
        n_int = vis.shape[0]
        n_ant = vis.shape[1]
        n_chans = vis.shape[3]
        n_stk = 4
        n_bls = n_ant * (n_ant - 1) / 2 + n_ant

        bls, ant_arr = coords.generateBaselineIds(n_ant)
        ant_arr0 = np.array(ant_arr) - 1  # Zero indexed
        flux = np.zeros([n_bls * n_int, n_chans * n_stk * 2], dtype='float32')

        try:
            assert vis.dtype == 'complex64'
        except AssertionError:
            raise RuntimeError('Vis data is not complex64, but is instead %s' % vis.dtype)

        for int_num in xrange(n_int):
            idx = int_num * n_bls
            vis_int = vis[int_num, ...]
            for ii in xrange(n_bls):
                ant1, ant2 = ant_arr0[ii]
                vv = vis_int[ant1, ant2, ...]
                xx = vv[:, 0, 0]
                yy = vv[:, 1, 1]
                xy = vv[:, 0, 1]
                yx = vv[:, 1, 0]
                flux[idx + ii] = np.column_stack((xx, yy, xy, yx)).flatten().view('float32')

        if remap:
            self.pp.h2("Remapping antennas")
            mapping = {
                "255A": "238A",
                "255B": "240B",
                "242B": "253B",
                "240B": "252B",
                "252A": "254A",
                "252B": "244B",
                "256B": "248B",
                "248B": "256B",
                "245B": "255B",
                "253B": "245B",
                "253A": "255A",
                "244B": "254B",
                "238A": "252A",
                "250B": "242B",
                "250A": "253A",
                "254B": "250B",
                "254A": "250A"
            }
            map_keys = set(mapping.keys())
            bls_all, ants = coords.generateBaselineIds(self.n_ant, autocorrs=True)

            for int_num in xrange(n_int):
                idx = int_num * n_bls
                vis_int = vis[int_num, ...]

                # For every antenna remapping
                for k in map_keys:
                    ant_old, pol_id = int(k[:3]), 0 if k[3] == 'B' else 1
                    ant_new = int(mapping[k][:3])

                    # Find affected baselines
                    bls_old = self.search_baselines(ant_old)
                    bls_new = self.search_baselines(ant_new)

                    # For every baseline affected
                    for bb in range(len(bls_old)):
                        # Find new antenna pair indexes
                        bl_old, bl_new = bls_old[bb], bls_new[bb]
                        bl_idx = bls_all.index(bl_old)

                        if bl_new >= 65536:
                            a1, a2 = (bl_new - 65536) / 2048 - 1, (bl_new - 65536) % 2048 - 1
                        else:
                            a1, a2 = bl_new / 256 - 1, bl_new % 256 - 1

                        # Grab all the visibility data for this baseline
                        xx = vis_int[a1, a2, :, 0, 0]
                        xy = vis_int[a1, a2, :, 0, 1]
                        yx = vis_int[a1, a2, :, 1, 0]
                        yy = vis_int[a1, a2, :, 1, 1]

                        # Now we need to figure out what we need to update
                        # Are we updating pol A or pol B? Ant1 or Ant2?
                        data = flux[idx + bl_idx]
                        sp = len(xx) * 2
                        if ant_new == a1:
                            if pol_id == 0:
                                data[0:sp] = xx.flatten().view('float32')
                                data[2 * sp:3 * sp] = xy.flatten().view('float32')
                            else:
                                data[1 * sp:2 * sp] = yy.flatten().view('float32')
                                data[3 * sp:] = yx.flatten().view('float32')
                        else:
                            if pol_id == 0:
                                data[0:sp] = xx.flatten().view('float32')
                                data[2 * sp:3 * sp] = yx.flatten().view('float32')
                            else:
                                data[1 * sp:2 * sp] = yy.flatten().view('float32')
                                data[3 * sp:] = xy.flatten().view('float32')

                        # Write this back into the right baseline
                        flux[idx + bl_idx] = data
        return flux

    def inspectFile(self, filename=None, filetype=None):
        """ Check file type, and load metadata corresponding

        filename (str): name of file. Alternatively, if a psrdada header dictionary
                        is passed, data will be loaded from shared memory. File type
                        is inferred from extension (unless filetype arg is also passed).
        filetype (str): Defaults to none. If passed, treat file as having an explicit
                        type. Useful for when extension does not match data.
        """
        # Check what kind of file to load

        if filetype is not None:
            self._inspectFile(filetype)

        else:
            if filename is None:
                pass
            elif type(filename) is tuple:
                # Tuple is header_dict and numpy data array
                head, data = filename[0], filename[1]
                return self.readDada(header_dict=head, data_arr=data, inspectOnly=True)
            else:
                file_ext = os.path.splitext(filename)[1][1:]
                self.filename = filename
                return self._inspectFile(file_ext)

    def _inspectFile(self, filetype):
        """ Lookup dictionary (case statement) for file types """
        return {
            'dada': self.inspectDada
        }.get(filetype, self.readError)()

    def inspectDada(self, n_int=None, header_dict=None, data_arr=None):
        """ Inspect a LEDA DADA file and return a dictionary describing the file contents.
            
            header_dict (dict): psrdada header. Defaults to None. If a dict is passed, then instead of
                                loading data from file, data will be loaded from data_arr
            data_arr (np.ndarray): data array. This should be a preformatted FLUX data array.
            """

        self.pp.h1("Inspecting DADA data")
        if type(header_dict) is dict:
            self.pp.h2("Inspecting from shared memory")
            d = HeaderDataUnit(header_dict, data_arr)
            flux = data_arr
            self.pp.h2("Generating baseline IDs")
            bls, ant_arr = self.generateBaselineIds(d.n_ant)
            bl_lower = []
            while len(bl_lower) < len(flux):
                bl_lower += bls
        else:
            self.pp.h2("Inspecting visibility data")
            d = dada.DadaReader(self.filename, n_int, inspectOnly=True)
            self.dada_header = d.header

        # Gather the metadata
        metadata = {'stokes': ['XX', 'YY', 'XY', 'YX'], 'correlator': d.header["INSTRUMENT"],
                    'instrument': d.header["INSTRUMENT"], 'telescope': d.header["TELESCOPE"], 'tInt': d.t_int}

        ## Basic setup
        ## Time offset
        dt_obj = datetime.strptime(d.header["UTC_START"], "%Y-%m-%d-%H:%M:%S")
        time_offset = d.t_offset
        dt_obj = dt_obj + timedelta(seconds=time_offset)
        #date_obs = dt_obj.strftime("%Y-%m-%dT%H:%M:%S")
        #dd_obs = dt_obj.strftime("%Y-%m-%d")
        metadata['tstart'] = float(dt_obj.strftime("%s.%f"))

        ## Frequency information
        metadata['nchan'] = d.n_chans
        metadata['reffreq'] = d.c_freq_mhz * 1e6
        metadata['chanbw'] = d.chan_bw_mhz * 1e6
        metadata['refpixel'] = d.n_chans / 2 + 1

        return metadata

    def _compute_lst_ha(self, src):
        """ Helper function for computing LST, HA, and RA from timestamp

        returns (ra_deg, dec_deg, lst_deg, ha_deg)
        """
        # First, compute LST
        tt_source = coords.parse_timestring(self.date_obs)
        ts_source = calendar.timegm(tt_source)
        lst_deg = self.computeSiderealTime(ts_source)

        # Find HA and DEC of source
        if src.upper() == 'ZEN':
            H, d     = 0, float(self.site.lat)
            dec_deg  = np.rad2deg(d)
            ra_deg   = lst_deg
            ha_deg   = 0
        else:
            try:
                src_names = ledafits_config.src_names
                src_ras   = ledafits_config.src_ras
                src_decs  = ledafits_config.src_decs
                idx = src_names.index(src.upper())
                self.pp.h2("Phasing to %s" % src_names[idx])
                ra_deg, dec_deg = src_ras[idx], src_decs[idx]

                # Now we have the RA and DEC, need to convert into hour angle
                ha_deg = lst_deg - ra_deg
            except ValueError:
                raise ValueError("Cannot phase to unknown source '%s'" % src)
        return ra_deg, dec_deg, lst_deg, ha_deg

    def setDefaultsLeda(self, n_uv_rows=None):
        """ set LEDA specific default values """
        if n_uv_rows is None:
            n_uv_rows = len(self.d_uv_data["BASELINE"])
            
        self.setDefaults(n_uv_rows)

        self.d_frequency["CH_WIDTH"] = self.s2arr(ledafits_config.CH_WIDTH)
        self.d_frequency["TOTAL_BANDWIDTH"] = self.s2arr(ledafits_config.SUB_BW)
        self.h_uv_data["TELESCOP"]          = ledafits_config.TELESCOP
        self.h_array_geometry["ARRNAM"]     = ledafits_config.ARRNAM

    def loadAntArr(self):
        """ Loads ANTENNA and ARRAY_GEOMETRY tables as set in leda_config """
        self.pp.h1("Loading ANTENNA and ARRAY_GEOMETRY from JSON")

        if self.telescope in ('LEDA', 'LWAOVRO', 'LWA-OVRO', 'LEDAOVRO', 'LEDA512', 'LEDA-OVRO'):
            self.pp.h3("Data appears to be from LWAOVRO")

            self.site = ledafits_config.ovro
            self.h_array_geometry = load_json(ledafits_config.json_h_array_geometry)
            self.d_array_geometry = load_json(ledafits_config.json_d_array_geometry)
            self.h_antenna = load_json(ledafits_config.json_h_antenna)
            self.d_antenna = load_json(ledafits_config.json_d_antenna)
            self.z_elength = load_json(ledafits_config.json_antenna_el_lens)

        elif self.telescope in ('LWA1', 'LWA-1', 'LWA-NM', 'LWANM', 'LEDA64', 'LEDA64-NM'):
            self.pp.h3("Data appears to be from LWA1")

            self.site = ledafits_config.lwa1
            self.h_array_geometry = load_json(ledafits_config.json_h_array_geometry_nm)
            self.d_array_geometry = load_json(ledafits_config.json_d_array_geometry_nm)
            self.h_antenna = load_json(ledafits_config.json_h_antenna_nm)
            self.d_antenna = load_json(ledafits_config.json_d_antenna_nm)
            self.z_elength = load_json(ledafits_config.json_antenna_el_lens_nm)
        else:
            raise ValueError("Cannot load array geometry for %s" % self.telescope)

    def computeSiderealTime(self, ts=None):
        """ Computes the LST for a given timestamp.

        ts (float): Timestamp to use. If none is given, DATE-OBS value
                    is used in lieu.

        returns LST in degrees
        """

        self.pp.h2("Computing LST from UTC")
        if ts is not None:
            dt_utc = datetime.utcfromtimestamp(ts)
        else:
            dt_utc = datetime.strptime(self.date_obs, "%Y-%m-%dT%H:%M:%S")

        site = self.site
        site.date = dt_utc
        lst, lst_deg = site.sidereal_time(), site.sidereal_time() / 2 / np.pi * 360
        self.pp.pp("UTC: %s" % dt_utc)
        self.pp.pp("LST: %s (%s)" % (lst, lst_deg))
        return lst_deg

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

        self.pp.h1("Generating UVW coordinates")
        ra_deg, dec_deg, lst_deg, ha_deg = self._compute_lst_ha(src)
        H = np.deg2rad(ha_deg)
        d = np.deg2rad(dec_deg)

        self.pp.pp("LST:        %2.3f deg" % lst_deg)
        self.pp.pp("Source RA:  %2.3f deg" % ra_deg)
        self.pp.pp("Source DEC: %2.3f deg" % dec_deg)
        self.pp.pp("HA:         %2.3f deg" % np.rad2deg(H))

        try:
            assert H < 2 * np.pi and d < 2 * np.pi
        except AssertionError:
            raise ValueError("HA and DEC are too large (may not be in radians).")

        # Recreate list of baselines
        self.pp.h2("Computing UVW coordinates for %s" % src)
        xyz = self.d_array_geometry['STABXYZ']
        if 257 in set(self.d_uv_data["BASELINE"]):
            bl_ids, ant_arr = coords.generateBaselineIds(self.n_ant)
            bl_vecs = coords.computeBaselineVectors(xyz)
        else:
            bl_ids, ant_arr = coords.generateBaselineIds(self.n_ant, autocorrs=False)
            bl_vecs = coords.computeBaselineVectors(xyz, autocorrs=False)

        n_iters = int(len(self.d_uv_data["BASELINE"]) / len(bl_ids))

        self.pp.h2("Generating timestamps")
        dd, tt = [], []
        for ii in range(n_iters):
            jd, jt = coords.convertToJulianTuple(self.date_obs)
            tdelta = self.t_int * ii / 86400.0  # In days
            jds = [jd for jj in range(len(ant_arr))]
            jts = [jt + tdelta for jj in range(len(ant_arr))]
            dd.append(jds)
            tt.append(jts)

        self.d_uv_data["DATE"] = np.array(dd, dtype='float64').ravel()
        self.d_uv_data["TIME"] = np.array(tt, dtype='float64').ravel()

        if use_stored:
            self.pp.h2("Loading stored values")
            self.loadUVW()
        else:

            uvw = coords.computeUVW(bl_vecs, H, d, conjugate=conjugate)

            # Fill with data
            # TODO: update this so that it can lock to zenith or phase to src
            uu, vv, ww = [], [], []

            for ii in range(n_iters):
                uu.append(uvw[:, 0])
                vv.append(uvw[:, 1])
                ww.append(uvw[:, 2])

            self.d_uv_data["UU"] = np.array(uu).ravel()
            self.d_uv_data["VV"] = np.array(vv).ravel()
            self.d_uv_data["WW"] = np.array(ww).ravel()

        if update_src:
            self.pp.h2("Updating SOURCE table")
            self.d_source["SOURCE"] = self.s2arr(src)
            self.d_source["RAEPO"] = self.s2arr(ra_deg)
            self.d_source["DECEPO"] = self.s2arr(dec_deg)
            self.source = src

    def dumpUVW(self, filename):
        """ Dump precomputed UVW coordinates to file 
        
        filename (str): name of output file (.json format)
        """
        d = {"UU": self.d_uv_data["UU"], "VV": self.d_uv_data["VV"], "WW": self.d_uv_data["WW"],
             "BASELINE": self.d_uv_data["BASELINE"]}

        self.pp.h2("Dumping UVW coords to %s" % filename)
        dump_json(d, filename)

    def loadUVW(self, filename=None):
        """ Load precomputed UVW coordinates from file
        
        filename (str): name of input file. If not set, uses default
                        from leda_config file.
        """

        self.pp.h2("Loading UVW coordinates from file")
        # TODO: Make this check between 64 and 512!
        if not filename:
            filename = ledafits_config.json_uvw_coordinates
        d = load_json(filename)
        self.d_uv_data["UU"] = d["UU"]
        self.d_uv_data["VV"] = d["VV"]
        self.d_uv_data["WW"] = d["WW"]
        self.d_uv_data["BASELINE"] = d["BASELINE"]

    def leda_set_value(self, key, value):
        """ Set values which are commonly incorrect from uvfits writer """

        if key == 'ARRNAM':
            self.h_array_geometry['ARRNAM'] = value
        if key == 'INTTIM':
            self.d_uv_data['INTTIM'][:] = value
        if key == 'TELESCOP':
            self.h_uv_data['TELESCOP'] = value

    def remove_miriad_baselines(self):
        """ Remove baseline data for all antennas with IDs > 255

        Miriad-type UVFITS files use the convention
            ant1*256+ant2 if ants < 255
            ant1*2048+ant2+65536 if ants >255
        The miriad convention screws up import into many reduction packages.
        """

        self.pp.h1("Removing MIRIAD baselines")
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
                self.pp.err(k)
                self.pp.err(self.d_uv_data[k])
                raise
            except ValueError:
                self.pp.err(k)
                self.pp.err(len(self.d_uv_data[k]))
                self.pp.err(len(ok_bls))
                raise

        #for k in self.d_antenna.keys():
        #    self.d_antenna[k] = self.d_antenna[k][0:self.n_ant]
        #    #    print len(self.d_antenna[k])

        #for k in self.d_array_geometry.keys():
        #    self.d_array_geometry[k] = self.d_array_geometry[k][0::self.n_ant]
        #    #    print len(self.d_array_geometry[k])

        self.pp.h2("Fixing NOPCAL (setting to zero)")
        self.h_antenna["NOPCAL"] = 0

        self.pp.h2("Setting INTTIME to %s" % self.t_int)
        self.d_uv_data["INTTIM"] = np.ones_like(self.d_uv_data["INTTIM"]) * self.t_int

    def flag_antenna(self, antenna_id, reason=None, severity=0):
        """ Flag antenna as bad

        antenna_id (int): ID of antenna to flag. Starts at 1.
        reason (str): Defaults to None; short (<24 char) reason for flagging
        severity (int): -1 Not assigned, 0 Known bad, 1 Probably bad, 2 Maybe bad
        """
        self.pp.h2("Flagging antenna %i" % antenna_id)

        flag_keywords = ['SOURCE_ID', 'ARRAY', 'ANTS', 'FREQID', 'BANDS', 'CHANS', 'PFLAGS', 'REASON', 'SEVERITY']
        try:
            self.d_flag["SOURCE_ID"]
        except KeyError:
            for k in flag_keywords:
                self.d_flag[k] = []

        if reason is None:
            reason = "Known bad antenna."

        flag_k_zeros = ['SOURCE_ID', 'ARRAY', 'FREQID']
        flag_k_one = ['BANDS']
        for k in flag_k_one:
            self.d_flag[k].append(1)
        for k in flag_k_zeros:
            self.d_flag[k].append(0)

        self.d_flag["ANTS"].append((antenna_id, 0))
        self.d_flag["REASON"].append(reason)
        self.d_flag["PFLAGS"].append((1, 1, 1, 1))
        self.d_flag["SEVERITY"].append(severity)
        self.d_flag["CHANS"].append((0, 4096))

    def phase_to_src(self, src='ZEN', generate_uvw=True):
        """ Apply phase corrections to phase to source.

        Generates new UVW coordinates, then applies geometric delay (W component)
        to phase flux data to the new phase center.

        Parameters
        ----------
        src (str): Source to phase to. Sources are three capital letters:
            ZEN: Zenith (RA will be computed from timestamps)
            CYG: Cygnus A
            CAS: Cassiopeia A
            TAU: Taurus A
            VIR: Virgo A
        generate_uvw (bool): Skip regeneration of UVW coords?

        """
        self.pp.h1("Phasing flux data to %s" % src)

        current_tgs = self.d_uv_data["WW"]
        if generate_uvw is True:
            self.generateUVW(src, update_src=True)
        freqs = self.formatFreqs()
        w = 2 * np.pi * freqs  # Angular freq
        # Note WW *is* the geometric delay tg
        new_tgs = self.d_uv_data["WW"]

        try:
            assert self.d_uv_data["FLUX"].dtype == 'float32'
        except AssertionError:
            raise RuntimeError("Unexpected data type for FLUX: %s" % str(self.d_uv_data["FLUX"].dtype))
        flux = self.d_uv_data["FLUX"].view('complex64')

        bls = set(self.d_uv_data["BASELINE"])
        if not 257 in bls:
            bls, ant_arr = coords.generateBaselineIds(self.n_ant, autocorrs=False)
        else:
            bls, ant_arr = coords.generateBaselineIds(self.n_ant, autocorrs=True)
        n_int = len(flux) / len(bls)

        for nn in range(n_int):
            for ii in range(len(bls)):
                # Compute phases for X and Y pol on antennas A and B
                tg = new_tgs[nn * len(bls) + ii] - current_tgs[nn * len(bls) + ii]
                #if ant1 < ant2:
                #    tg *= -1    # Compensate for geometry
                p = np.exp(-1j * w * tg)  # Needs to be -ve as compensating delay
                phase_corrs = np.column_stack((p, p, p, p)).flatten()
                flux[nn * len(bls) + ii] = flux[nn * len(bls) + ii] * phase_corrs

            # Now we have applied geometric delays, we need to
            # convert from viewing as complex to viewing as floats
            assert flux.dtype == 'complex64'
            self.d_uv_data["FLUX"] = flux.view('float32')

    def apply_cable_delays(self, debug=True):
        """ Apply antenna cable delays

        Each cable introduces a phase shift of
            phi = 2 pi f t
        Visibility is VpVq*, so we need to apply
            exp(-i  (phip - phiq))
        to compensate for cable delay
        """

        self.pp.h1("Applying cable delays")
        #t0 = time.time()
        # Load antenna Electrical Lengths
        sol = ledafits_config.SPEED_OF_LIGHT
        try:
            els = self.z_elength["EL"]
        except:
            raise RuntimeError("No cable delay data for telescope '%s'" % self.telescope)
        els = np.array(els)
        tdelts = els / sol

        if debug:
            print "X-POL (ns)  \tY-POL (ns)"
            for line in tdelts:
                print "%2.2f   \t%2.2f" % (line[0] * 1e9, line[1] * 1e9)

        # Generate frequency array from metadata
        freqs = self.formatFreqs()
        # Compute phase delay for each antenna pair
        try:
            assert self.d_uv_data["FLUX"].dtype == 'float32'
        except AssertionError:
            raise RuntimeError("Unexpected data type for FLUX: %s" % str(self.d_uv_data["FLUX"].dtype))

        # Convert the data to complex values
        flux = self.d_uv_data["FLUX"].view('complex64')

        # Pre-compute the phasing information
        bls, ant_arr = coords.generateBaselineIds(self.n_ant)
        w = 2 * np.pi * freqs  # Angular freq
        delay_corrs = np.zeros((4, len(bls), len(freqs)), dtype=flux.dtype)
        for ii in range(len(bls)):
            ant1, ant2 = ant_arr[ii]
            bl = bls[ii]
            td1, td2 = tdelts[ant1 - 1, :], tdelts[ant2 - 1, :]

            # Compute phases for X and Y pol on antennas A and B
            pxa, pya, pxb, pyb = w * td1[0], w * td1[1], w * td2[0], w * td2[1]

            # Corrections require negative sign (otherwise reapplying delays)
            delay_corrs[0, ii, :] = np.exp(1j * (pxa - pxb))  # XX
            delay_corrs[1, ii, :] = np.exp(1j * (pya - pyb))  # YY
            delay_corrs[2, ii, :] = np.exp(1j * (pxa - pyb))  # XY
            delay_corrs[3, ii, :] = np.exp(1j * (pya - pxb))  # YX

        n_int = len(flux) / len(bls)
        for nn in range(n_int):
            for ii in range(len(bls)):
                e_xx = delay_corrs[0, ii, :].flatten()
                e_yy = delay_corrs[1, ii, :].flatten()
                e_xy = delay_corrs[2, ii, :].flatten()
                e_yx = delay_corrs[3, ii, :].flatten()

                phase_corrs = np.column_stack((e_xx, e_yy, e_xy, e_yx)).flatten()
                flux[nn * len(bls) + ii] = flux[nn * len(bls) + ii] * phase_corrs

        assert flux.dtype == 'complex64'
        self.d_uv_data["FLUX"] = flux.view('float32')
        
    def extractTotalPower(self, antenna_id, timestamps=False):
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
