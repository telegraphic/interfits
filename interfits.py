#! /usr/bin/env python
# encoding: utf-8
"""
interfits.py
============

reads a uvfits file and loads it into InterFits class

"""

import sys, os, re
import pyfits as pf, numpy as np

import pylab as plt

from pyFitsidi import *

def h1(headstr):
    """ Print a string as a header """
    print '\n', headstr
    underline = ''
    for i in range(len(headstr)): underline += '-'
    print underline

class VerificationError(Exception):
    """ Custom data verification exception """
    pass
    

class InterFits(object):
    def __init__(self, filename=None):
        self.filename = filename
        
        # Set up dictionaries to store data
        self.h_antenna = {}
        self.h_array_geometry = {}
        self.h_source = {}
        self.h_uv_data = {}
        self.h_frequncy = {}
        
        self.d_antenna = {}
        self.d_array_geometry = {}
        self.d_source = {}
        self.d_uv_data = {}
        self.d_frequency = {}

        self.stokes_codes = {
            1  : 'Stokes I',
            2  : 'Stokes Q',
            3  : 'Stokes U',
            4  : 'Stokes V',
            -1 : 'RR',
            -2 : 'LL',
            -3 : 'RL',
            -4 : 'LR',
            -5 : 'XX',
            -6 : 'YY',
            -7 : 'XY',
            -8 : 'YX'
        }
        
        # Check what kind of file to load
        if filename:
            regex = '([0-9A-Za-z-_]+).uvfits'
            match = re.search(regex, filename)
            if match: self.readFile('uvfits')
            
            regex = '([0-9A-Za-z-_]+).fitsidi'
            match = re.search(regex, filename)
            if match: self.readFile('fitsidi')
            
            regex = '([0-9A-Za-z-_]+).hdf'
            match = re.search(regex, filename)
            if match: self.readFile('hdf5')
    
    def readFile(self, filetype):
        """ Lookup dictionary (case statement) for file types """
        return {
             'uvfits'  : self.readUvfits,
             'fitsidi' : self.readFitsidi,
             'hdf5'    : self.readHdf5
              }.get(filetype)()
        
    def searchKeys(self, pattern, header):
        """ Search through a header, returning a list of matching keys 
        
        pattern: str
            regular expression pattern to match
        header: pyfits header
            header of a pyfits HDU to search within
        """
        
        matches =[]
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
        
        
        Numeric Codes for Stokes Parameters
    
        Code    Parameter
        ----    ---------
        1       I
        2       Q
        3       U
        4       V
        -1      RR
        -2      LL
        -3      RL
        -4      LR
        -5      XX
        -6      YY
        -7      XY
        -8      YX
        ----    ---------
    
        """
        
        h1("Opening uvfits data")
        
        self.fits     = pf.open(self.filename)
        self.uvdata   = self.fits[0].data
        self.uvhead   = self.fits[0].header
        
        self.antdata  = self.fits[1].data
        self.anthead  = self.fits[1].header
        
        # Load basic metadata
        self.telescope  = self.uvhead['TELESCOP'].strip()
        self.instrument = self.uvhead['INSTRUME'].strip()
        self.source     = self.uvhead['OBJECT'].strip()
        self.date_obs   = self.uvhead['DATE-OBS'].strip()
        self.n_ant      = self.antdata.shape[0]
        
        print self.fits.info()
        print "Telescope:  %s"%self.telescope
        print "Instrument: %s"%self.instrument
        print "Object:     %s"%self.source
        print "Date obs:   %s"%self.date_obs
        
        # Load array geometry data
        h1("Loading array geometry")
        ag_keywords = ['ARRAYX','ARRAYY','ARRAYZ']
        ag_data     = ['ANNAME', 'STABXYZ', 'NOSTA', 'MNTSTA', 'STAXOF']
        for k in ag_keywords: self.h_array_geometry[k] = self.anthead[k]
        for k in ag_data:     self.d_array_geometry[k] = self.antdata[k]
        
        # Load antenna table data
        h1("Loading antenna table")
        an_keywords = ['NOPCAL']
        an_data     = ['POLTYA','POLAA','POLCALA','POLTYB','POLAB','POLCALB']
        for k in an_keywords: self.h_antenna[k] = self.anthead[k]
        for k in an_data:     self.d_antenna[k] = self.antdata[k]
        
        # Load frequency table data
        # This is the first of the non-straightforward conversions
        # Need to find frequency index. It's probably CTYPE4 but safer to look it up
        h1("Loading frequency table")
        ctypes = self.searchKeys('CTYPE\d', self.uvhead)
        c = [(ct, self.uvhead[ct], self.uvhead['NAXIS%s'%ct.lstrip('CTYPE')]) for ct in ctypes]
        freq_cid = c[[x[1]=='FREQ' for x in c].index(True)][0].lstrip('CTYPE')
        
        ch_width = float(self.uvhead['CDELT%s'%freq_cid])
        bw       = ch_width * int(self.uvhead['NAXIS%s'%freq_cid])
        self.d_frequency['CH_WIDTH']        = ch_width
        self.d_frequency['TOTAL_BANDWIDTH'] = bw
        
        # Load source table data
        h1("Loading source table")
        self.d_source['SOURCE'] = self.source
        self.d_source['RAEPO']  = self.uvhead['OBSRA']
        self.d_source['DECEPO'] = self.uvhead['OBSDEC']
        
        
        # Load UV-DATA
        h1("Loading UV-data")
        uv_datacols = ['UU','VV','WW','BASELINE','DATE']
        for k in uv_datacols: self.d_uv_data[k] = self.uvdata[k]
        
        s = self.uvdata['DATA'].shape
        
        # Find stokes axis type and values
        for ct in ctypes:
            if self.uvhead[ct].strip() == 'STOKES':
                stokes_axid = int(ct.lstrip('CTYPE'))
                break
        
        stokes_axis_len = int(self.uvhead['NAXIS%i'%stokes_axid])
        stokes_code     = int(self.uvhead['CRVAL%i'%stokes_axid])
        stokes_delt     = int(self.uvhead['CDELT%i'%stokes_axid])
        stokes_vals  = range(stokes_code, stokes_code + stokes_delt* stokes_axis_len, stokes_delt)
        self.stokes_axis = [self.stokes_codes[i] for i in stokes_vals]
        
        # Should have 7 axes, likely
        if len(s) != 7:
            if len(s) == 6:
                #print "Reshaping uv-data..."
                new_shape = (s[0], s[1], s[2], 1, s[3], s[4],s[5])
                self.d_uv_data['DATA'] = self.uvdata['DATA'].reshape(new_shape)
            else:
                print "ERROR: Data axis not understood, incorrect length"
                print c
                raise
        else:
            self.d_uv_data['DATA'] = self.uvdata['DATA']
        #print c
    
    
    def readFitsidi(self):
        """ TODO """
        pass
    
    def readHdf5(self):
        """ TODO """
        pass
        
    def export_fitsidi(self, filename_out, config_xml):
        """ Export data as FITS IDI 
        
        filename_out: str
            output filename
        config_xml: str
            path to config file
        
        """
        
        h1('Creating Primary HDU')
        hdu = make_primary(config=config_xml)
        print hdu.header.ascardlist()
        print hdu
        
        h1('\nCreating ARRAY_GEOMETRY')
        tbl_array_geometry = make_array_geometry(config=config_xml, num_rows=32)
        print tbl_array_geometry.header

        h1('\nCreating ANTENNA')
        tbl_antenna = make_array_geometry(config=config_xml, num_rows=32)
        print tbl_antenna.header
        
        h1('\nCreating FREQUECY')
        tbl_frequency = make_frequency(config=config_xml, num_rows=32)
        print tbl_frequency.header
    
    def get_baseline_ids(self, ref_ant, triangle='upper'):
        """ Retrieve baseline ids """
        bl_ids = []
        
        if triangle =='upper':
            # Get values in array
            if ref_ant > 1:
                for i in range(1, ref_ant):
                    bl_ids.append(i*256 + ref_ant)
            
            # Get row
            bl_min, bl_max = (256*ref_ant, 256*ref_ant + self.n_ant)
            for b in range(bl_min, bl_max+1): bl_ids.append(b)
            
            return bl_ids
        else:
            print "Lower triangle not supported yet"
            raise     
    
    def verify_baseline_order(self):
        """ Check baseline IDs are in order """
        
        print "Verification: Checking uv_data baseline order..."
        
        bls =  [int(b) for b in self.d_uv_data['BASELINE']]
        
        # Generate lower an upper triangular matrices
        bl_lower, bl_upper = [], []
        for i in range(1, self.n_ant+1):
            for j in range(1, self.n_ant+1):
                if j >= i:
                    bl_lower.append(256*i + j)
                elif j <= i:
                    bl_upper.append(256*j + i)
        
        # Check every baseline is right, over all dumps
        n_bls    = self.n_ant*(self.n_ant-1)/2 + self.n_ant
        n_dumps = len(bls) / n_bls
        
        lower_t, upper_t = False, False
        for i in range(n_dumps):
            if bls[i*n_bls:(i+1)*n_bls] == bl_lower:
                #print "LOWER"
                upper_t = False
                lower_t = True
            if bls[i*n_bls:(i+1)*n_bls] == bl_upper:
                #print "UPPER"
                lower_t = False
                upper_t = True
            
            if not lower_t and not upper_t:
                raise VerificationError("Baseline order neither upper or lower triangular.")
        
        if lower_t: 
            print "Verification: OK. Baselines in lower triangular order."
        if upper_t:
            print "Verification: OK. Baselines in upper triangular order."
        
        return True       
        
    def verify(self):
        """ Run a series of diagnostics to test data validity """
        h1("Data verification")
        self.verify_baseline_order()


if __name__ == '__main__':
    
    #filename = 'Test_0.uvfits'
    filename = 'band1.uvfits'
    uv = InterFits(filename)
    
    uv.export_fitsidi('test.fitsidi', 'config/leda64.xml')
    
    #uv.verify()
