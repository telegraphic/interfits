#! /usr/bin/env python
# encoding: utf-8
"""
uvfits_reader.py
================

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
    

class UvFits(object):
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
        
        # Should have 7 axes, likely
        if len(s) != 7:
            if len(s) == 6:
                print "Reshaping uv-data..."
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
        
        lower_t, upper_t = True, True
        for i in range(n_dumps):
            if bls[i*n_bls:(i+1)*n_bls] == bl_lower:
                #print "LOWER"
                upper_t = False
            if bls[i*n_bls:(i+1)*n_bls] == bl_upper:
                #print "UPPER"
                lower_t = False
            
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
    
    def plot_single_baseline(self, ant1, ant2, axis=0):
        """ Plot single baseline 
        
        ant1: int
            antenna number of first antenna
        ant2: int
            antenna number of second antenna
        """
        bls = self.d_uv_data['BASELINE']
        
        if ant1 > ant2: bl_id = 256*ant2 + ant1
        else: bl_id = 256*ant1 + ant2
        
        x_data    = self.d_uv_data['DATA'][bls == bl_id,0,0,0,:,axis]  # Baselines, freq and stokes
        x         = x_data[:,:,0] + 1j * x_data[:,:,1]
        
        fig = plt.figure(figsize=(4*3.5,3*3.2))
        figtitle = '%s %s: %s -- %s\n'%(self.telescope, self.instrument, self.source, self.date_obs)
        figtitle += "Baseline %i %i"%(ant1, ant2)
        fig.suptitle(figtitle, fontsize=18)
        ax = plt.subplot(121)
        plt.imshow(np.abs(x))
        plt.title("Amplitude")
        plt.xlabel("Frequency channel")
        plt.ylabel("Time")
        ax.set_aspect(x.shape[1] / x.shape[0])
        plt.colorbar(orientation='horizontal')
        
        ax = plt.subplot(122)
        plt.imshow(np.angle(x))
        plt.title("Phase")
        plt.xlabel("Frequency channel")
        plt.ylabel("Time")
        ax.set_aspect(x.shape[1] / x.shape[0])
        cbar = plt.colorbar(orientation='horizontal')
        cbar.set_ticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi-0.05])
        cbar.set_ticklabels(["$-\pi$","$-\pi/2$",0,"$\pi/2$","$\pi$"])
        
        plt.tight_layout()
        plt.show()
        
    
    def plot_visibilities(self, ref_ant=1, axis=0, plot_type='phase'):
        """ Plot visibilities
        
        ref_ant: int
            reference antenna
        axis: int
            which axis (stokes parameter) to view
        plot_type: str
            Plot type. Can be amplitude 'amp', phase 'phase'.
        """
        
        bls = self.d_uv_data['BASELINE']

        # Extract the relevant baselines using a truth array
        bls = bls.tolist()
        bl_ids = self.get_baseline_ids(ref_ant)
        bl_truths = np.array([(b in bl_ids) for b in bls])
        
        x_data    = self.d_uv_data['DATA'][bl_truths,0,0,0,:,axis]  # Baselines, freq and stokes
        x_cplx    = x_data[:,:,0] + 1j * x_data[:,:,1]
        
        
        # Plot the figure
        n_x, n_y = 4, 8
        fig = plt.figure(figsize=(4*3.5,3*3.2))
        figtitle = '%s %s: %s -- %s'%(self.telescope, self.instrument, self.source, self.date_obs)
        for i in range(n_x):
            for j in range(n_y):
                ax = fig.add_subplot(n_x, n_y, i*n_y + j +1)
                ax.set_title("%s %s"%(ref_ant, i*n_y + j +1))
                #ax.set_title("%s %s"%(i, j))
                x = x_cplx[i*n_y+j::self.n_ant]

                #img.set_interpolation('nearest') # Bicubic etc
                #img.set_cmap('jet')               # summer, hot, spectral, YlGnBu
                
                                
                if plot_type == 'phase':
                    img = ax.imshow(np.angle(x), vmin=-np.pi, vmax=np.pi)
                elif plot_type == 'amp':
                    img = ax.imshow(np.abs(x))
                else:
                    print "Error: plot_type %s not understood"%plot_type
                    raise

                if i == n_x-1:
                    ax.set_xlabel('Freq')
                if j == 0:
                    ax.set_ylabel('Time')
                ax.set_aspect(x.shape[1] / x.shape[0])
        
        if plot_type == 'phase':
            # Add phase colorbar
            cax = fig.add_axes([0.925,0.08,0.015,0.8])
            cbar = fig.colorbar(img, cax=cax)
            #cbar.set_label('Phase [rad]')
            cbar.set_ticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
            cbar.set_ticklabels(["$-\pi$","$-\pi/2$",0,"$\pi/2$","$\pi$"])
                    
            plt.subplots_adjust(left=0.05, right=0.9, top=0.9, bottom=0.05)
        else:
            plt.tight_layout()
        plt.suptitle(figtitle)
        plt.show()

    def plot_autocorrs(self, axis=0):
        """ Plot autocorrelations for all antennas
        """
        
        bls = self.d_uv_data['BASELINE']

        # Extract the relevant baselines using a truth array
        bls = bls.tolist()
        bl_ids = [256*i + i for i in range(1,33)]
        bl_truths = np.array([(b in bl_ids) for b in bls])
        
        print self.d_uv_data['DATA'].shape
        #x_data    = self.d_uv_data['DATA'][bl_truths,0,0,:,0,axis]  # Baselines, freq and stokes
        x_data    = self.d_uv_data['DATA'][bl_truths, 0,0,0,:,axis,:]
        x_cplx    = x_data[:,:,0] + 1j * x_data[:,:,1]
        
        print x_cplx.shape
        
        # Plot the figure
        n_x, n_y = 4, 8
        fig = plt.figure(figsize=(4*3.5,3*3.2))
        figtitle = '%s %s: %s -- %s'%(self.telescope, self.instrument, self.source, self.date_obs)
        for i in range(n_x):
            for j in range(n_y):
                ax = fig.add_subplot(n_x, n_y, i*n_y + j +1)
                ax.set_title("%s"%( i*n_y + j +1))
                #ax.set_title("%s %s"%(i, j))
                
                
                x = x_cplx[i*n_y+j::self.n_ant]
                
                x_pow     = np.abs(x)
                x_avg     = np.average(x_pow, axis=0)
                x_max     = np.max(x_pow, axis=0)
                x_med     = np.median(x_pow, axis=0)
                x_min     = np.min(x_pow, axis=0)
                
                #img.set_interpolation('nearest') # Bicubic etc
                #img.set_cmap('jet')               # summer, hot, spectral, YlGnBu
                
                ax.plot(x_avg)
                ax.plot(x_min)
                ax.plot(x_max)
                
                if i == n_x-1:
                    ax.set_xlabel('Freq')
                if j == 0:
                    ax.set_ylabel('Amplitude')
                #ax.set_aspect(x.shape[1] / x.shape[0])
        

        plt.tight_layout()
        plt.suptitle(figtitle)
        plt.show()



if __name__ == '__main__':
    
    filename = 'Test_0.uvfits'
    uv = UvFits(filename)
    
    uv.export_fitsidi('test.fitsidi', 'config/leda64.xml')
    
    uv.verify()
    
    #uv.plot_autocorrs(axis=1)
    #uv.plot_visibilities(ref_ant=5, axis=0, plot_type='phase')
    #uv.plot_single_baseline(1, 16)