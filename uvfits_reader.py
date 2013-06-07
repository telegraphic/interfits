# encoding: utf-8
"""
uvfits_reader.py
================

reads a uvfits file and loads it into InterFits class

"""

import sys, os, re
import pyfits as pf, numpy as np

import pylab as plt

class UvFits(object):
    def __init__(self, filename):
        self.filename = filename    
        self.read_uvfits()
        
    def searchKeys(self, pattern):
        """ Search through a header, returning a list of matching keys """
        
        matches =[]
        keys = self.uvhead.keys()
        for k in keys:
            m = re.search(pattern, k)
            if m: 
                matches.append(m.string)
        matches.sort()
        return matches
        
    def read_uvfits(self):
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
        
        self.fits     = pf.open(filename)
        self.uvdata   = self.fits[0].data
        self.uvhead   = self.fits[0].header
        
        self.antdata  = self.fits[1].data
        self.anthead  = self.fits[1].header
        
        # Load basic metadata
        self.telescope  = self.uvhead['TELESCOP'].strip()
        self.instrument = self.uvhead['INSTRUME'].strip()
        self.source     = self.uvhead['OBJECT'].strip()
        self.date_obs   = self.uvhead['DATE-OBS'].strip()
        
        self.n_ant    = self.antdata.shape[0]
        
        ctypes = self.searchKeys('CTYPE\d')
        caxes  = [int(ct.lstrip('CTYPE')) for ct in ctypes]
    
        print [(ct, self.uvhead[ct], self.uvhead['NAXIS%s'%ct.lstrip('CTYPE')]) for ct in ctypes]
    
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
    
    
    def plot_single_baseline(self, ant1, ant2, axis=0):
        """ Plot single baseline 
        
        ant1: int
            antenna number of first antenna
        ant2: int
            antenna number of second antenna
        """
        bls = self.uvdata['BASELINE']
        
        if ant1 > ant2: bl_id = 256*ant2 + ant1
        else: bl_id = 256*ant1 + ant2
        
        x_data    = self.uvdata['DATA'][bls == bl_id,0,0,:,axis]  # Baselines, freq and stokes
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
        
        bls = self.uvdata['BASELINE']

        # Extract the relevant baselines using a truth array
        bls = bls.tolist()
        bl_ids = self.get_baseline_ids(ref_ant)
        bl_truths = np.array([(b in bl_ids) for b in bls])
        
        # OLD METHOD THAT ONLY WORKS FOR  ANTENNA 1
        #bl_min, bl_max = (256*ref_ant, 256*ref_ant + self.n_ant)   
        #bl_truths = np.all((bls<=bl_max, bls>=bl_min), axis=0)
        
        x_data    = self.uvdata['DATA'][bl_truths,0,0,:,axis]  # Baselines, freq and stokes
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

if __name__ == '__main__':
    
    filename = 'band1.uvfits'
    uv = UvFits(filename)
    
    #uv.plot_visibilities(ref_ant=5, axis=0, plot_type='phase')
    uv.plot_single_baseline(1, 16)