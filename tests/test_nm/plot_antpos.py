#!/usr/bin/env python
"""
plot_antpos.py -- Plot antenna positions (and UVW coordinates)

Plots the antenna positions for a given Interfits compatible file
(UVFITS, FITS-IDI, etc). Can also plot the UVW coordinates in the 
file.

Antenna positions are stored in FITS files under STABXYZ, and are
in geocentric coordinates (not lat-long-elev = topocentric). UVW
coordinates are also stored in the file, under UU, VV, and WW.

"""

import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import sys
from optparse import OptionParser

from test_main import *

def plot_ants(x,y,z, fig, c='#cc0000', m='o', three_d=False):
    
    xl, yl, zl, ='X', 'Y', 'Z'
    
    if three_d:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c=c)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_zlabel(zl)
    
    else:
        ax1 = fig.add_subplot(211)
        ax1.plot(x, y, marker=m, lw=0, c=c)
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        
        ax2 = fig.add_subplot(223)
        ax2.plot(x, z, marker=m, lw=0, c=c)
        ax2.set_xlabel("X")
        ax2.set_ylabel("Z")
        
        ax3 = fig.add_subplot(224)
        ax3.plot(y, z, marker=m, lw=0, c=c)
        ax3.set_xlabel("Y")
        ax3.set_ylabel("Z")         

def plot_uvw(x,y,z, fig, c='#cc0000'):
    
    ax1 = fig.add_subplot(211)
    ax1.plot(u, v, marker='o', lw=0, c=c)
    ax1.set_xlabel("U")
    ax1.set_ylabel("V")
    
    ax2 = fig.add_subplot(223)
    ax2.plot(u, w, marker='o', lw=0, c=c)
    ax2.set_xlabel("U")
    ax2.set_ylabel("W")

    ax3 = fig.add_subplot(224)
    ax3.plot(v, w, marker='o', lw=0, c=c)
    ax3.set_xlabel("V")
    ax3.set_ylabel("W")        

if __name__ == '__main__':

    p = OptionParser()
    p.set_usage('dada2uvfits.py <filename_in> [options]')
    p.set_description(__doc__)
    p.add_option("-u", "--uvw", dest="plot_uvw", action='store_true', 
                 help="Plots UVW positions in addition to antenna positions.")
    p.add_option("-p", "--print", dest="print_antpos", action='store_true', 
                 help="Print STABXYZ to screen.")
                                  
    (options, args) = p.parse_args(sys.argv[1:])
        
    c = [
            [0.7, 0, 0, 0.6],
            [0, 0.7, 0, 0.6],
            [0, 0, 0.7, 0.6]
        ]
    
    m = ['s', 'd', 'o']
    ii = 0
    fig = plt.figure(1)
    
    if options.plot_uvw:
        fig2 = plt.figure(2)
    
    for f in args:
        try:
            l = LedaFits(f)
            #l.phase_to_src('CYG')
            xyz = l.d_array_geometry["STABXYZ"]
            x, y, z = np.split(xyz, 3, axis=1)
            
            print f
            print "------------------"
            if options.print_antpos:
                print xyz
            
            plot_ants(x, y, z, fig, c=c[ii], m=m[ii])
            
            if options.plot_uvw:
                uvd = l.d_uv_data
                u, v, w = uvd["UU"], uvd["VV"], uvd["WW"]                
                u, v, w = u * 1e9, v * 1e9, w * 1e9
                
                plot_uvw(u, v, w, fig2, c=c[ii])
            
            ii +=1 
                
        except:
            raise
            
    plt.show()