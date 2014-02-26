#! /usr/bin/env python
# encoding: utf-8
"""
generate_from_dada
==================

Generate FITS-IDI files from dada. Test. 
"""

from test_main import *
import pylab as plt

if __name__ == '__main__':
    #uvw = LedaFits('uvf-zen.fitsidi')
    uvw = LedaFits('nm-zen.fitsidi')
    uvw.leda_set_value("ARRNAM", "LWA1")
    uvw.leda_set_value("TELESCOP", "LEDA64NM")
    p_vec = uvw._compute_pointing_vec("CYG")
     
    # Extract data for baseline 1-13
    for aa in range(32):
        ant1, ant2 = 13, aa + 1
        if aa < ant1:
            ws = -1
            bl_id = ant2 * 256 + ant1
        else:
            ws = 1
            bl_id = ant1 * 256 + ant2
        print bl_id
        bls = uvw.d_uv_data["BASELINE"]
        
        ant_locs = uvw.d_array_geometry["STABXYZ"]
        freqs = uvw.formatFreqs()
        w  = 2 * np.pi * freqs # Angular freq
        d = uvw.formatStokes()[0, bls == bl_id][0]
        
        # Compute compensating phase
        bl_vec   = ant_locs[ant1-1] - ant_locs[ant2-1]
        tg = np.dot(bl_vec, p_vec) / ledafits_config.SPEED_OF_LIGHT
    
        print "\nBL: %i Ants: %i, %i"%(bl_id, ant1, ant2)
        print "tg: %2.2f, bl_vec: %s"%(tg * 1e9, str(bl_vec))
        phs = np.exp(1j * w * tg * ws)
        plt.figure(1)
        plt.subplot(8, 4, aa+1)
        plt.plot(np.angle(phs))
        plt.plot(np.angle(d))
        plt.title("%i %i"%(ant1, ant2))
        
        plt.figure(2)
        plt.subplot(8, 4, aa+1)
        phs = np.exp(-1j * w * tg * ws)
        plt.plot(np.angle(phs * d))
        plt.title("%i %i"%(ant1, ant2))
    plt.show()
    
    