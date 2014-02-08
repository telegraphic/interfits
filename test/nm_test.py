#! /usr/bin/env python
# encoding: utf-8
"""
generate_from_dada
==================

Generate FITS-IDI files from dada. Test. 
"""

from test_main import *

def generate_fitsidi(filename_in, filename_out=None):
    """ Generate a fitsidi file from a dada file 
    
    filename_in (str): name of .dada file to open
    filename_out (str): output filename. Defaults to input with
                        .fitsidi extension appended
    """
    if filename_out is None:
        filename_out = os.path.split(filename_in)[0] + '.fitsidi'
    
    uvw = LedaFits(filename_in)
    uvw.loadAntArr()
    uvw.generateUVW(src='ZEN', use_stored=False, update_src=True)
    #uvw.flag_antenna(84)
    #uvw.apply_cable_delays()
    #uvw.remove_miriad_baselines()
    uvw.exportFitsidi(filename_out)
    
if __name__ == '__main__':
    import sys, os
    
    try:
        filename_in  = 'nm-2014-01-29-22.dada'
        filename_out = 'nm.fitsidi'
    except:
        print "ERROR: you must enter a filename"
        print "USAGE: python generate_fitsidi.py <filename_in> <filename_out>"
        exit()
    
    generate_fitsidi(filename_in, filename_out)

