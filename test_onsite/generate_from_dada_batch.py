#! /usr/bin/env python
# encoding: utf-8
"""
generate_from_dada
==================

Generate FITS-IDI files from dada. Test. 
"""

from test_main import *

def generate_fitsidi(filename_in, filename_out=None, src=None):
    """ Generate a fitsidi file from a dada file 
    
    filename_in (str): name of .dada file to open
    filename_out (str): output filename. Defaults to input with
                        .fitsidi extension appended
    """
    if filename_out is None:
        filename_out = os.path.split(filename_in)[0] + '.fitsidi'
    if src is None:
        src = 'ZEN'
    
    uvw = LedaFits(filename_in)
    uvw.loadAntArr()
    uvw.generateUVW(src=src, use_stored=False, update_src=True)
    uvw.remove_miriad_baselines()
    uvw.exportFitsidi(filename_out)
    
if __name__ == '__main__':
    import sys, os
    
    try:
        filename_in  = sys.argv[1]
    except:
        print "ERROR: you must enter a filename"
        print "USAGE: python generate_fitsidi.py <filename_in>"
        exit()
    
    src = 'CYG'
    
    for ii in range(1,12):
        for jj in (1, 2):
            try:
                fileroot = '/nfs/data/ledaovro%i/data%i/one/'%(ii, jj)
                filepath = os.path.join(fileroot, filename_in)
                filename_out = 'ledaovro%i_%i_%s.fitsidi'%(ii, jj, src)
                generate_fitsidi(filepath, filename_out, src=src)
            except:
                print "ERROR: Cannot process %s"%filepath

