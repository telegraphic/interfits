#! /usr/bin/env python
# encoding: utf-8
"""
generate_from_dada
==================

Generate FITS-IDI files from dada. Test. 
"""

from test_main import *

if __name__ == '__main__':
    uvw = LedaFits('nm-zen.fitsidi')
    uvw.generateUVW(src='CYG', use_stored=False, update_src=True)
    uvw.phase_to_src(src='CYG', debug=True)