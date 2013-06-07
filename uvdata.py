# encoding: utf-8
"""
interfits.py
=========

Created by Danny Price on 2013-06-06.

Class for uv-data interchange, between uvfits, measurement sets, fits-idi and hdf

"""

import sys, os
import pyfits as pf, numpy as np

class InterFits(object):
    """ Class for uv data interchange """
    def __init__(self):
        uvdata = UvData()
    
    def __repr__(self):
        return "InterFits file"

class UvData(object):
    """ UV data table """
    def __init__(self):
        pass    

    def __repr__(self):
        return "UV Data"        