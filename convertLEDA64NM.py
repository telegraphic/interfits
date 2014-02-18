#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy
from datetime import datetime

from ledafits import LedaFits


def main(args):
	filenames = args
	for filename in filenames:
		# Load in the file
		uvw = LedaFits(filename)
		
		# Build the output name
		obsDate = datetime.strptime(uvw.date_obs, "%Y-%m-%dT%H:%M:%S")
		outname = "%s_%s_%s.FITS_1" % (uvw.instrument, uvw.telescope, obsDate.strftime("%Y%m%d%H%M%S"))
		print "Converting '%s' to '%s'..." % (filename, outname)
		
		# Add in the metadata
		uvw.loadAntArr()
		uvw.generateUVW(src='ZEN', use_stored=False, update_src=True)
		uvw.apply_cable_delays()
		#uvw.phase_to_src(src='ZEN', debug=False)
		uvw.exportFitsidi(outname)
		uvw.verify()


if __name__ == "__main__":
	main(sys.argv[1:])