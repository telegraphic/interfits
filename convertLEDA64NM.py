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
		obsFreq = int((uvw.formatFreqs()).mean() / 1e6)
		outname = "%s_%s_%s_%iMHz.FITS_1" % (uvw.instrument, uvw.telescope, obsDate.strftime("%Y%m%d%H%M%S"), obsFreq)
		print "Converting '%s' to '%s'..." % (os.path.basename(filename), outname)
		
		# Add in the UVW coordinates
		uvw.generateUVW(src='ZEN', use_stored=False, update_src=True)
		
		# Apply the cable delays
		uvw.apply_cable_delays()
		
		## Phase to zenith
		#uvw.phase_to_src(src='ZEN', debug=False)
		
		# Verify
		uvw.verify()
		
		# Save as FITS IDI
		uvw.exportFitsidi(outname)


if __name__ == "__main__":
	main(sys.argv[1:])