#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for converting LEDA64-NM Dada files into FITS IDI.  This script is 
also capable of grouping/combining files that were taken at the same time.
"""

import os
import re
import sys
import numpy
from datetime import datetime

from interfits.ledafits import LedaFits

_annameRE = re.compile('^.*?(?P<id>\d{1,3})$')


def main(args):
	filenames = args
	
	# Inspect the files to try and figure out what is what
	metadataList = []
	for filename in filenames:
		uvw = LedaFits()
		metadataList.append( (filename, uvw.inspectFile(filename)) )
		
	# Group the files by start time and save the filenames and frequency ranges
	groups = []
	for filename,metadata in metadataList:
		tStart = metadata['tstart']
		freqStart = metadata['reffreq'] + 0.0
		freqStop  = metadata['reffreq'] + metadata['nchan']*metadata['chanbw']
		
		## See if this file represents the start of a new group or not
		new = True
		for group in groups:
			if tStart == group[0]:
				new = False
				group[1].append(filename)
				group[2].append((freqStart,freqStop))
				break
				
		## A new group has been found
		if new:
			group = [tStart, [filename,], [(freqStart,freqStop),]]
			groups.append(group)
			
	# Report
	print "Got %i files with groupings:" % len(filenames)
	validity = []
	for i,group in enumerate(groups):
		## Sort the group by frequency
		freqs = []
		for start,stop in group[2]:
			freqs.append(start)
		freqOrder = [j[0] for j in sorted(enumerate(freqs), key=lambda x:x[1])]
		group[1] = [group[1][j] for j in freqOrder]
		group[2] = [group[2][j] for j in freqOrder]
		
		## Report and validate
		print "  Group #%i" % (i+1,)
		print "    -> start time %s (%.2f)" % (datetime.utcfromtimestamp(group[0]), group[0])
		valid = True
		for j,(name,(start,stop)) in enumerate(zip(group[1], group[2])):
			### Check for frequency continuity
			try:
				freqDiff = start - oldStop
			except NameError:
				freqDiff = 0
			oldStop = stop
			if freqDiff != 0:
				valid = False
				
			### Report on this file
			print "      %i: %s from %.2f to %.2f MHz" % (j+1, os.path.basename(name), start/1e6, stop/1e6)
		validity.append(valid)
		print "    -> valid set? %s" % valid
		
		## Reset the validity between groups
		del oldStop
	print " "
	
	# Combine
	for i,(valid,group) in enumerate(zip(validity,groups)):
		print "Combining group #%i..." % (i+1,)
		
		## Jump over invalid groups
		if not valid:
			print "  -> invalid, skipping"
			continue
			
		## Read in the files
		uvws = []
		for filename in group[1]:
			uvws.append( LedaFits(filename) )
			
		## Make a note of lowest frequency value
		freq_min = numpy.min(uvws[0].formatFreqs())
		
		## Concatenate together the various FLUX sets
		if len(uvws) > 1:
			timeBL = uvws[0].d_uv_data["FLUX"].shape[0]
			freqPolComp = uvws[0].d_uv_data["FLUX"].shape[1]
			new_uv_data = numpy.zeros((timeBL, freqPolComp*len(group[1])), dtype=uvws[0].d_uv_data["FLUX"].dtype)
			for i,uvw in enumerate(uvws):
				new_uv_data[:,i*freqPolComp:(i+1)*freqPolComp] = 1.0*uvw.d_uv_data["FLUX"]
			uvws[0].d_uv_data["FLUX"] = new_uv_data
			
		## Overwrite frequency axis keywords so that we can export UV_DATA table correctly
		uvws[0].h_common["REF_FREQ"] = freq_min
		uvws[0].h_common["REF_PIXL"] = 1
		uvws[0].h_common["NO_CHAN"]  *= len(group[1])
		uvws[0].h_params["NCHAN"] = uvws[0].h_common["NO_CHAN"] 
		uvws[0].d_frequency["TOTAL_BANDWIDTH"]  *= len(group[1])
		
		## Remove the other LedaFits instances since we only need the first one now
		while len(uvws) > 1:
			del uvws[-1]
			
		## Load in the current antenna mapping to associate antenna IDs with stand names
		ant_ids = uvws[0].d_array_geometry['NOSTA']
		ant_nms = [int(_annameRE.match(nm).group('id')) for nm in uvws[0].d_array_geometry['ANNAME']]
		
		## Figure out which antenna ID each switching outrigger is
		### Stand numbers for the switching outriggers
		outriggers = [35, 257, 259]
		
		ant_ids_outriggers = []
		for outrigger in outriggers:
			ant_ids_outriggers.append( ant_ids[ant_nms.index(outrigger)] )
			
		## Extract
		import pylab
		pols = uvws[0].stokes_axis
		freqs = uvws[0].formatFreqs()
		for stand,ant_id in zip(outriggers, ant_ids_outriggers):
			print "  -> Extracting data for Stand #%i..." % stand
			timestamps, data = uvws[0].extractTotalPower(ant_id, timestamps=True)
			tRel = timestamps - timestamps[0]
			tRel *= 24*3600	# Timestamps in seconds
			data = data.real	# Autocorrelations should be real
			
			### Plot mean power over time
			pylab.figure()
			for i in xrange(data.shape[0]):
				pylab.plot(tRel, data[i,:,:].mean(axis=-1), linestyle='-', marker='o', label='%i-%s' % (stand, pols[i]))
			pylab.xlabel('Time [s]')
			pylab.ylabel('Mean PSD [arb.]')
			pylab.legend()
			pylab.draw()
			
			### Plot the spectra
			pylab.figure()
			for i in xrange(data.shape[0]):
				pylab.subplot(2, 2, i+1)
				for j in xrange(data.shape[1]):
					pylab.plot(freqs/1e6, data[i,j,:])
				pylab.title("%i-%s" % (stand, pols[i]))
				pylab.xlabel('Frequency [MHz]')
				pylab.ylabel('PSD [arb.]')
			pylab.draw()
		pylab.show()


if __name__ == "__main__":
	main(sys.argv[1:])
