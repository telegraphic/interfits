#! /usr/bin/env python
# encoding: utf-8
"""
dada2uvfits.py
Convert LEDA-512 OV data in dada format into uvfits format
created: 08 July 2013

lat: 34.070, lng: -107.628, elev: 2133.6 m 
"""

import re, sys, datetime, os, subprocess, shutil
from optparse import OptionParser
from colorama import Fore, Back, Style

from datetime import timedelta
import ephem
import numpy as np

from test_main import *


# Correlator setup - should not need to change
SUBBAND_BWIDTH  = 14.4     # Bandwidth in MHz
SUBBAND_NCHANS  = 600       # Number of channels in each subband
START_FREQ      = 51.6 - 14.4/2 # Start frequency of correlator in MHz
N_SCANS         = 100        # Number of scans in a dada file

(latitude, longitude, elevation) = ('34.07', '-107.628', 2133.6)

def h1(text):
    """ Print text as a fancy header """
    
    print(Fore.GREEN + '\n%s'%text)
    line = ""
    for ii in range(len(text)):
        line += "-"
    print(Fore.GREEN + line)
    print(Fore.WHITE)

def printRed(text):
    """ Print something in red """
    print(Fore.RED + text + Fore.WHITE)
    
def callSubprocess(args, test=False):
    """ Run a subprocess call with some printing and stuff """
    print(Fore.MAGENTA),
    for arg in args: 
        print arg,
    print(Fore.WHITE)
    if not test: 
        subprocess.call(args)

def computeLstFromFilename(filename):
    """ Print the sidereal times of a dada file. Now reads DADA header """
    # Create OVRO observer
    
    ov = ephem.Observer()
    ov.lon  = longitude
    ov.lat  = latitude
    ov.elev = elevation

    d = dada.DadaSubBand(filename, n_int=0)
    
    dt_obj = datetime.strptime(d.header["UTC_START"], "%Y-%m-%d-%H:%M:%S")                
    tsamp  = float(d.header["TSAMP"]) * 1e-6 # Sampling time per channel, in microseconds 
    navg   = int(d.header["NAVG"])           # Number of averages per integration         
    int_tim = tsamp * navg                   # Integration time is tsamp * navg           
                                                                                          
    byte_offset = int(d.header["OBS_OFFSET"])                                             
    bytes_per_avg = int(d.header["BYTES_PER_AVG"])                                        
    num_int = byte_offset / bytes_per_avg                                                 
    time_offset = num_int * int_tim                                                       
    
    pat = '(\d+)-(\d+)-(\d+)[-_](\d\d)[:h](\d\d)[:m](\d\d)$'

    match = re.search(pat, d.header["UTC_START"])
    if match:
        # Convert re match to integers, apart from file extension
        #(y, m, d, hh, mm, ss) = [int(m) for m in match.groups()[:-1]]
        dt = dt_obj + timedelta(seconds=time_offset)
        ov.date = dt
        lst = ov.sidereal_time()
        date_str = "%04d%02d%02d"%(dt.year,dt.month,dt.day)
        time_str = "%02d%02d%02d"%(dt.hour,dt.minute,dt.second)
        lst_str  = str(float(lst) / 2 / np.pi * 24)
        #print lst
        #print lst_str  
        #lst = str(lst).split(":")
        #lst_str  = "%s%s%s"%(lst[0], lst[1], lst[2].split(".")[0])
        
        printRed( "UTC START:   %s"%d.header["UTC_START"]         )
        printRed( "TIME OFFSET: %s"%timedelta(seconds=time_offset))
        printRed( "NEW START:   (%s, %s)"%(date_str, time_str)    )
        
        return date_str, time_str, lst_str
    else:
        print filename
        raise Exception("FilenameToSiderealError")

def findAndReplace(replace_dict, filename_in, filename_out, test_mode=False):
    """ Open a file, search for substitutions, output to a new file"""
    f = open(filename_in)
    lines = f.readlines()
    f.close()
    
    print "Generating %s"%filename_out
    for i in range(len(lines)):
        for k in replace_dict.keys():
            (lines[i], n_subs) = re.subn("\$"+k+"\$", str(replace_dict[k]), lines[i])
        
            if n_subs > 0:
                print " \'%s\' -> \'%s\'"%(k, replace_dict[k])
    
    if not test_mode:
        print "Writing to file %s"%filename_out           
        fo = open(filename_out, "w")
        fo.writelines(lines)
        fo.close()

def generateHeader(param_dict, filename_out,  test_mode=False, template="uvfits_headers/header.tpl"):
    """ Generate header for corr2uvfits """
    findAndReplace(param_dict, template ,filename_out, test_mode)

if __name__ == "__main__":
    # Option parsing to allow command line arguments to be parsed
    p = OptionParser()
    p.set_usage('dada2uvfits.py <filename_in> [options]')
    p.set_description(__doc__)
    p.add_option("-l", "--locktoinit", dest="lock_to_init", action='store_true', 
                 help="Locks phase centre to initial hour angle rather than tracking in RA.")
    p.add_option("-b", "--bandid", dest="band_id", type="int", default=1,
                 help="sub-band ID, defaults to 1 (starts at 1 not 0)")
    p.add_option("-z", "--test", dest="test", action='store_true', 
                 help="Turn on test mode (do not run subprocesses)")
    p.add_option("-F", "--field", dest="field_name", type="str", default="Zenith",
                 help="Name of field")
    p.add_option("-T", "--acctoread", dest="acc_to_read", type="int", default=99999,
                 help="Number of accumulations to read from INPUT file.")
                                  
    (options, args) = p.parse_args(sys.argv[1:])
    
    try:
        filename_dada = args[0]
        filename_la = filename_dada+'_%s.LA'%(options.band_id - 1)
        filename_lc = filename_dada+'_%s.LC'%(options.band_id - 1)
    except IndexError:
        print "Error: you must pass a filename."
        print "use -h for help"
        exit()
    
    h1("Generating headers")
    date_str, time_str, lst_str = computeLstFromFilename(filename_dada)
    fileroot_out = "%s_b%s_d%s_utc%s"%(options.field_name, options.band_id, date_str, time_str)
    if not os.path.exists(fileroot_out) and not options.test:
        try:
            os.mkdir(fileroot_out)
        except:
            raise
    
    # New method: Read from dada header
    d = dada.DadaSubBand(filename_dada, n_int=0)
    tsamp  = float(d.header["TSAMP"]) * 1e-6 # Sampling time per channel, in microseconds 
    navg   = int(d.header["NAVG"])           # Number of averages per integration         
    int_tim = tsamp * navg                   # Integration time is tsamp * navg
    byte_offset = int(d.header["FILE_SIZE"])                                             
    bytes_per_avg = int(d.header["BYTES_PER_AVG"])                                        
    num_int = byte_offset / bytes_per_avg                                                 
    time_offset = num_int * int_tim
    
    header_params = {
        'FIELDNAME'  : options.field_name,
        'N_CHANS'    : int(d.header["NCHAN"]),
        'N_SCANS'    : num_int,
        'INT_TIME'   : int_tim,
        'FREQCENT'   : float(d.header["CFREQ"]),
        'BANDWIDTH'  : float(d.header["BW"]),
        'RA_HRS'     : lst_str,
        'DEC_DEGS'   : latitude,
        'DATE'       : date_str,
        'TIME'       : time_str
    }
    
    generateHeader(header_params, os.path.join(fileroot_out, 'header.txt'),  options.test)
    
    ######
    # Run lconvert (currently doesn't support arguments!)
    ######
    
    h1("Running lconvert")
        
    lconvert_args = ['./lconvert64NM', 
        filename_dada
        ]
    
    callSubprocess(lconvert_args, options.test)
    
    ######
    # Run corr2uvfits on output files
    ######
    h1("Running corr2uvfits")    
    if options.lock_to_init: 
        corr2uvfits_flags = '-l'
    else: 
        corr2uvfits_flags = ''
    corr2uvfits_args = ['./corr2uvfits', 
        '-A', '%s,%s'%(longitude,latitude),
        '-a', filename_la,
        '-c', filename_lc,
        '-o', os.path.join(fileroot_out, fileroot_out+'.uvfits'), 
        '-H', os.path.join(fileroot_out, 'header.txt'),
        '-I', os.path.join('uvfits_headers', 'instr_config.tpl'),
        '-S', os.path.join('uvfits_headers', 'antenna_locations.tpl'), 
        '-f', str(1),
        corr2uvfits_flags
        ]
        
    callSubprocess(corr2uvfits_args, test=options.test)
    
    print "\nTidying up..."
    os.remove(filename_la)
    os.remove(filename_lc)
    print "DONE!"