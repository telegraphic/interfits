LEDA-512 DADA -> UVFITS CONVERTER
---------------------------------

dada2uvfits.py is a script to generate a uvfits file from LEDA-512 dada files.

USAGE: 
python dada2uvfits.py <filename> [options]
Options:
  -h, --help            show this help message and exit
  -l, --locktoinit      Locks phase centre to initial hour angle rather than
                        tracking in RA.
  -b BAND_ID, --bandid=BAND_ID
                        sub-band ID, defaults to 1 (starts at 1 not 0)
  -z, --test            Turn on test mode (do not run subprocesses)
  -F FIELD_NAME, --field=FIELD_NAME
                        Name of field. Default is Zenith.
  -t ACC_PER_FILE, --accperfile=ACC_PER_FILE
                        Number of accumulations per OUTPUT file.
  -T ACC_TO_READ, --acctoread=ACC_TO_READ
                        Number of accumulations to read from INPUT file.

Note: based upon corr2uvfits, written by Randall Wayth in 2008 for MWA.

REQUIREMENTS:
cfitsio v3: installed in a location where the compiler/linker will see it.
python modules: numpy, pyephem, colorama

CONTENTS:
SLALIB_C: 				SLA library. (not to be redistributed outside the MWA project)
uvifts.c uvfits.h: 		reader/writer of uvfits files
corr2uvfits.c: 			reads raw binary correlation files, reformats, calls uvfits writer
Makefile

BUILDING:
make the program by typing "make". This will create a standalone program called "corr2uvfits".
type ./corr2uvfits
for a summary of command line arguments.


