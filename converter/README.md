Readme
======

**dada2uvfits.py** is a script to generate a *.uvfits* file from a LEDA-512 *.dada* file.
Briefly, this script:
* Converts from *.dada* to *.LA* and *.LC* using **lconvert**
* Generates a header file from the *.dada* filestamp
* Converts from *.LA* and *.LC* into *.uvfits*. This is done with **corr2uvfits**, 
  which was written by Randall Wayth (back in 2008 for MWA).

Usage
-----
From the command line, run:

```
python dada2uvfits.py <filename> [options]
```

Options are:

```python
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
```

Requirements
------------
**cfitsio v3:** installed in a location where the compiler/linker will see it.

**SLALIB_C:**  A directory of SLALIB which is not to be redistributed (and so not in this repo)

**python modules:** numpy, pyephem, colorama

Install
-------
cd into the converter directory and run

```
  make
```

Which will build **corr2uvfits** and **lconvert**


