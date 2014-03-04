InterFits
=========

Interfits is a python package for converting and reading radio interferometric data. It can read
UV-FITS and FITS-IDI, and also reads two FITS-IDI inspired data formats: JSON-IDI and HDF5-IDI.
InterFits can export data tables in JSON-IDI, FITS-IDI or HDF5-IDI.

Additionally, the LedaFits class reads the raw file formats used in the LEDA project (http://www.cfa.harvard.edu/LEDA/),
and contains methods for combining and computing observational metadata (e.g. LST and UVW coordinate calculation).

InterFits comes with a basic data viewer, called QtUv. This can open and display the raw data for any data type that
can be read by InterFits.

Requirements
------------

* Numpy > 1.6  - http://www.numpy.org/
* pyfits > 2.3 - http://www.stsci.edu/institute/software_hardware/pyfits
* lxml - http://lxml.de/
* h5py - http://www.h5py.org/
* pyephem - http://rhodesmill.org/pyephem/

For the included scripts, PyQt4 or PySide, pandas, ujson, LSL, and matplotlib are also required.


