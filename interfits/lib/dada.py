# -*- coding: utf-8 -*-

"""
dada.py
=======

Module for reading in data in the Dada format.
"""

import numpy as np
import glob
import os

#from interfits.lib.timeit import timeit

__version__ = '0.0'
__all__ = ['DadaReader', 'lookup_warn', '__version__', '__all__']


def lookup_warn(table, key, default=None):
    try:
        return table[key]
    except KeyError:
        if default is not None:
            print "#Warning: No key '%s'; using default value of %s" \
                  % (key, default)
            return default
        else:
            print "#Warning: No key '%s'" % key
            return None


class DadaReader(object):
    """ Dada file reader for raw LEDA correlator data.

    Reads the header of a dada file and converts data into
    an ndimensional visibility matrix.

    Parameters
    ----------
    filename: str
        name of dada file to open
    n_int: int
        number of integrations to read. If None, will only read header
    inspectOnly: bool
        If inspectOnly, will only read header and will not unpack data.
    """
    DEFAULT_HEADER_SIZE = 4096

    def __init__(self, filename, n_int=None, inspectOnly=False):
        self.filename = filename
        self.read_header()

        # Load entire file unless n_int is specified
        if not inspectOnly:
            if n_int is None:
                file_size_hdr = int(self.header["FILE_SIZE"])
                file_size_dsk = os.path.getsize(filename)
                file_size     = min([file_size_hdr, file_size_dsk])
                bpa           = int(self.header["BYTES_PER_AVG"])
                
                n_int = file_size / bpa

            self.compute_matrix_indexes()
            self.read_data(0, n_int=n_int)
            self.datestamp = self.header['UTC_START']

            
    def __repr__(self):
       return str(self.header)

    def read_header(self, header_size=None, extension=None):
        """ Read dada file header """
        if header_size is None:
            header_size = self.DEFAULT_HEADER_SIZE
        if extension is None:
            extension = 'dada'
        byte_offset = 0
        #self.datestamp = datestamp
        self.extension = extension
        self.header_size = header_size
        filename = self.filename
        f = open(filename, 'rb')
        headerstr = f.read(header_size)
        self.parse_header(headerstr)
        f.seek(0, 2)
        self.bytes_per_file = f.tell() - self.header_size
        f.close()

    #@timeit
    def parse_header(self, headerstr):
        """ Parse dada header and form useful quantities """
        header = {}
        for line in headerstr.split('\n'):
            try:
                key, value = line.split()
            except ValueError:
                break
            key = key.strip()
            value = value.strip()
            header[key] = value
        self.header = header

        # Load commonly used values into variables
        self.c_freq_mhz     = float(lookup_warn(header, 'CFREQ', 0.))
        self.bandwidth_mhz  = float(lookup_warn(header, 'BW', 1.))
        self.n_chans = int(header["NCHAN"])
        self.n_pol   = int(header["NPOL"])
        self.n_ant   = int(header["NSTATION"])
        try:
            self.chan_bw_mhz    = float(header["CHAN_WIDTH"])
        except:
            self.chan_bw_mhz    = self.bandwidth_mhz / self.n_chans

        # Calculate number of integrations within this file
        # File may not be complete, hence file_size_dsk is read too
        file_size_hdr = int(header["FILE_SIZE"])
        file_size_dsk = os.path.getsize(self.filename)
        file_size     = min([file_size_hdr, file_size_dsk])
        bpa           = int(header["BYTES_PER_AVG"])
        self.n_int = file_size / bpa

        # Calculate integration time per accumulation
        tsamp      = float(header["TSAMP"]) * 1e-6   # Sampling time per channel, in microseconds
        navg       = int(header["NAVG"])             # Number of averages per integration
        int_tim    = tsamp * navg                    # Integration time is tsamp * navg
        self.t_int = int_tim

        # Calculate the time offset since the observation started
        byte_offset = int(header["OBS_OFFSET"])
        bytes_per_avg = int(header["BYTES_PER_AVG"])
        num_int_since_obs_start = byte_offset / bytes_per_avg
        time_offset_since_obs_start = num_int_since_obs_start * int_tim
        self.t_offset = time_offset_since_obs_start

        try:
            self.n_station = int(header['NSTAND'])
        except KeyError:
            self.n_station = int(lookup_warn(header, 'NSTATION', 32))
        self.n_input = self.n_station * self.n_pol
        self.n_dim = int(header['NDIM'])
        self.n_bit = int(header['NBIT'])
        self.dtype = \
            np.float32 if self.n_bit == 32 else \
                np.int16 if self.n_bit == 16 else \
                    np.int8 if self.n_bit == 8 else \
                        None
        self.navg = int(lookup_warn(header, 'NAVG', 25 * 8192))
        self.bytes_per_avg = int(lookup_warn(header, 'BYTES_PER_AVG', 10444800))
        self.data_order = lookup_warn(header, 'DATA_ORDER',
                                      'REG_TILE_TRIANGULAR_2x2')

    #@timeit
    def compute_matrix_indexes(self):
        """ Compute the matrix indexes """
        if self.data_order == 'REG_TILE_TRIANGULAR_2x2':
            reg_rows = 2
            reg_cols = 2
            self.matlen = self.reg_tile_triangular_matlen(reg_rows, reg_cols)
            # Build lookup table to map matrix idx --> row/col
            self.matrows = np.zeros(self.matlen, dtype=np.uint32)
            self.matcols = np.zeros(self.matlen, dtype=np.uint32)
            for i in xrange(self.matlen):
                row, col = self.reg_tile_triangular_coords(i, reg_rows, reg_cols)
                self.matrows[i] = row
                self.matcols[i] = col
        else:
            raise KeyError("Unsupported data order '%s'" % self.data_order)

    #@timeit
    def read_data(self, first_int, n_int=1):
        """
        Returns the specified integrations as a numpy array with shape:
        (nint, nchans, nstation, nstation, npol, npol), dtype=complex64

        Parameters
        ----------
        first_int: int
            Number of integrations to skip from start of file (i.e. offset)
        n_int: int
            Number of integrations to read
        """
        byte_offset = first_int * self.bytes_per_avg
        nbytes = n_int * self.bytes_per_avg

        file_idx = byte_offset // self.bytes_per_file
        file_offset = byte_offset % self.bytes_per_file
        file_byte_label = file_idx * self.bytes_per_file
        filename = self.filename

        print "#Reading   ", filename
        f = open(filename, 'rb')
        f.seek(self.header_size + file_offset)
        # Note: We load as raw bytes to allow arbitrary file boundaries
        data = np.fromfile(f, dtype=np.uint8, count=nbytes)
        f.close()

        # Transform data into a visibilty matrix
        data = self.transform_raw_data(data, n_int)
        return data

    #@timeit
    def transform_raw_data(self, data, n_int, fill_conjugate=False):
        """ Transform dada data into a useful visibility matrix

        Parameters
        ----------
        data: np.ndarray
            data array to transform into a visibility matrix
        n_int:
            number of integrations in data array
        fill_conjugate: bool
            Default False. Will compute upper half (conjugate) of full visibilty
            matrix if set to True.

        Notes
        -----
        This function dominates data processing time. Unless you need a full visibility
        matrix, do not run this with fill_conjugate=True, as this approximately doubles
        the amount of processing time.

        TODO: This may break if system endianness is different
        TODO: Add support for outputting in upper/lower triangular format
        """
        #print "#Converting", self.filename

        data = data.view(dtype=self.dtype).astype(np.float32)
        # Note: The real and imag components are stored separately
        data = data.reshape((n_int, 2, self.n_chans, self.matlen))

        # Scatter values into new full matrix
        fullmatrix = np.zeros((n_int, self.n_chans, self.n_input, self.n_input),
                              dtype=np.complex64)

        if not fill_conjugate:
            # Note cols then rows and conjugation (-1j) -- this is required
            data = data[..., 0, :, :] + np.complex64(-1j) * data[..., 1, :, :]
            fullmatrix[..., self.matcols, self.matrows] = data
        else:
            # Fill out the other (conjugate) triangle
            # Note rows then cols and no conjugation -- in contrast to above
            data = data[..., 0, :, :] + np.complex64(1j) * data[..., 1, :, :]
            fullmatrix[..., self.matrows, self.matcols] = data
            tri_inds = np.arange(self.n_input * (self.n_input + 1) / 2, dtype=np.uint32)
            rows, cols = self.triangular_coords(tri_inds)
            fullmatrix[..., cols, rows] = np.conj(fullmatrix[..., rows, cols])

        # Reorder so that pol products change fastest
        fullmatrix = fullmatrix.reshape(n_int, self.n_chans,
                                        self.n_station, self.n_pol,
                                        self.n_station, self.n_pol)
        fullmatrix = fullmatrix.transpose([0, 2, 4, 1, 3, 5])
        self.data = fullmatrix


    def triangular_coords(self, matrix_idx):
        """ Generate triangular coordinates. Returns row, col for each matrix index """
        row = (-0.5 + np.sqrt(0.25 + 2 * matrix_idx)).astype(np.uint32)
        col = matrix_idx - row * (row + 1) / 2
        return row, col

    def reg_tile_triangular_matlen(self, reg_rows, reg_cols):
        """ Compute the matrix length of regular tile triangle data """
        return (self.n_station / reg_rows + 1) * \
               (self.n_station / reg_cols / 2) * self.n_pol ** 2 * reg_rows * reg_cols

    def reg_tile_triangular_coords(self, matrix_idx, reg_rows, reg_cols):
        """ Compute row, col for a given matrix index, with given register tile sizes """
        npol = self.n_pol
        reg_tile_nbaseline = (self.n_station / reg_rows + 1) * (self.n_station / reg_cols / 2)
        rem = matrix_idx
        reg_col = rem / (reg_rows * reg_tile_nbaseline * npol * npol)
        rem %= (reg_rows * reg_tile_nbaseline * npol * npol)
        reg_row = rem / (reg_tile_nbaseline * npol * npol)
        rem %= (reg_tile_nbaseline * npol * npol)
        tile_row, tile_col = self.triangular_coords(rem / (npol * npol))
        rem %= (npol * npol)
        pol_col = rem / npol
        rem %= npol
        pol_row = rem

        row = pol_col + npol * (reg_row + reg_cols * tile_row)
        col = pol_row + npol * (reg_col + reg_rows * tile_col)

        return row, col
