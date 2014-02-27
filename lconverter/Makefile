# Makefile for corr2uvfits.

CFLAGS=-g -O -Wall -D_FILE_OFFSET_BITS=64 -L.
#CFITSIO_INCS=$(shell pkg-config --silence-errors --cflags cfitsio)
#CFITSIO_LIBS=$(shell pkg-config --silence-errors --libs cfitsio)
CFITSIO_INCS=-I../cfitsio/include
CFITSIO_LIBS=-L../cfitsio/lib
INCS=$(shell python -c "if len('${INCLUDE}')>0:print ' '.join(['-I ' + s for s in '${INCLUDE}'.split(':')])") -L${CFITSLIB} -I${CFITSINC}

TARGETS=corr2uvfits lconvert
SRCDIR=CONV2UVFITS


all: $(TARGETS)

corr2uvfits: corr2uvfits.c uvfits.c libsla.a 
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o corr2uvfits corr2uvfits.c uvfits.c $(CFITSIO_LIBS) -lcfitsio -lsla -lm
	
lconvert:  Lconverter512OV.c Lconverter64NM.c
	cc -o lconvert512OV Lconverter512OV.c
	cc -o lconvert64NM Lconverter64NM.c

libsla.a:
	cd SLALIB_C ; make
	rm -f SLALIB_C/*.o
clean:
	rm -f *.o $(TARGETS) libsla.a SLALIB_C/*.o 
	rm -rf corr2uvfits.dSYM
