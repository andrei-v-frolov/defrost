################################################################
# $Id: Makefile,v 1.1 2007/09/12 22:00:24 frolov Exp $
# DEFROST Makefile (Intel Fortran compiler)
################################################################

# Fortran compiler (adjust for your machine, -r8 and -fpp are mandatory)
FC = ifort
FFLAGS = -O3 -ipo -xT -r8 -pc80 -fpp -parallel
LDFLAGS = -static-intel

# FFTW library (set THREADS for threaded FFTW)
THREADS = 4
FFTWINC = 
FFTWLIB = -lfftw3

ifdef THREADS
FFLAGS += -DFFTWTHREADS=$(THREADS)
FFTWLIB += -lfftw3_threads
endif

# SILO library (native VisIt data format)
SILOINC = 
SILOLIB = -L/usr/local/lib64 -lsilo -lhdf5

FFLAGS += -DSILO


################################################################

all: defrost

clean:
	rm -f *.bak gmon.out core

defrost: defrost.f90
	$(FC) $(FFLAGS) $(FFTWINC) $(SILOINC) $^ -o $@ $(LDFLAGS) $(FFTWLIB) $(SILOLIB)
