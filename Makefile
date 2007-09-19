################################################################
# $Id: Makefile,v 1.2 2007/09/19 23:14:02 frolov Exp $
# DEFROST Makefile (Intel Fortran compiler)
################################################################

# Fortran compiler (adjust for your machine, -r8 and -fpp are mandatory)
FC = ifort
FFLAGS = -O3 -ipo -xT -r8 -pc80 -fpp -parallel
LDFLAGS = -static-intel

# Uncomment to use dynamically allocated arrays
# FFLAGS += -DDYNAMIC_ARRAYS

# FFTW library (set THREADS for threaded FFTW)
THREADS = 4
FFTWINC = 
FFTWLIB = -lfftw3

ifdef THREADS
FFLAGS += -DFFTWTHREADS=$(THREADS)
FFTWLIB += -lfftw3_threads
endif

# SILO library (native VisIt data format)
SILO = 1

ifdef SILO
SILOINC = 
SILOLIB = -L/usr/local/lib64 -lsilo -lhdf5

FFLAGS += -DSILO
endif


################################################################

all: defrost

clean:
	rm -f *.bak gmon.out core

defrost: defrost.f90
	$(FC) $(FFLAGS) $(FFTWINC) $(SILOINC) $^ -o $@ $(LDFLAGS) $(FFTWLIB) $(SILOLIB)
