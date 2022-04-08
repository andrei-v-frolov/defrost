################################################################
# $Id: Makefile,v 1.3 2012/04/20 23:19:39 frolov Stab $
# DEFROST Makefile (Intel Fortran compiler)
################################################################

# Fortran compiler (adjust for your machine, -r8 and -fpp are mandatory)
FC = ifort
FFLAGS = -xT -O3 -ipo -r8 -pc80 -fpp
LDFLAGS = -static-intel

# Uncomment for parallel build
THREADS = 4

# Uncomment to use dynamically allocated arrays
# FFLAGS += -DDYNAMIC_ARRAYS

# Changing memory model might be advisable instead
# FFLAGS += -mcmodel=large
# LDFLAGS = -shared-intel

# FFTW library (set THREADS for threaded FFTW)
FFTWINC = 
FFTWLIB = -lfftw3

ifdef THREADS
FFLAGS += -parallel -DTHREADS=$(THREADS)
FFTWLIB += -lfftw3_threads
endif

# SILO library (native VisIt data format)
SILO = 1

ifdef SILO
SILOINC = -I/usr/local/include
SILOLIB = -L/usr/local/lib64 -lsilo -lhdf5

FFLAGS += -DSILO
endif


################################################################

all: defrost

clean:
	rm -f *.bak gmon.out core

defrost: defrost.f90
	$(FC) $(FFLAGS) $(FFTWINC) $(SILOINC) $^ -o $@ $(LDFLAGS) $(FFTWLIB) $(SILOLIB)
