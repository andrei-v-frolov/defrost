D E F R O S T   ---   R E A D M E
=================================

##### Reheating code doing something... #####
<http://www.sfu.ca/physics/cosmology/defrost>

Copyright (C) 2007-2012 Andrei Frolov <<frolov@sfu.ca>>  
Distributed under the terms of the GNU General Public License  
If you use this code for your research, please cite
[arXiv:0809.4904] (http://arxiv.org/abs/0809.4904)

***

Included in this distribution are:
	
	README		- this file
	Makefile	- build script
	defrost.f90	- DEFROST source code

***

QUICK START:
============

Download and install required libraries (FFTW, optionally SILO/HDF5).
Edit `Makefile` for your Fortran compiler settings and library locations.
Edit `defrost.f90` specifying your reheating model (number of scalar
fields, initial conditions, potential and its derivatives), configure
the solver (box size, spatial and time resolutions, etc) and output
settings. Build `defrost` binary by saying `make` and run it, optionally
saving the log `defrost > LOG`.

**CAUTION:** Being a 3D PDE solver, DEFROST needs a lot of resources.
Active memory footprint while running 256^3 simulation is around 512Mb,
and the execution speed on a dual Xeon 5160 machine (64 bit, 4 cores at
3.0GHz, quad-channel DDR2 RAM) is 0.32s per time step (for 2 field
model, second order integrator). If you enable 3D output, a lot of data
will get dumped to your disk. The size of a single 256^3 frame
containing all variables is about 1Gb.

USING DEFROST:
==============

#### 0) PREREQUISITES: ####

DEFROST is written in Fortran-90 (with Fortran preprocessor used).
The code is developed and optimized for Intel Fortran compiler (`ifort`).
GNU Fortran compilers (`g77`, `g95`, and `gfortran`) are *not* supported.

DEFROST uses [FFTW library] (http://www.fftw.org/) for Fourier transforms.
Download and compile FFTW version 3.0 or later (tested with FFTW 3.1.2),
including Fortran wrappers and threaded FFT support. Fedora users can
install pre-compiled binaries by saying `yum install fftw fftw-devel`.

DEFROST output is instrumented for interactive 3D visualization with
[LLNL VisIt] (http://www.llnl.gov/visit/). DEFROST can dump 3D data in
native VisIt format (using SILO library) or raw box-of-values format
(simple to read and write, but missing some metadata and requiring many
files per frame as opposed to a single SILO database).

If you would like to use VisIt native SILO format (recommended), you
will need to install SILO and HDF5 libraries. Download and install
[HDF5 library] (http://www.hdfgroup.org/) first, which SILO uses as a
back-end. Fedora users can install pre-compiled binaries by saying
`yum install hdf5 hdf5-devel`. Download and install [SILO library]
(https://wci.llnl.gov/codes/silo/) next. Getting older versions of
SILO to build against shared HDF5 library required [patching]
(http://www.sfu.ca/physics/cosmology/defrost/silo/README.html),
but this is not a problem anymore.

On a Mac, all the required libraries are available from [MacPorts]
(http://www.macports.org/).

#### 1) CONFIGURE MAKEFILE: ####

The provided `Makefile` is configured for optimized parallel build using
Intel Fortran compiler, threaded FFTs and SILO output. Unfortunately,
GNU Fortran compiler produces (very) much slower code and cannot handle
parallelism, so the support for it was dropped. Compiling with other
Fortran compilers should be possible, but was not tested.

Adjust `FFLAGS` to optimize for your machine (in particular, `-x` flag
selects target CPU). As DEFROST uses Fortran preprocessor, `-fpp` flag
is necessary for it to compile. Floating point numbers are set to double
precision by default (`-r8` flag), but the code can be compiled for
single precision (`-r4`) to decrease memory footprint, or quad precision
(`-r16`) to improve roundoff errors (at 10x runtime penalty). Regardless
of solver precision, the bindings to FFTW and SILO will use double
precision floats.

DEFROST solver code is fully parallel, and can be build as such using
either automatic parallelizer (`-parallel`) or OpenMP (`-openmp`).
Automatic parallelizer will produce slightly more efficient code.
Threaded FFTs will be used if `THREADS` variable is set; comment it out
if you don't want that (or you don't have libfftw3_threads).

By default, DEFROST uses automatic array variables, which get allocated
on stack. This can easily exceed default process limits set by the shell
(especially if OpenMP is used). If you get immediate segmentation fault
on startup without any error messages, try overriding stack size limit
(`ulimit -s unlimited` in bash), or use dynamically allocated arrays (at
small performance penalty) by defining `DYNAMIC_ARRAYS` in Makefile.
Changing compiler memory model by using `-mcmodel=large` flag might be
better instead.

Finally, give locations to search for header files and libraries (using
`-I...` and `-L...` flags correspondingly), if your FFTW and SILO/HDF5
libraries are not in the standard locations. If you do not have (or want
to use) SILO libraries, comment out `SILO` variable, and DEFROST will
compile for raw output instead.

#### 2) MODIFY SOURCE CODE: ####

There is a number of parameters defined in DEFROST source code which
control the solver and the output behavior, and are intended to be
modified by user. They are documented inline. I tried to separate
reheating model from numerical implementation as much as practical, but
if you want to change the model, you will have to edit the source code
directly. Version 1.x code is described in [arXiv:0809.4904]
(http://arxiv.org/abs/0809.4904). Version 2.x uses new (and awesome)
symplectic PDE integrators we developed, and will be documented in
upcoming publication.

The 3D field evolution equations are discretized on a `n^3` uniform grid
with spacing `dx`. Time step `dt` has to be small enough to satisfy
Courant's conditions, and to resolve oscillations of all the fields.
Parameter `m` (padded grid size) deserves special explanation. Inner
evolution loop of DEFROST is highly optimized, and runs almost at the
memory bandwidth speed. If array access pattern stride is commensurate
with cache page size, memory access performance can be affected (this is
known as "cache thrashing"). In the worst case, the effect can be as
much as 2x slowdown. Padding array by a few empty elements avoids it.
Try setting `m+1` prime, or a product of a few large primes.

DEFROST has flexible output options (controlled by `output$...`
parameters) and lazy calculation policy (only stuff you want gets
calculated). Disabling unneeded output can save execution time and disk
space. Select quantities {`bov`/`psd`/`cdf`} and variables
{`fld`/`vel`/`rho`/`prs`/`pot`} you want output, and select curve data
format {`gnu`/`vis`}. The output is dumped every `nt` time steps, and 3D
image can be downsampled (using nearest neighbour) to save disk space.
By default, dilution due to expansion is scaled out. All output (except
execution log to stdout) can be disabled by setting `output` parameter
to `.false.`. For fast logging of expansion history *only*, calls to
`checkpt()` inside the main evolution loop can be omitted.

For performance reasons, scalar field potential and its derivatives had
to be inlined separately inside the evolution loop in `...step()`
routines. In addition, human-readable potential string is output in
`head()` routine. These are defined in the model section using
preprocessor macros. Don't forget that potential enters indirectly into
initial conditions (both for homogeneous field components and
fluctuations), so they will need to be adjusted as well.

#### 3) BUILD AND RUN: ####

This is easy. Say `make` and `defrost`.

DEFROST will print expansion history and average density and pressure on
`stdout` as it goes along. Say `defrost > LOG` if you want to keep it.

PLOTTING THE RESULTS:
=====================

#### 1) FIELD STATISTICS: ####

DEFROST can output spectra and distributions of select variables
(enabled by default). The file `CDF` contains cumulative distribution
function, and file `PSD` contains power spectral density (logarithmic in
base 10, in arbitrary units). Both files are in plain text format, with
brief information header (lines starting with comment symbol '#'),
followed by data. Values for all the fields appear on the same line in
order specified in the header, separated by white space. Blocks of
values at different time (frames) are separated by two empty lines.

This data format is directly used by [gnuplot] (http://www.gnuplot.info/),
and should be easy to convert for your favourite plotting program.
Some gnuplot examples:

    plot 'LOG' using 1:2 with lines                 # expansion history
    plot 'LOG' using 1:($5/$4) with lines           # equation of state
    plot 'LOG' using 1:(($4/3-$3**2)*$2**2) w l     # flatness violation

    plot 'PSD' i 0:512:16 u 2:((10**$3)*$2**4) w l  # inflaton spectra

#### 2) 3D VISUALIZATIONS: ####

DEFROST can output 3D snapshots of select variables periodically during
the evolution (disabled by default due to high storage requirements, set
`output$bov` to enable). The output is instrumented for interactive 3D
visualization with [LLNL VisIt] (http://www.llnl.gov/visit/). Check
examples on DEFROST website or the image gallery on VisIt website to see
what you can do with it. If you compiled DEFROST with SILO support, you
can just open `frame-*.silo` database to access DEFROST simulation data.
Otherwise, you will have to open databases for each variable separately
(they are named `xxx-*.bov`). Set `output$vis` to export statistics to
VisIt as well.

#### 3) MAKING MOVIES: ####

VisIt can make a movie of your 3D visualization. Set up your plots, tidy
up using `Controls > Annotation...`, and go to `File > Save Movie...`
wizard. Movie encoder shipped with VisIt distribution is not very good,
so I recommend saving movie frames as images (in a lossless format such
as PNG), and encoding using stand-alone encoder. I highly recommend
[mplayer/mencoder] (http://www.mplayerhq.hu/) with lavc codec (builtin)
or [x264 codec] (http://www.videolan.org/developers/x264.html). If you
have deep pockets and want GUI and lots of presets, [Adobe Media Encoder]
(http://www.adobe.com/devnet/flash/quickstart/video_encoder.html) (ships
with CS5.5) might be an option.

Mencoder is a very capable program, but can be intimidating to use
because of command line interface and multitude of options. Grab
[my presets] (http://www.sfu.ca/physics/cosmology/defrost/mencoder.conf)
to see encoding settings I use. The encoding itself is best done in two
passes:

    mencoder -profile x264 -x264encopts pass=1 mf://"/tmp/movie_*.png" -o test.avi
    mencoder -profile x264 -x264encopts pass=2 mf://"/tmp/movie_*.png" -o test.avi

After encoding is done, test the playback with `mplayer test.avi`. If
you want to play the movie from other applications (e.g. to include it
in a PowerPoint presentation), you will need to install system-wide
MPEG4/H264 codecs. I recommend [Combined Community Codec Pack]
(http://cccp-project.net/) on Windows, and [Perian] (http://perian.org/)
for QuickTime support on Mac.

FINAL COMMENTS:
===============

Bug reports, patches, suggestions, flames, etc. should go to the author.
Have fun! ^_^

***

Andrei Frolov <<frolov@sfu.ca>>  -- $Id: README,v 1.4 2012/04/20 22:51:33 frolov Stab $
