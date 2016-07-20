# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/include/

extra_sources=
if BUILD_COMMS_MPI
  extra_sources+=communicator/Communicator_mpi.cc
endif

if BUILD_COMMS_SHMEM
  extra_sources+=communicator/Communicator_shmem.cc
endif

if BUILD_COMMS_NONE
  extra_sources+=communicator/Communicator_none.cc
endif

#
# Libraries
#
include Make.inc
include Eigen.inc

lib_LIBRARIES = libGrid.a

libGrid_a_SOURCES = $(CCFILES) $(extra_sources)

fftwdir = $(prefix)/lib/
fftw_DATA = libfftw3.a

#
# Include files
#
otherincludedir = $(includedir)/Grid
nobase_otherinclude_HEADERS =$(HFILES) $(EFILES) fftw3.h
