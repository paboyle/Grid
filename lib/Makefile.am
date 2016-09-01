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

lib_LTLIBRARIES = libGrid.la

libGrid_la_SOURCES             = $(CCFILES) $(extra_sources)
libGrid_ladir                  = $(pkgincludedir)
nobase_dist_pkginclude_HEADERS = $(HFILES) $(eigen_files) Config.h
