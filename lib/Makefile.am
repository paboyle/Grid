extra_sources=
if BUILD_COMMS_MPI
  extra_sources+=communicator/Communicator_mpi.cc
  extra_sources+=communicator/Communicator_base.cc
endif

if BUILD_COMMS_MPI3
  extra_sources+=communicator/Communicator_mpi3.cc
  extra_sources+=communicator/Communicator_base.cc
endif

if BUILD_COMMS_SHMEM
  extra_sources+=communicator/Communicator_shmem.cc
  extra_sources+=communicator/Communicator_base.cc
endif

if BUILD_COMMS_NONE
  extra_sources+=communicator/Communicator_none.cc
  extra_sources+=communicator/Communicator_base.cc
endif

#
# Libraries
#
include Make.inc
include Eigen.inc

lib_LIBRARIES = libGrid.a

libGrid_a_SOURCES              = $(CCFILES) $(extra_sources)
libGrid_adir                   = $(pkgincludedir)
nobase_dist_pkginclude_HEADERS = $(HFILES) $(eigen_files) Config.h
