extra_sources=
extra_headers=

if BUILD_COMMS_MPI3
  extra_sources+=communicator/Communicator_mpi3.cc
  extra_sources+=communicator/Communicator_base.cc
  extra_sources+=communicator/SharedMemoryMPI.cc
  extra_sources+=communicator/SharedMemory.cc
endif

if BUILD_COMMS_NONE
  extra_sources+=communicator/Communicator_none.cc
  extra_sources+=communicator/Communicator_base.cc
  extra_sources+=communicator/SharedMemoryNone.cc
  extra_sources+=communicator/SharedMemory.cc
endif

if BUILD_HDF5
  extra_sources+=serialisation/Hdf5IO.cc 
  extra_headers+=serialisation/Hdf5IO.h
  extra_headers+=serialisation/Hdf5Type.h
endif

#
# Libraries
#
include Make.inc
include Eigen.inc

lib_LIBRARIES = libGrid.a

CCFILES += $(extra_sources)
HFILES  += $(extra_headers)

libGrid_a_SOURCES              = $(CCFILES)
libGrid_adir                   = $(pkgincludedir)
nobase_dist_pkginclude_HEADERS = $(HFILES) $(eigen_files) Config.h
