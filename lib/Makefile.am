# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/

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

lib_LIBRARIES = libGrid.a
libGrid_a_SOURCES = $(CCFILES) $(extra_sources)


#	qcd/action/fermion/PartialFractionFermion5D.cc\	\
#
# Include files
#
nobase_include_HEADERS=$(HFILES)

