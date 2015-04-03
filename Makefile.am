# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/

#
# Libraries
#
lib_LIBRARIES = libGrid.a
libGrid_a_SOURCES = Grid_init.cc

#
# Include files
#
include_HEADERS = Grid_config.h\
	Grid.h\
	Grid_simd.h\
	Grid_vComplexD.h\
	Grid_vComplexF.h\
	Grid_vRealD.h\
	Grid_vRealF.h\
	Grid_Cartesian.h\
	Grid_Lattice.h\
	Grid_Communicator.h\
	Grid_QCD.h\
	Grid_aligned_allocator.h\
	Grid_cshift.h\
	Grid_cshift_common.h\
	Grid_cshift_fake.h\
	Grid_cshift_mpi.h\
	Grid_cshift_none.h\
	Grid_math_types.h


#
# Test code
#
bin_PROGRAMS = Grid_main

extra_sources=
if BUILD_COMMS_MPI
  extra_sources+=Grid_mpi.cc
endif
if BUILD_COMMS_FAKE
  extra_sources+=Grid_fake.cc
endif
if BUILD_COMMS_NONE
  extra_sources+=Grid_fake.cc
endif

Grid_main_SOURCES = \
	Grid_main.cc\
	$(extra_sources)

Grid_main_LDADD = libGrid.a
