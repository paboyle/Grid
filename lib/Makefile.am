# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/


extra_sources=
if BUILD_COMMS_MPI
  extra_sources+=Grid_communicator_mpi.cc
  extra_sources+=Grid_stencil_common.cc
endif
if BUILD_COMMS_NONE
  extra_sources+=Grid_communicator_fake.cc
  extra_sources+=Grid_stencil_common.cc
endif

#
# Libraries
#
lib_LIBRARIES = libGrid.a
libGrid_a_SOURCES = Grid_init.cc $(extra_sources)

#
# Include files
#
include_HEADERS = Grid_config.h\
	Grid.h\
	Grid_QCD.h\
	Grid_aligned_allocator.h\
	Grid_cartesian.h\
	Grid_cartesian_base.h\
	Grid_cartesian_full.h\
	Grid_cartesian_red_black.h\
	Grid_communicator.h\
	Grid_comparison.h\
	Grid_config.h\
	Grid_cshift.h\
	Grid_cshift_common.h\
	Grid_cshift_mpi.h\
	Grid_cshift_none.h\
	Grid_lattice.h\
	Grid_math.h\
	Grid_math_arith.h\
	Grid_math_tensors.h\
	Grid_math_traits.h\
	Grid_predicated.h\
	Grid_simd.h\
	Grid_stencil.h\
	Grid_summation.h\
	Grid_vComplexD.h\
	Grid_vComplexF.h\
	Grid_vInteger.h\
	Grid_vRealD.h\
	Grid_vRealF.h



