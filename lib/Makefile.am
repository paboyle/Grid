# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/

extra_sources=
if BUILD_COMMS_MPI
  extra_sources+=communicator/Grid_communicator_mpi.cc
  extra_sources+=stencil/Grid_stencil_common.cc
endif
if BUILD_COMMS_NONE
  extra_sources+=communicator/Grid_communicator_fake.cc
  extra_sources+=stencil/Grid_stencil_common.cc
endif

#
# Libraries
#
lib_LIBRARIES = libGrid.a
libGrid_a_SOURCES = Grid_init.cc $(extra_sources)

#
# Include files
#
include_HEADERS =\
	Grid.h\
	Grid_aligned_allocator.h\
	Grid_cartesian.h\
	Grid_communicator.h\
	Grid_comparison.h\
	Grid_config.h\
	Grid_cshift.h\
	Grid_lattice.h\
	Grid_math.h\
	Grid_simd.h\
	Grid_stencil.h\
	Grid_summation.h\
	Grid_where.h

nobase_include_HEADERS=\
	cartesian/Grid_cartesian_base.h\
	cartesian/Grid_cartesian_full.h\
	cartesian/Grid_cartesian_red_black.h\
	cshift/Grid_cshift_common.h\
	cshift/Grid_cshift_mpi.h\
	cshift/Grid_cshift_none.h\
	lattice/Grid_lattice_arith.h\
	lattice/Grid_lattice_conformable.h\
	lattice/Grid_lattice_coordinate.h\
	lattice/Grid_lattice_local.h\
	lattice/Grid_lattice_peekpoke.h\
	lattice/Grid_lattice_reality.h\
	lattice/Grid_lattice_reduction.h\
	lattice/Grid_lattice_rng.h\
	lattice/Grid_lattice_trace.h\
	lattice/Grid_lattice_transfer.h\
	lattice/Grid_lattice_transpose.h\
	math/Grid_math_arith.h\
	math/Grid_math_arith_add.h\
	math/Grid_math_arith_mac.h\
	math/Grid_math_arith_mul.h\
	math/Grid_math_arith_scalar.h\
	math/Grid_math_arith_sub.h\
	math/Grid_math_inner.h\
	math/Grid_math_outer.h\
	math/Grid_math_peek.h\
	math/Grid_math_poke.h\
	math/Grid_math_reality.h\
	math/Grid_math_tensors.h\
	math/Grid_math_trace.h\
	math/Grid_math_traits.h\
	math/Grid_math_transpose.h\
	qcd/Grid_QCD.h\
	simd/Grid_vComplexD.h\
	simd/Grid_vComplexF.h\
	simd/Grid_vInteger.h\
	simd/Grid_vRealD.h\
	simd/Grid_vRealF.h

