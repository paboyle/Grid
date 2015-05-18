# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/

extra_sources=
if BUILD_COMMS_MPI
  extra_sources+=communicator/Grid_communicator_mpi.cc
endif

if BUILD_COMMS_NONE
  extra_sources+=communicator/Grid_communicator_none.cc
endif

#
# Libraries
#
lib_LIBRARIES = libGrid.a
libGrid_a_SOURCES =\
	 Grid_init.cc\
	stencil/Grid_stencil_common.cc\
	qcd/Grid_qcd_dirac.cc\
	qcd/Grid_qcd_wilson_dop.cc\
	algorithms/approx/Zolotarev.cc\
	algorithms/approx/Remez.cc\
	$(extra_sources)

#
# Include files
#
nobase_include_HEADERS = algorithms/approx/bigfloat.h\
	algorithms/approx/Chebyshev.h\
	algorithms/approx/Remez.h\
	algorithms/approx/Zolotarev.h\
	algorithms/iterative/ConjugateGradient.h\
	algorithms/iterative/NormalEquations.h\
	algorithms/iterative/SchurRedBlack.h\
	algorithms/LinearOperator.h\
	algorithms/SparseMatrix.h\
	cartesian/Grid_cartesian_base.h\
	cartesian/Grid_cartesian_full.h\
	cartesian/Grid_cartesian_red_black.h\
	communicator/Grid_communicator_base.h\
	cshift/Grid_cshift_common.h\
	cshift/Grid_cshift_mpi.h\
	cshift/Grid_cshift_none.h\
	Grid.h\
	Grid_algorithms.h\
	Grid_aligned_allocator.h\
	Grid_cartesian.h\
	Grid_communicator.h\
	Grid_comparison.h\
	Grid_cshift.h\
	Grid_extract.h\
	Grid_lattice.h\
	Grid_math.h\
	Grid_simd.h\
	Grid_stencil.h\
	Grid_threads.h\
	lattice/Grid_lattice_arith.h\
	lattice/Grid_lattice_base.h\
	lattice/Grid_lattice_comparison.h\
	lattice/Grid_lattice_conformable.h\
	lattice/Grid_lattice_coordinate.h\
	lattice/Grid_lattice_ET.h\
	lattice/Grid_lattice_local.h\
	lattice/Grid_lattice_overload.h\
	lattice/Grid_lattice_peekpoke.h\
	lattice/Grid_lattice_reality.h\
	lattice/Grid_lattice_reduction.h\
	lattice/Grid_lattice_rng.h\
	lattice/Grid_lattice_trace.h\
	lattice/Grid_lattice_transfer.h\
	lattice/Grid_lattice_transpose.h\
	lattice/Grid_lattice_where.h\
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
	parallelIO/GridNerscIO.h\
	qcd/Grid_qcd.h\
	qcd/Grid_qcd_2spinor.h\
	qcd/Grid_qcd_dirac.h\
	qcd/Grid_qcd_wilson_dop.h\
	simd/Grid_vComplexD.h\
	simd/Grid_vComplexF.h\
	simd/Grid_vInteger.h\
	simd/Grid_vRealD.h\
	simd/Grid_vRealF.h

