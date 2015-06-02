# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_builddir)/lib

#
# Test code
#
bin_PROGRAMS = Grid_main \
	Grid_simd \
	Grid_gamma  \
	Grid_cshift \
	Grid_cshift_red_black \
	Grid_stencil \
	Grid_nersc_io \
	Grid_rng \
	Grid_remez \
	Grid_rng_fixed \
	Grid_wilson_evenodd \
	Grid_wilson_cg_unprec \
	Grid_wilson_cg_prec \
	Grid_wilson_cg_schur \
	Grid_dwf_even_odd\
	Grid_dwf_cg_unprec\
	Grid_dwf_cg_prec\
	Grid_dwf_cg_schur

test: 
	for f in $bin_PROGRAMS
	do
	./$f > $f.log
	done

Grid_main_SOURCES = Grid_main.cc
Grid_main_LDADD = -lGrid

Grid_rng_SOURCES = Grid_rng.cc
Grid_rng_LDADD = -lGrid

Grid_rng_fixed_SOURCES = Grid_rng_fixed.cc
Grid_rng_fixed_LDADD = -lGrid

Grid_remez_SOURCES = Grid_remez.cc
Grid_remez_LDADD = -lGrid

Grid_nersc_io_SOURCES = Grid_nersc_io.cc
Grid_nersc_io_LDADD = -lGrid

Grid_cshift_SOURCES = Grid_cshift.cc
Grid_cshift_LDADD = -lGrid

Grid_cshift_red_black_SOURCES = Grid_cshift_red_black.cc
Grid_cshift_red_black_LDADD = -lGrid

Grid_gamma_SOURCES = Grid_gamma.cc
Grid_gamma_LDADD = -lGrid

Grid_stencil_SOURCES = Grid_stencil.cc
Grid_stencil_LDADD = -lGrid

Grid_simd_SOURCES = Grid_simd.cc
Grid_simd_LDADD = -lGrid

#Grid_simd_new_SOURCES = Grid_simd_new.cc
#Grid_simd_new_LDADD = -lGrid

Grid_wilson_evenodd_SOURCES = Grid_wilson_evenodd.cc
Grid_wilson_evenodd_LDADD = -lGrid

Grid_wilson_cg_unprec_SOURCES = Grid_wilson_cg_unprec.cc
Grid_wilson_cg_unprec_LDADD = -lGrid

Grid_wilson_cg_prec_SOURCES = Grid_wilson_cg_prec.cc
Grid_wilson_cg_prec_LDADD = -lGrid

Grid_wilson_cg_schur_SOURCES = Grid_wilson_cg_schur.cc
Grid_wilson_cg_schur_LDADD = -lGrid

Grid_dwf_even_odd_SOURCES = Grid_dwf_even_odd.cc
Grid_dwf_even_odd_LDADD = -lGrid

Grid_dwf_cg_unprec_SOURCES = Grid_dwf_cg_unprec.cc
Grid_dwf_cg_unprec_LDADD = -lGrid

Grid_dwf_cg_prec_SOURCES = Grid_dwf_cg_prec.cc
Grid_dwf_cg_prec_LDADD = -lGrid

Grid_dwf_cg_schur_SOURCES = Grid_dwf_cg_schur.cc
Grid_dwf_cg_schur_LDADD = -lGrid

