# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_builddir)/lib

#
# Test code
#
bin_PROGRAMS = Test_main \
	Test_simd \
	Test_gamma  \
	Test_cshift \
	Test_cshift_red_black \
	Test_stencil \
	Test_nersc_io \
	Test_rng \
	Test_remez \
	Test_rng_fixed \
	Test_wilson_evenodd \
	Test_wilson_cg_unprec \
	Test_wilson_cg_prec \
	Test_wilson_cg_schur \
	Test_dwf_even_odd\
	Test_dwf_cg_unprec\
	Test_dwf_cg_prec\
	Test_dwf_cg_schur\
	Test_many_evenodd\
	Test_many_cg



Test_main_SOURCES = Test_main.cc
Test_main_LDADD = -lGrid

Test_rng_SOURCES = Test_rng.cc
Test_rng_LDADD = -lGrid

Test_rng_fixed_SOURCES = Test_rng_fixed.cc
Test_rng_fixed_LDADD = -lGrid

Test_remez_SOURCES = Test_remez.cc
Test_remez_LDADD = -lGrid

Test_nersc_io_SOURCES = Test_nersc_io.cc
Test_nersc_io_LDADD = -lGrid

Test_cshift_SOURCES = Test_cshift.cc
Test_cshift_LDADD = -lGrid

Test_cshift_red_black_SOURCES = Test_cshift_red_black.cc
Test_cshift_red_black_LDADD = -lGrid

Test_gamma_SOURCES = Test_gamma.cc
Test_gamma_LDADD = -lGrid

Test_stencil_SOURCES = Test_stencil.cc
Test_stencil_LDADD = -lGrid

Test_simd_SOURCES = Test_simd.cc
Test_simd_LDADD = -lGrid

#Test_simd_new_SOURCES = Test_simd_new.cc
#Test_simd_new_LDADD = -lGrid

Test_wilson_evenodd_SOURCES = Test_wilson_evenodd.cc
Test_wilson_evenodd_LDADD = -lGrid

Test_wilson_cg_unprec_SOURCES = Test_wilson_cg_unprec.cc
Test_wilson_cg_unprec_LDADD = -lGrid

Test_wilson_cg_prec_SOURCES = Test_wilson_cg_prec.cc
Test_wilson_cg_prec_LDADD = -lGrid

Test_wilson_cg_schur_SOURCES = Test_wilson_cg_schur.cc
Test_wilson_cg_schur_LDADD = -lGrid

Test_dwf_even_odd_SOURCES = Test_dwf_even_odd.cc
Test_dwf_even_odd_LDADD = -lGrid

Test_dwf_cg_unprec_SOURCES = Test_dwf_cg_unprec.cc
Test_dwf_cg_unprec_LDADD = -lGrid

Test_dwf_cg_prec_SOURCES = Test_dwf_cg_prec.cc
Test_dwf_cg_prec_LDADD = -lGrid

Test_dwf_cg_schur_SOURCES = Test_dwf_cg_schur.cc
Test_dwf_cg_schur_LDADD = -lGrid

Test_many_evenodd_SOURCES = Test_many_evenodd.cc
Test_many_evenodd_LDADD = -lGrid

Test_many_cg_SOURCES = Test_many_cg.cc
Test_many_cg_LDADD = -lGrid
