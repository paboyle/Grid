# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_builddir)/lib

#
# Test code
#
bin_PROGRAMS = Grid_main Grid_stencil Grid_nersc_io Grid_cshift Grid_gamma  Grid_simd Grid_rng Grid_remez Grid_rng_fixed

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

Grid_gamma_SOURCES = Grid_gamma.cc
Grid_gamma_LDADD = -lGrid

Grid_stencil_SOURCES = Grid_stencil.cc
Grid_stencil_LDADD = -lGrid

Grid_simd_SOURCES = Grid_simd.cc
Grid_simd_LDADD = -lGrid
