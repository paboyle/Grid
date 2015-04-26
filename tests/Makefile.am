# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_srcdir)/lib

#
# Test code
#
bin_PROGRAMS = Grid_main Grid_stencil Grid_nersc_io Grid_cshift Grid_gamma Grid_wilson Grid_simd

Grid_main_SOURCES = Grid_main.cc
Grid_main_LDADD = -lGrid

Grid_nersc_io_SOURCES = Grid_nersc_io.cc
Grid_nersc_io_LDADD = -lGrid

Grid_cshift_SOURCES = Grid_cshift.cc
Grid_cshift_LDADD = -lGrid

Grid_gamma_SOURCES = Grid_gamma.cc
Grid_gamma_LDADD = -lGrid

Grid_stencil_SOURCES = Grid_stencil.cc
Grid_stencil_LDADD = -lGrid

Grid_wilson_SOURCES = Grid_wilson.cc
Grid_wilson_LDADD = -lGrid

Grid_simd_SOURCES = Grid_simd.cc
Grid_simd_LDADD = -lGrid
