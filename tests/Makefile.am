# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_srcdir)/lib

#
# Test code
#
bin_PROGRAMS = Grid_main test_Grid_stencil Grid_nersc_io Grid_cshift

Grid_main_SOURCES = Grid_main.cc
Grid_main_LDADD = -lGrid

Grid_nersc_io_SOURCES = Grid_nersc_io.cc
Grid_nersc_io_LDADD = -lGrid

Grid_cshift_SOURCES = Grid_cshift.cc
Grid_cshift_LDADD = -lGrid

test_Grid_stencil_SOURCES = test_Grid_stencil.cc
test_Grid_stencil_LDADD = -lGrid
