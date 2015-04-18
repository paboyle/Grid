# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_srcdir)/lib

#
# Test code
#
bin_PROGRAMS = Grid_main test_Grid_stencil

Grid_main_SOURCES = Grid_main.cc
Grid_main_LDADD = -lGrid

test_Grid_stencil_SOURCES = test_Grid_stencil.cc
test_Grid_stencil_LDADD = -lGrid

