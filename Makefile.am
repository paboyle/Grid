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
include_HEADERS = Grid.h\
	Grid_vComplexD.h\
	Grid_vComplexF.h\
	Grid_vRealD.h\
	Grid_vRealF.h\
	Grid_Cartesian.h\
	Grid_Lattice.h\
	Grid_config.h

#
# Test code
#
bin_PROGRAMS = Grid_main
Grid_main_SOURCES = \
	Grid_main.cc\
	Grid_fake.cc

Grid_main_LDADD = libGrid.a
