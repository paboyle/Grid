# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_builddir)/lib

#
# Test code
#
bin_PROGRAMS = Grid_wilson Grid_comms Grid_memory_bandwidth Grid_su3 Grid_wilson_cg_unprec

Grid_wilson_SOURCES = Grid_wilson.cc
Grid_wilson_LDADD = -lGrid

Grid_wilson_cg_unprec_SOURCES = Grid_wilson_cg_unprec.cc
Grid_wilson_cg_unprec_LDADD = -lGrid

Grid_comms_SOURCES = Grid_comms.cc
Grid_comms_LDADD = -lGrid

Grid_su3_SOURCES = Grid_su3.cc
Grid_su3_LDADD = -lGrid

Grid_memory_bandwidth_SOURCES = Grid_memory_bandwidth.cc
Grid_memory_bandwidth_LDADD = -lGrid

