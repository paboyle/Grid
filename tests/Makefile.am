# additional include paths necessary to compile the C++ library

SUBDIRS = 
if BUILD_CHROMA_REGRESSION
  SUBDIRS+= qdpxx
endif

AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_builddir)/lib

if BUILD_ZMM
  bin_PROGRAMS=Test_zmm
endif

include Make.inc
