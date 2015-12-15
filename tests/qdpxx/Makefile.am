# additional include paths necessary to compile the C++ library

AM_CXXFLAGS = -I$(top_srcdir)/lib `chroma-config --cxxflags`
AM_LDFLAGS = -L$(top_builddir)/lib `chroma-config --ldflags` `chroma-config --libs`


include Make.inc
