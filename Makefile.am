# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/
SUBDIRS = lib tests benchmarks docs


if DOXYGEN_DOC
directory = $(top_srcdir)/docs/

doxyfile:
        (cd $(directory) && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1; 
endif 
