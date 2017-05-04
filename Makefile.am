# additional include paths necessary to compile the C++ library
SUBDIRS = lib benchmarks tests extras

include $(top_srcdir)/doxygen.inc

bin_SCRIPTS=grid-config

tests: all
	$(MAKE) -C tests tests

.PHONY: tests doxygen-run doxygen-doc $(DX_PS_GOAL) $(DX_PDF_GOAL)

AM_CXXFLAGS += -I$(top_builddir)/include
ACLOCAL_AMFLAGS = -I m4
