# additional include paths necessary to compile the C++ library
SUBDIRS = lib benchmarks tests extras

include $(top_srcdir)/doxygen.inc

bin_SCRIPTS=grid-config

BUILT_SOURCES = version.h

version.h:
	echo "`git log -n 1 --format=format:"#define GITHASH \\"%H:%d\\"%n" HEAD`" > $(srcdir)/lib/version.h

.PHONY: bench check tests doxygen-run doxygen-doc $(DX_PS_GOAL) $(DX_PDF_GOAL)

tests-local: all
bench-local: all
check-local: all

AM_CXXFLAGS += -I$(top_builddir)/include

ACLOCAL_AMFLAGS = -I m4
