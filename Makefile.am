# additional include paths necessary to compile the C++ library
SUBDIRS = Grid Hadrons benchmarks tests

include $(top_srcdir)/doxygen.inc

bin_SCRIPTS=grid-config

BUILT_SOURCES = version.h

version.h:
	if [ `git status --porcelain | wc -l` -gt 0 ]; then a="uncommited changes"; else a="no uncommitted changes"; fi; echo "`git log -n 1 --format=format:"#define GITHASH \\"%H:%d $$a\\"%n" HEAD`" > $(srcdir)/Grid/version.h

.PHONY: bench check tests doxygen-run doxygen-doc $(DX_PS_GOAL) $(DX_PDF_GOAL)

tests-local: all
bench-local: all
check-local: all

AM_CXXFLAGS += -I$(top_builddir)/include

ACLOCAL_AMFLAGS = -I m4
