# additional include paths necessary to compile the C++ library
SUBDIRS = lib benchmarks tests

.PHONY: tests

tests:
	$(MAKE) -C tests tests

AM_CXXFLAGS += -I$(top_builddir)/include
ACLOCAL_AMFLAGS = -I m4
