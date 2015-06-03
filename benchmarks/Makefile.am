# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_builddir)/lib

#
# Test code
#
bin_PROGRAMS = \
	Benchmark_comms \
	Benchmark_memory_bandwidth \
	Benchmark_su3 \
	Benchmark_wilson \
	Benchmark_dwf

Benchmark_comms_SOURCES = Benchmark_comms.cc
Benchmark_comms_LDADD = -lGrid

Benchmark_memory_bandwidth_SOURCES = Benchmark_memory_bandwidth.cc
Benchmark_memory_bandwidth_LDADD = -lGrid

Benchmark_su3_SOURCES = Benchmark_su3.cc Benchmark_su3_test.cc Benchmark_su3_expr.cc
Benchmark_su3_LDADD = -lGrid

Benchmark_wilson_SOURCES = Benchmark_wilson.cc
Benchmark_wilson_LDADD = -lGrid

Benchmark_dwf_SOURCES = Benchmark_dwf.cc
Benchmark_dwf_LDADD = -lGrid

