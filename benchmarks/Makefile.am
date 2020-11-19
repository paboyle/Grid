include Make.inc

bench-local: all
	./Benchmark_su3
	./Benchmark_memory_bandwidth
	./Benchmark_wilson
	./Benchmark_dwf --dslash-unroll

OBJECTS=$(shell find ../Grid/ -name '*.o') 

Benchmark_dwf_fp32.nolib:	Benchmark_dwf_fp32.o
	$(CXXLINK)  Benchmark_dwf_fp32.o $(OBJECTS) $(LIBS)

