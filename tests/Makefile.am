SUBDIRS = . core forces hmc solver debug smearing 

if BUILD_CHROMA_REGRESSION
  SUBDIRS+= qdpxx
endif

include Make.inc

check-local: all
	./Test_simd
	./Test_cshift
	./Test_stencil
	./Test_cayley_even_odd_vec
	./Test_dwf_mixedcg_prec



