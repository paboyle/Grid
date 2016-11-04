SUBDIRS = core forces hmc solver debug	

if BUILD_CHROMA_REGRESSION
  SUBDIRS+= qdpxx
endif

.PHONY: subtests

include Make.inc

subtests:
	for d in $(SUBDIRS); do make -C $${d} tests; done
