SUBDIRS = core forces hmc solver debug smearing

if BUILD_CHROMA_REGRESSION
  SUBDIRS+= qdpxx
endif

.PHONY: subtests

include Make.inc

subtests:
	for d in $(SUBDIRS); do $(MAKE) -C $${d} tests; done
