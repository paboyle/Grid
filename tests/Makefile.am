SUBDIRS = core forces hmc solver debug	

if BUILD_CHROMA_REGRESSION
  SUBDIRS+= qdpxx
endif

AM_LDFLAGS += -L$(top_builddir)/lib

include Make.inc
