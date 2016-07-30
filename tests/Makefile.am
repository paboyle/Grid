SUBDIRS = core forces hmc solver debug	

if BUILD_CHROMA_REGRESSION
  SUBDIRS+= qdpxx
endif

AM_LDFLAGS = -L$(top_builddir)/lib

if USE_LAPACK
AM_CXXFLAGS += -DUSE_LAPACK
if USE_LAPACK_LIB
#if test "X${ac_LAPACK}X" != XyesX 
AM_CXXFLAGS += -I$(ac_LAPACK)/include
AM_LDFLAGS += -L$(ac_LAPACK)/lib
#fi
endif
endif

include Make.inc
