AM_CXXFLAGS = -I$(top_srcdir)/programs -I$(top_srcdir)/lib
AM_LDFLAGS  = -L$(top_builddir)/lib

bin_PROGRAMS = Hadrons

## general sources
Hadrons_SOURCES =   \
    Application.cc  \
    Environment.cc  \
    FermionAction.cc\
    Global.cc       \
    Hadrons.cc      \
    Module.cc

## fermion actions
Hadrons_SOURCES +=  \
    AWilson.cc

## general modules
Hadrons_SOURCES +=  \
    MQuark.cc       \
    MSource.cc

## solver modules
Hadrons_SOURCES +=  \
    SRBPrecCG.cc

## contraction modules
Hadrons_SOURCES +=  \
    CMeson.cc

Hadrons_LDADD   = -lGrid
