AM_CXXFLAGS = -I$(top_srcdir)/programs -I$(top_srcdir)/lib
AM_LDFLAGS  = -L$(top_builddir)/lib

bin_PROGRAMS = Hadrons

# general sources
Hadrons_SOURCES =   \
    Application.cc  \
    Environment.cc  \
    Global.cc       \
    Hadrons.cc      \
    Module.cc

# general modules
Hadrons_SOURCES +=  \
    MQuark.cc       \
    MSource.cc

# fermion actions
Hadrons_SOURCES +=  \
    AWilson.cc

# contraction modules
Hadrons_SOURCES +=  \
    CMeson.cc

# gauge modules
Hadrons_SOURCES +=  \
    GLoad.cc        \
    GRandom.cc      \
    GUnit.cc

# solver modules
Hadrons_SOURCES +=  \
    SolRBPrecCG.cc

Hadrons_LDADD   = -lGrid
