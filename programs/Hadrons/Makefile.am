AM_CXXFLAGS = -I$(top_srcdir)/programs -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_builddir)/lib

bin_PROGRAMS = Hadrons

Hadrons_SOURCES =   \
    Application.cc  \
    Environment.cc  \
    Global.cc       \
    Hadrons.cc      \
    Module.cc       \
    ModuleFactory.cc

Hadrons_SOURCES +=  \
    MQuark.cc

Hadrons_SOURCES +=  \
    CMeson.cc

Hadrons_LDADD   = -lGrid
