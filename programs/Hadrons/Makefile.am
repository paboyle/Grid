AM_CXXFLAGS = -I$(top_srcdir)/programs -I../$(top_srcdir)/programs -I$(top_srcdir)/lib
AM_LDFLAGS  = -L$(top_builddir)/lib

bin_PROGRAMS = Hadrons

# general sources
Hadrons_SOURCES =   \
    Application.cc  \
    Environment.cc  \
    Global.cc       \
    Hadrons.cc      \
    Module.cc

# modules
include modules.inc
Hadrons_SOURCES += $(modules)

Hadrons_LDADD   = -lGrid
