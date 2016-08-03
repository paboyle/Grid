AM_CXXFLAGS += -I$(top_srcdir)/programs -I../$(top_srcdir)/programs

bin_PROGRAMS = Hadrons

include modules.inc

Hadrons_SOURCES =   \
    $(modules)      \
    Application.cc  \
    Environment.cc  \
    Global.cc       \
    Hadrons.cc      \
    Module.cc

Hadrons_LDADD   = -lGrid
