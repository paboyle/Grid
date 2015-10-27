AM_CXXFLAGS = -I$(top_srcdir)/programs -I$(top_srcdir)/lib
AM_LDFLAGS = -L$(top_builddir)/lib

bin_PROGRAMS = Hadrons

Hadrons_SOURCES =   \
    Application.cc  \
    Application.hpp \
    Global.hpp      \
    Global.cc       \
    Hadrons.cc

Hadrons_LDADD   = -lGrid
