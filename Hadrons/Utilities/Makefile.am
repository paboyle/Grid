bin_PROGRAMS = HadronsXmlRun HadronsFermionEP64To32

HadronsXmlRun_SOURCES = HadronsXmlRun.cc
HadronsXmlRun_LDADD   = ../libHadrons.a ../../Grid/libGrid.a

HadronsFermionEP64To32_SOURCES  = EigenPackCast.cc
HadronsFermionEP64To32_CXXFLAGS = $(AM_CXXFLAGS) -DFIN=WilsonImplD::FermionField -DFOUT=WilsonImplF::FermionField
HadronsFermionEP64To32_LDADD    = ../libHadrons.a ../../Grid/libGrid.a
