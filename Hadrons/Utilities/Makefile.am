AM_LDFLAGS += -L../../Hadrons

bin_PROGRAMS = HadronsXmlRun HadronsFermionEP64To32

HadronsXmlRun_SOURCES = HadronsXmlRun.cc
HadronsXmlRun_LDADD   = -lHadrons -lGrid

HadronsFermionEP64To32_SOURCES  = EigenPackCast.cc
HadronsFermionEP64To32_CXXFLAGS = $(AM_CXXFLAGS) -DFIN=WilsonImplD::FermionField -DFOUT=WilsonImplF::FermionField
HadronsFermionEP64To32_LDADD    = -lHadrons -lGrid
