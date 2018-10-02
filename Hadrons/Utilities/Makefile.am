AM_LDFLAGS += -L../../Hadrons

bin_PROGRAMS = HadronsXmlRun FermionEigenPackCastToSingle

HadronsXmlRun_SOURCES = HadronsXmlRun.cc
HadronsXmlRun_LDADD   = -lHadrons -lGrid

FermionEigenPackCastToSingle_SOURCES  = EigenPackCast.cc
FermionEigenPackCastToSingle_CXXFLAGS = $(AM_CXXFLAGS) -DFIN=WilsonImplD::FermionField -DFOUT=WilsonImplF::FermionField
FermionEigenPackCastToSingle_LDADD    = -lHadrons -lGrid
