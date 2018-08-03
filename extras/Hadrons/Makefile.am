lib_LIBRARIES = libHadrons.a
bin_PROGRAMS  = HadronsXmlRun

include modules.inc

libHadrons_a_SOURCES = \
    $(modules_cc)      \
    Application.cc     \
    Environment.cc     \
	Exceptions.cc      \
    Global.cc          \
    Module.cc		   \
	VirtualMachine.cc
libHadrons_adir = $(pkgincludedir)/Hadrons
nobase_libHadrons_a_HEADERS = \
	$(modules_hpp)            \
	AllToAllVectors.hpp       \
	Application.hpp           \
	EigenPack.hpp             \
	Environment.hpp           \
	Exceptions.hpp            \
	Factory.hpp               \
	GeneticScheduler.hpp      \
	Global.hpp                \
	Graph.hpp                 \
	Module.hpp                \
	Modules.hpp               \
	ModuleFactory.hpp         \
	Solver.hpp                \
	VirtualMachine.hpp

HadronsXmlRun_SOURCES = HadronsXmlRun.cc
HadronsXmlRun_LDADD   = libHadrons.a -lGrid
