lib_LIBRARIES = libHadrons.a
bin_PROGRAMS  = HadronsXmlRun HadronsXmlSchedule

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
	Application.hpp           \
	Environment.hpp           \
	Exceptions.hpp            \
	Factory.hpp               \
	GeneticScheduler.hpp      \
	Global.hpp                \
	Graph.hpp                 \
	LanczosUtils.hpp          \
	Module.hpp                \
	Modules.hpp               \
	ModuleFactory.hpp         \
	VirtualMachine.hpp

HadronsXmlRun_SOURCES = HadronsXmlRun.cc
HadronsXmlRun_LDADD   = libHadrons.a -lGrid

HadronsXmlSchedule_SOURCES = HadronsXmlSchedule.cc
HadronsXmlSchedule_LDADD   = libHadrons.a -lGrid
