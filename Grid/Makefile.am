extra_sources=
extra_headers=

if BUILD_COMMS_MPI3
  extra_sources+=communicator/Communicator_mpi3.cc
  extra_sources+=communicator/Communicator_base.cc
  extra_sources+=communicator/SharedMemoryMPI.cc
  extra_sources+=communicator/SharedMemory.cc
endif

if BUILD_COMMS_NONE
  extra_sources+=communicator/Communicator_none.cc
  extra_sources+=communicator/Communicator_base.cc
  extra_sources+=communicator/SharedMemoryNone.cc
  extra_sources+=communicator/SharedMemory.cc
endif

if BUILD_HDF5
  extra_sources+=serialisation/Hdf5IO.cc 
  extra_headers+=serialisation/Hdf5IO.h
  extra_headers+=serialisation/Hdf5Type.h
endif


all: version-cache Version.h

version-cache:
	@if [ `git status --porcelain | grep -v '??' | wc -l` -gt 0 ]; then\
		a="uncommited changes";\
	else\
		a="clean";\
	fi;\
	echo "`git log -n 1 --format=format:"#define GITHASH \\"%H:%d $$a\\"%n" HEAD`" > vertmp;\
	if [ -e version-cache ]; then\
		d=`diff vertmp version-cache`;\
		if [ "$${d}" != "" ]; then\
			mv vertmp version-cache;\
			rm -f Version.h;\
		fi;\
	else\
		mv vertmp version-cache;\
		rm -f Version.h;\
	fi;\
	rm -f vertmp

Version.h: version-cache
	cp version-cache Version.h

.PHONY: version-cache

#
# Libraries
#
include Make.inc
include Eigen.inc

extra_sources+=$(WILS_FERMION_FILES)
extra_sources+=$(STAG_FERMION_FILES)
if BUILD_ZMOBIUS
  extra_sources+=$(ZWILS_FERMION_FILES)
endif
if BUILD_GPARITY
  extra_sources+=$(GP_FERMION_FILES)
endif
if BUILD_FERMION_REPS
  extra_sources+=$(ADJ_FERMION_FILES)
  extra_sources+=$(TWOIND_FERMION_FILES)
endif

lib_LIBRARIES = libGrid.a

CCFILES += $(extra_sources)
HFILES  += $(extra_headers) Config.h Version.h

libGrid_a_SOURCES              = $(CCFILES)
libGrid_adir                   = $(includedir)/Grid
nobase_dist_pkginclude_HEADERS = $(HFILES) $(eigen_files) $(eigen_unsupp_files)
