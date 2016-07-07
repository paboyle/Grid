
bin_PROGRAMS += Test_nersc_io Test_serialisation


Test_nersc_io_SOURCES=Test_nersc_io.cc
Test_nersc_io_LDADD=-lGrid


Test_serialisation_SOURCES=Test_serialisation.cc
Test_serialisation_LDADD=-lGrid

