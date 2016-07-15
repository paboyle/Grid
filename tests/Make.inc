
bin_PROGRAMS += Test_cshift Test_simd Test_stencil


Test_cshift_SOURCES=Test_cshift.cc
Test_cshift_LDADD=-lGrid


Test_simd_SOURCES=Test_simd.cc
Test_simd_LDADD=-lGrid


Test_stencil_SOURCES=Test_stencil.cc
Test_stencil_LDADD=-lGrid

