armclang 20.0 VLA

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-precision=double --enable-comms=none --enable-openmp CXX=g++-10.0.1 CC=gcc-10.0.1 CXXFLAGS="-std=c++11 -march=armv8-a+sve -msve-vector-bits=512 -fno-gcse -DA64FX -DA64FXASM -DDSLASHINTRIN" LDFLAGS=-static GRID_LDFLAGS=-static MPI_CXXLDFLAGS=-static


gcc 10.0.1 VLA

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-precision=double --enable-comms=none --enable-openmp CXX=armclang++ CC=armclang CXXFLAGS="-std=c++11 -fno-unroll-loops -mllvm -vectorizer-min-trip-count=2 -march=armv8-a+sve -DA64FX -DA64FXASM -DDSLASHINTRIN" LDFLAGS=-static GRID_LDFLAGS=-static MPI_CXXLDFLAGS=-static

should remove "-fno-strict-aliasing" for gcc 10
