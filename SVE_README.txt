gcc 10.0.1 VLA

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-precision=double --enable-comms=none --enable-openmp CXX=g++-10.0.1 CC=gcc-10.0.1 CXXFLAGS="-std=c++11 -march=armv8-a+sve -msve-vector-bits=512 -fno-gcse -DA64FX -DA64FXASM -DDSLASHINTRIN" LDFLAGS=-static GRID_LDFLAGS=-static MPI_CXXLDFLAGS=-static


armclang 20.0 VLA

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-precision=double --enable-comms=none --enable-openmp CXX=armclang++ CC=armclang CXXFLAGS="-std=c++11 -fno-unroll-loops -mllvm -vectorizer-min-trip-count=2 -march=armv8-a+sve -DARMCLANGCOMPAT -DA64FX -DA64FXASM -DDSLASHINTRIN" LDFLAGS=-static GRID_LDFLAGS=-static MPI_CXXLDFLAGS=-static

must use armclang 20.0 with ARMCLANGCOMPAT applied, otherwise Benchmark_wilson gives wrong result


armclang 20.1 VLA

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-precision=double --enable-comms=none --enable-openmp CXX=armclang++ CC=armclang CXXFLAGS="-std=c++11 -mcpu=a64fx -DA64FX -DA64FXASM -DDSLASHINTRIN" LDFLAGS=-static GRID_LDFLAGS=-static MPI_CXXLDFLAGS=-static

must use armclang 20.1 with ARMCLANGCOMPAT applied, otherwise Benchmark_wilson gives wrong result
Test_simd build error caused by -mcpu=a64fx ?



Fujitsu FCC

../configure --with-lime=$HOME/grid-a64fx/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-precision=double --enable-comms=none --enable-openmp --with-mpfr=/home/users/gre/gre-1/grid-a64fx/mpfr-build/install CXX=FCC CC=fcc CXXFLAGS="-Nclang -Kfast -DA64FX -DA64FXASM -DDSLASHINTRIN"





what about "-fno-strict-aliasing" in general?

