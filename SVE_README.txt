* gcc 10.1 prebuild, QPACE4 interactive login w/ MPI

scl enable gcc-toolset-10 bash
module load mpi/openmpi-aarch64

../configure --enable-simd=A64FX --enable-comms=mpi3 --enable-shm=shmget CXX=mpicxx CC=mpicc


================================== deprecated ================================================

* gcc 10.1 prebuild, QPACE4 interactive login

scl enable gcc-toolset-10 bash

../configure --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=none --enable-openmp CXX=g++ CC=gcc CXXFLAGS="-std=c++11 -march=armv8-a+sve -msve-vector-bits=512 -fno-gcse -DA64FXFIXEDSIZE -DA64FXASM -DDSLASHINTRIN"

* gcc 10.1 prebuild w/ MPI, QPACE4 interactive login

scl enable gcc-toolset-10 bash
module load mpi/openmpi-aarch64

../configure --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=mpi-auto --enable-shm=shmget --enable-openmp CXX=mpicxx CC=mpicc CXXFLAGS="-std=c++11 -march=armv8-a+sve -msve-vector-bits=512 -fno-gcse -DA64FXFIXEDSIZE -DA64FXASM -DDSLASHINTRIN"

------------------------------------------------------------------------------

* armclang 20.2 (qp4)

../configure --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=none --enable-openmp CXX=armclang++ CC=armclang CXXFLAGS="-std=c++11 -mcpu=a64fx -DA64FX -DARMCLANGCOMPAT -DA64FXASM -DDSLASHINTRIN"

------------------------------------------------------------------------------

* gcc 10.0.1 VLA (merlin)

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=none --enable-openmp CXX=g++-10.0.1 CC=gcc-10.0.1 CXXFLAGS="-std=c++11 -march=armv8-a+sve -msve-vector-bits=512 -fno-gcse -DA64FX -DA64FXASM -DDSLASHINTRIN" LDFLAGS=-static GRID_LDFLAGS=-static MPI_CXXLDFLAGS=-static


* gcc 10.0.1 fixed-size ACLE (merlin)

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=none --enable-openmp CXX=g++-10.0.1 CC=gcc-10.0.1 CXXFLAGS="-std=c++11 -march=armv8-a+sve -msve-vector-bits=512 -fno-gcse -DA64FXFIXEDSIZE -DA64FXASM -DDSLASHINTRIN"


* gcc 10.0.1 fixed-size ACLE (fjt) w/ MPI

export OMPI_CC=gcc-10.0.1
export OMPI_CXX=g++-10.0.1
export MPICH_CC=gcc-10.0.1
export MPICH_CXX=g++-10.0.1

$ ../configure --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=mpi3 --enable-openmp CXX=mpiFCC CC=mpifcc CXXFLAGS="-std=c++11 -march=armv8-a+sve -msve-vector-bits=512 -fno-gcse -DA64FXFIXEDSIZE -DA64FXASM -DDSLASHINTRIN -DTOFU -I/opt/FJSVxtclanga/tcsds-1.2.25/include/mpi/fujitsu -lrt" LDFLAGS="-L/opt/FJSVxtclanga/tcsds-1.2.25/lib64 -lrt"

--------------------------------------------------------

* armclang 20.0 VLA (merlin)

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=none --enable-openmp CXX=armclang++ CC=armclang CXXFLAGS="-std=c++11 -fno-unroll-loops -mllvm -vectorizer-min-trip-count=2 -march=armv8-a+sve -DARMCLANGCOMPAT -DA64FX -DA64FXASM -DDSLASHINTRIN" LDFLAGS=-static GRID_LDFLAGS=-static MPI_CXXLDFLAGS=-static

TODO check ARMCLANGCOMPAT


* armclang 20.1 VLA (merlin)

../configure --with-lime=/home/men04359/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=none --enable-openmp CXX=armclang++ CC=armclang CXXFLAGS="-std=c++11 -mcpu=a64fx -DARMCLANGCOMPAT -DA64FX -DA64FXASM -DDSLASHINTRIN" LDFLAGS=-static GRID_LDFLAGS=-static MPI_CXXLDFLAGS=-static

TODO check ARMCLANGCOMPAT


* armclang 20.1 VLA (fjt cluster)

../configure --with-lime=$HOME/local --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=none --enable-openmp CXX=armclang++ CC=armclang CXXFLAGS="-std=c++11 -mcpu=a64fx -DARMCLANGCOMPAT -DA64FX -DA64FXASM -DDSLASHINTRIN -DTOFU"

TODO check ARMCLANGCOMPAT


* armclang 20.1 VLA w/MPI (fjt cluster)

../configure --with-lime=$HOME/local --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=mpi3 --enable-openmp CXX=mpiFCC CC=mpifcc CXXFLAGS="-std=c++11 -mcpu=a64fx -DA64FX -DA64FXASM -DDSLASHINTRIN -DTOFU -I/opt/FJSVxtclanga/tcsds-1.2.25/include/mpi/fujitsu -lrt" LDFLAGS="-L/opt/FJSVxtclanga/tcsds-1.2.25/lib64"

No ARMCLANGCOMPAT -> still correct ?

--------------------------------------------------------

* Fujitsu fcc

../configure --with-lime=$HOME/grid-a64fx/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=none --enable-openmp --with-mpfr=/home/users/gre/gre-1/grid-a64fx/mpfr-build/install CXX=FCC CC=fcc CXXFLAGS="-Nclang -Kfast -DA64FX -DA64FXASM -DDSLASHINTRIN"


* Fujitsu fcc w/ MPI

../configure --with-lime=$HOME/grid-a64fx/lime/c-lime --without-hdf5 --enable-gen-simd-width=64 --enable-simd=GEN --enable-comms=mpi --enable-openmp --with-mpfr=/home/users/gre/gre-1/grid-a64fx/mpfr-build/install CXX=mpiFCC CC=mpifcc CXXFLAGS="-Nclang -Kfast -DA64FX -DA64FXASM -DDSLASHINTRIN -DTOFU"
