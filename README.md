# Grid [![Build Status](https://travis-ci.org/paboyle/Grid.svg?branch=master)](https://travis-ci.org/paboyle/Grid)
Data parallel C++ mathematical object library

Last update 2015/7/30

This library provides data parallel C++ container classes with internal memory layout
that is transformed to map efficiently to SIMD architectures. CSHIFT facilities
are provided, similar to HPF and cmfortran, and user control is given over the mapping of
array indices to both MPI tasks and SIMD processing elements.

* Identically shaped arrays then be processed with perfect data parallelisation.
* Such identically shapped arrays are called conformable arrays.

The transformation is based on the observation that Cartesian array processing involves
identical processing to be performed on different regions of the Cartesian array.

The library will both geometrically decompose into MPI tasks and across SIMD lanes.
Local vector loops are parallelised with OpenMP pragmas.

Data parallel array operations can then be specified with a SINGLE data parallel paradigm, but
optimally use MPI, OpenMP and SIMD parallelism under the hood. This is a significant simplification
for most programmers.

The layout transformations are parametrised by the SIMD vector length. This adapts according to the architecture.
Presently SSE4 (128 bit) AVX, AVX2 (256 bit) and IMCI and AVX512 (512 bit) targets are supported (ARM NEON on the way).

These are presented as 

     vRealF, vRealD, vComplexF, vComplexD 

internal vector data types. These may be useful in themselves for other programmers.
The corresponding scalar types are named

     RealF, RealD, ComplexF, ComplexD

MPI, OpenMP, and SIMD parallelism are present in the library.

   You can give `configure' initial values for configuration parameters
by setting variables in the command line or in the environment.  Here
are examples:

     ./configure CXX=clang++ CXXFLAGS="-std=c++11 -O3 -msse4" --enable-simd=SSE4

     ./configure CXX=clang++ CXXFLAGS="-std=c++11 -O3 -mavx" --enable-simd=AVX

     ./configure CXX=clang++ CXXFLAGS="-std=c++11 -O3 -mavx2" --enable-simd=AVX2

     ./configure CXX=icpc CXXFLAGS="-std=c++11 -O3 -mmic" --enable-simd=AVX512 --host=none
     
Note: Before running configure it could be necessary to execute the script 
       
       script/filelist


     
For developers:
Use reconfigure_script in the scripts/ directory to create the autotools environment 

