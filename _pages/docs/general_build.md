---
title : "Documentation"
author_profile: false
excerpt: "Building on Intel and AMD targets"
header:
  overlay_color: "#5DADE2"
permalink: /docs/general_build/
sidebar:
  nav : docs
---
{% include base_path %}
The information included in this page has been updated on *March 2018* and it is valid for the [release version 0.7.0](https://github.com/paboyle/Grid/tree/release/v0.7.0).

{% include toc icon="gears" title="Contents" %}

### Building for the Intel Knights Landing

The following configuration is recommended for the [Intel Knights Landing](http://ark.intel.com/products/codename/48999/Knights-Landing) platform:

``` text
../configure --enable-precision=double\
             --enable-simd=KNL        \
             --enable-comms=mpi-auto  \
             --with-gmp=<path>        \
             --with-mpfr=<path>       \
             --enable-mkl             \
             CXX=icpc MPICXX=mpiicpc
```

where `<path>` is the UNIX prefix where GMP and MPFR are installed. If you are working on a Cray machine that does not use the `mpiicpc` wrapper, please use:

``` text
../configure --enable-precision=double\
             --enable-simd=KNL        \
             --enable-comms=mpi       \
             --with-gmp=<path>        \
             --with-mpfr=<path>       \
             --enable-mkl             \
             CXX=CC CC=cc
```

### Building for the Intel Haswell


The following configuration is recommended for the [Intel Haswell platform](https://ark.intel.com/products/codename/42174/Haswell):

```bash
  ../configure --enable-precision=double\
             --enable-simd=AVX2       \
             --enable-comms=mpi-auto \
             --enable-mkl             \
             CXX=icpc MPICXX=mpiicpc
```

The MKL flag enables use of BLAS and FFTW from the Intel Math Kernels Library.

If gmp and mpfr are NOT in standard places (`/usr/`) these flags may be needed:

```bash
               --with-gmp=<path>        \
               --with-mpfr=<path>       
```

where `<path>` is the UNIX prefix where GMP and MPFR are installed. 

If you are working on a Cray machine that does not use the `mpiicpc` wrapper, please use:

```bash
  ../configure --enable-precision=double\
             --enable-simd=AVX2       \
             --enable-comms=mpi      \
             --enable-mkl             \
             CXX=CC CC=cc
```

If using the Intel MPI library, threads should be pinned to NUMA domains using:

```bash
        export I_MPI_PIN=1
```
This is the default.

### Building for the Intel Skylake

The following configuration is recommended for the [Intel Skylake platform](https://ark.intel.com/products/codename/37572/Skylake):

```bash
  ../configure --enable-precision=double\
             --enable-simd=AVX512     \
             --enable-comms=mpi-auto  \
             --enable-mkl             \
             CXX=mpiicpc
```

The MKL flag enables use of BLAS and FFTW from the Intel Math Kernels Library.

If gmp and mpfr are NOT in standard places (`/usr/`) these flags may be needed:

```bash
               --with-gmp=<path>        \
               --with-mpfr=<path>       \
```

where `<path>` is the UNIX prefix where GMP and MPFR are installed. 

If you are working on a Cray machine that does not use the `mpiicpc` wrapper, please use:

``` bash
  ../configure --enable-precision=double\
             --enable-simd=AVX512     \
             --enable-comms=mpi       \
             --enable-mkl             \
             CXX=CC CC=cc
```

If using the Intel MPI library, threads should be pinned to NUMA domains using:

```bash
        export I_MPI_PIN=1
```

This is the default. 

### Building for the AMD Epyc

The [AMD EPYC](https://www.amd.com/en/products/epyc) is a multichip module comprising 32 cores spread over four distinct chips each with 8 cores.
So, even with a single socket node there is a quad-chip module. Dual socket nodes with 64 cores total
are common. Each chip within the module exposes a separate NUMA domain.
There are four NUMA domains per socket and we recommend one MPI rank per NUMA domain.
MPI-3 is recommended with the use of four ranks per socket,
and 8 threads per rank. 

The following configuration is recommended for the AMD EPYC platform:

```bash
  ../configure --enable-precision=double\
             --enable-simd=AVX2       \
             --enable-comms=mpi3 \
             CXX=mpicxx 
```

If `gmp` and `mpfr` are NOT in standard places (`/usr/`) these flags may be needed::

```bash
               --with-gmp=<path>        \
               --with-mpfr=<path>       
```

where `<path>` is the UNIX prefix where GMP and MPFR are installed. 

Using MPICH and g++ v4.9.2, best performance can be obtained using explicit GOMP_CPU_AFFINITY flags for each MPI rank.
This can be done by invoking MPI on a wrapper script omp_bind.sh to handle this. 

It is recommended to run 8 MPI ranks on a single dual socket AMD EPYC, with 8 threads per rank using MPI3 and
shared memory to communicate within this node:

```bash
  mpirun -np 8 ./omp_bind.sh ./Benchmark_dwf --mpi 2.2.2.1 --dslash-unroll --threads 8 --grid 16.16.16.16 --cacheblocking 4.4.4.4 
```

Where omp_bind.sh does the following:

```bash
  #!/bin/bash

  numanode=` expr $PMI_RANK % 8 `
  basecore=`expr $numanode \* 16`
  core0=`expr $basecore + 0 `
  core1=`expr $basecore + 2 `
  core2=`expr $basecore + 4 `
  core3=`expr $basecore + 6 `
  core4=`expr $basecore + 8 `
  core5=`expr $basecore + 10 `
  core6=`expr $basecore + 12 `
  core7=`expr $basecore + 14 `

  export GOMP_CPU_AFFINITY="$core0 $core1 $core2 $core3 $core4 $core5 $core6 $core7"
  echo GOMP_CUP_AFFINITY $GOMP_CPU_AFFINITY

  $@
```

### Build setup for laptops, other compilers, non-cluster builds

Many versions of `g++` and `clang++` work with Grid, and involve merely replacing `CXX` (and `MPICXX`),
and omit the `enable-mkl` flag. 

Single node, non MPI builds are enabled with:

```bash
  --enable-comms=none
```

FFTW support that is not in the default search path may then enabled with:

```bash
  --with-fftw=<installpath>
```

BLAS will not be compiled in by default, and Lanczos will default to Eigen diagonalisation.



### Notes

- [GMP](https://gmplib.org/) is the GNU Multiple Precision Library.
- [MPFR](http://www.mpfr.org/) is a C library for multiple-precision floating-point computations with correct rounding.
- Both libaries are necessary for the RHMC support. 




{% include paginator.html %}