# Grid [![Teamcity status](http://ci.cliath.ph.ed.ac.uk/app/rest/builds/aggregated/strob:(buildType:(affectedProject(id:GridBasedSoftware_Grid)),branch:name:develop)/statusIcon.svg)](http://ci.cliath.ph.ed.ac.uk/project.html?projectId=GridBasedSoftware_Grid&tab=projectOverview) 

**Data parallel C++ mathematical object library.**

License: GPL v2.

Last update June 2017.

_Please do not send pull requests to the `master` branch which is reserved for releases._



### Description
This library provides data parallel C++ container classes with internal memory layout
that is transformed to map efficiently to SIMD architectures. CSHIFT facilities
are provided, similar to HPF and cmfortran, and user control is given over the mapping of
array indices to both MPI tasks and SIMD processing elements.

* Identically shaped arrays then be processed with perfect data parallelisation.
* Such identically shaped arrays are called conformable arrays.

The transformation is based on the observation that Cartesian array processing involves
identical processing to be performed on different regions of the Cartesian array.

The library will both geometrically decompose into MPI tasks and across SIMD lanes.
Local vector loops are parallelised with OpenMP pragmas.

Data parallel array operations can then be specified with a SINGLE data parallel paradigm, but
optimally use MPI, OpenMP and SIMD parallelism under the hood. This is a significant simplification
for most programmers.

The layout transformations are parametrised by the SIMD vector length. This adapts according to the architecture.
Presently SSE4, ARM NEON (128 bits) AVX, AVX2, QPX (256 bits), IMCI and AVX512 (512 bits) targets are supported.

These are presented as `vRealF`, `vRealD`, `vComplexF`, and `vComplexD` internal vector data types. 
The corresponding scalar types are named `RealF`, `RealD`, `ComplexF` and `ComplexD`.

MPI, OpenMP, and SIMD parallelism are present in the library.
Please see [this paper](https://arxiv.org/abs/1512.03487) for more detail.


### Compilers

Intel ICPC v16.0.3 and later

Clang v3.5 and later (need 3.8 and later for OpenMP)

GCC   v4.9.x (recommended)

GCC   v6.3 and later

### Specific machine compilation instructions - Summit, Tesseract

[The Wiki contains specific instructions for some Summit, Tesseract and GPU compilation](https://github.com/paboyle/Grid/wiki)

### Important: 

Some versions of GCC appear to have a bug under high optimisation (-O2, -O3).

The safety of these compiler versions cannot be guaranteed at this time. Follow Issue 100 for details and updates.

GCC   v5.x

GCC   v6.1, v6.2

### Bug report

_To help us tracking and solving more efficiently issues with Grid, please report problems using the issue system of GitHub rather than sending emails to Grid developers._

When you file an issue, please go though the following checklist:

1. Check that the code is pointing to the `HEAD` of `develop` or any commit in `master` which is tagged with a version number. 
2. Give a description of the target platform (CPU, network, compiler). Please give the full CPU part description, using for example `cat /proc/cpuinfo | grep 'model name' | uniq` (Linux) or `sysctl machdep.cpu.brand_string` (macOS) and the full output the `--version` option of your compiler.
3. Give the exact `configure` command used.
4. Attach `config.log`.
5. Attach `grid.config.summary`.
6. Attach the output of `make V=1`.
7. Describe the issue and any previous attempt to solve it. If relevant, show how to reproduce the issue using a minimal working example.

### Required libraries
Grid requires:

[GMP](https://gmplib.org/), 

[MPFR](http://www.mpfr.org/) 

Bootstrapping grid downloads and uses for internal dense matrix (non-QCD operations) the Eigen library.

Grid optionally uses:

[HDF5](https://support.hdfgroup.org/HDF5/)  

[LIME](http://usqcd-software.github.io/c-lime/) for ILDG and SciDAC file format support. 

[FFTW](http://www.fftw.org) either generic version or via the Intel MKL library.

LAPACK either generic version or Intel MKL library.


### Quick start
First, start by cloning the repository:

``` bash
git clone https://github.com/paboyle/Grid.git
```

Then enter the cloned directory and set up the build system:

``` bash
cd Grid
./bootstrap.sh
```

Now you can execute the `configure` script to generate makefiles (here from a build directory):

``` bash
mkdir build; cd build
../configure --enable-simd=AVX --enable-comms=mpi-auto --prefix=<path>
```

where `--enable-simd=` set the SIMD type, `--enable-
comms=`, and `<path>` should be replaced by the prefix path where you want to
install Grid. Other options are detailed in the next section, you can also use `configure
--help` to display them. Like with any other program using GNU autotool, the
`CXX`, `CXXFLAGS`, `LDFLAGS`, ... environment variables can be modified to
customise the build.

Finally, you can build, check, and install Grid:

``` bash
make; make check; make install
```

To minimise the build time, only the tests at the root of the `tests` directory are built by default. If you want to build tests in the sub-directory `<subdir>` you can execute:

``` bash
make -C tests/<subdir> tests
```
If you want to build all the tests at once just use `make tests`.

### Build configuration options

- `--prefix=<path>`: installation prefix for Grid.
- `--with-gmp=<path>`: look for GMP in the UNIX prefix `<path>`
- `--with-mpfr=<path>`: look for MPFR in the UNIX prefix `<path>`
- `--with-fftw=<path>`: look for FFTW in the UNIX prefix `<path>`
- `--enable-lapack[=<path>]`: enable LAPACK support in Lanczos eigensolver. A UNIX prefix containing the library can be specified (optional).
- `--enable-mkl[=<path>]`: use Intel MKL for FFT (and LAPACK if enabled) routines. A UNIX prefix containing the library can be specified (optional).
- `--enable-numa`: enable NUMA first touch optimisation
- `--enable-simd=<code>`: setup Grid for the SIMD target `<code>` (default: `GEN`). A list of possible SIMD targets is detailed in a section below.
- `--enable-gen-simd-width=<size>`: select the size (in bytes) of the generic SIMD vector type (default: 32 bytes).
- `--enable-comms=<comm>`: Use `<comm>` for message passing (default: `none`). A list of possible SIMD targets is detailed in a section below.
- `--enable-rng={sitmo|ranlux48|mt19937}`: choose the RNG (default: `sitmo `).
- `--disable-timers`: disable system dependent high-resolution timers.
- `--enable-chroma`: enable Chroma regression tests.
- `--enable-doxygen-doc`: enable the Doxygen documentation generation (build with `make doxygen-doc`)

### Possible communication interfaces

The following options can be use with the `--enable-comms=` option to target different communication interfaces:

| `<comm>`       | Description                                                   |
| -------------- | ------------------------------------------------------------- |
| `none`         | no communications                                             |
| `mpi[-auto]`   | MPI communications                                            |
| `mpi3[-auto]`  | MPI communications using MPI 3 shared memory                  |
| `shmem `       | Cray SHMEM communications                                     |

For the MPI interfaces the optional `-auto` suffix instructs the `configure` scripts to determine all the necessary compilation and linking flags. This is done by extracting the informations from the MPI wrapper specified in the environment variable `MPICXX` (if not specified `configure` will scan though a list of default names). The `-auto` suffix is not supported by the Cray environment wrapper scripts. Use the standard versions instead.  

### Possible SIMD types

The following options can be use with the `--enable-simd=` option to target different SIMD instruction sets:

| `<code>`    | Description                            |
| ----------- | -------------------------------------- |
| `GEN`       | generic portable vector code           |
| `SSE4`      | SSE 4.2 (128 bit)                      |
| `AVX`       | AVX (256 bit)                          |
| `AVXFMA`    | AVX (256 bit) + FMA                    |
| `AVXFMA4`   | AVX (256 bit) + FMA4                   |
| `AVX2`      | AVX 2 (256 bit)                        |
| `AVX512`    | AVX 512 bit                            |
| `NEONv8`    | [ARM NEON](http://infocenter.arm.com/help/index.jsp?topic=/com.arm.doc.den0024a/ch07s03.html) (128 bit)                     |
| `QPX`       | IBM QPX (256 bit)                      |

Alternatively, some CPU codenames can be directly used:

| `<code>`    | Description                            |
| ----------- | -------------------------------------- |
| `KNL`       | [Intel Xeon Phi codename Knights Landing](http://ark.intel.com/products/codename/48999/Knights-Landing) |
| `SKL`       | [Intel Skylake with AVX512 extensions](https://ark.intel.com/products/codename/37572/Skylake#@server) |
| `BGQ`       | Blue Gene/Q                            |

#### Notes:
- We currently support AVX512 for the Intel compiler and GCC (KNL and SKL target). Support for clang will appear in future versions of Grid when the AVX512 support in the compiler will be more advanced.
- For BG/Q only [bgclang](http://trac.alcf.anl.gov/projects/llvm-bgq) is supported. We do not presently plan to support more compilers for this platform.
- BG/Q performances are currently rather poor. This is being investigated for future versions.
- The vector size for the `GEN` target can be specified with the `configure` script option `--enable-gen-simd-width`.

### Build setup for Intel Knights Landing platform

The following configuration is recommended for the Intel Knights Landing platform:

``` bash
../configure --enable-simd=KNL        \
             --enable-comms=mpi-auto  \
             --enable-mkl             \
             CXX=icpc MPICXX=mpiicpc
```
The MKL flag enables use of BLAS and FFTW from the Intel Math Kernels Library.

If you are working on a Cray machine that does not use the `mpiicpc` wrapper, please use:

``` bash
../configure --enable-simd=KNL        \
             --enable-comms=mpi       \
             --enable-mkl             \
             CXX=CC CC=cc
```

If gmp and mpfr are NOT in standard places (/usr/) these flags may be needed:
``` bash
               --with-gmp=<path>        \
               --with-mpfr=<path>       \
```
where `<path>` is the UNIX prefix where GMP and MPFR are installed. 

Knight's Landing with Intel Omnipath adapters with two adapters per node 
presently performs better with use of more than one rank per node, using shared memory 
for interior communication. This is the mpi3 communications implementation. 
We recommend four ranks per node for best performance, but optimum is local volume dependent.

``` bash
../configure --enable-simd=KNL        \
             --enable-comms=mpi3-auto \
             --enable-mkl             \
             CC=icpc MPICXX=mpiicpc 
```

### Build setup for Intel Haswell Xeon platform

The following configuration is recommended for the Intel Haswell platform:

``` bash
../configure --enable-simd=AVX2       \
             --enable-comms=mpi3-auto \
             --enable-mkl             \
             CXX=icpc MPICXX=mpiicpc
```
The MKL flag enables use of BLAS and FFTW from the Intel Math Kernels Library.

If gmp and mpfr are NOT in standard places (/usr/) these flags may be needed:
``` bash
               --with-gmp=<path>        \
               --with-mpfr=<path>       \
```
where `<path>` is the UNIX prefix where GMP and MPFR are installed. 

If you are working on a Cray machine that does not use the `mpiicpc` wrapper, please use:

``` bash
../configure --enable-simd=AVX2       \
             --enable-comms=mpi3      \
             --enable-mkl             \
             CXX=CC CC=cc
```
Since Dual socket nodes are commonplace, we recommend MPI-3 as the default with the use of 
one rank per socket. If using the Intel MPI library, threads should be pinned to NUMA domains using
```
        export I_MPI_PIN=1
```
This is the default.

### Build setup for Intel Skylake Xeon platform

The following configuration is recommended for the Intel Skylake platform:

``` bash
../configure --enable-simd=AVX512     \
             --enable-comms=mpi3      \
             --enable-mkl             \
             CXX=mpiicpc
```
The MKL flag enables use of BLAS and FFTW from the Intel Math Kernels Library.

If gmp and mpfr are NOT in standard places (/usr/) these flags may be needed:
``` bash
               --with-gmp=<path>        \
               --with-mpfr=<path>       \
```
where `<path>` is the UNIX prefix where GMP and MPFR are installed. 

If you are working on a Cray machine that does not use the `mpiicpc` wrapper, please use:

``` bash
../configure --enable-simd=AVX512     \
             --enable-comms=mpi3      \
             --enable-mkl             \
             CXX=CC CC=cc
```
Since Dual socket nodes are commonplace, we recommend MPI-3 as the default with the use of 
one rank per socket. If using the Intel MPI library, threads should be pinned to NUMA domains using
``` 
        export I_MPI_PIN=1
```
This is the default. 

#### Expected Skylake Gold 6148 dual socket (single prec, single node 20+20 cores) performance using NUMA MPI mapping): 

mpirun -n 2 benchmarks/Benchmark_dwf --grid 16.16.16.16 --mpi 2.1.1.1 --cacheblocking 2.2.2.2 --dslash-asm --shm 1024 --threads 18 

TBA


### Build setup for AMD EPYC / RYZEN

The AMD EPYC is a multichip module comprising 32 cores spread over four distinct chips each with 8 cores.
So, even with a single socket node there is a quad-chip module. Dual socket nodes with 64 cores total
are common. Each chip within the module exposes a separate NUMA domain.
There are four NUMA domains per socket and we recommend one MPI rank per NUMA domain.
MPI-3 is recommended with the use of four ranks per socket,
and 8 threads per rank. 

The following configuration is recommended for the AMD EPYC platform.

``` bash
../configure --enable-simd=AVX2       \
             --enable-comms=mpi3 \
             CXX=mpicxx 
```

If gmp and mpfr are NOT in standard places (/usr/) these flags may be needed:
``` bash
               --with-gmp=<path>        \
               --with-mpfr=<path>       \
```
where `<path>` is the UNIX prefix where GMP and MPFR are installed. 

Using MPICH and g++ v4.9.2, best performance can be obtained using explicit GOMP_CPU_AFFINITY flags for each MPI rank.
This can be done by invoking MPI on a wrapper script omp_bind.sh to handle this. 

It is recommended to run 8 MPI ranks on a single dual socket AMD EPYC, with 8 threads per rank using MPI3 and
shared memory to communicate within this node:

mpirun -np 8 ./omp_bind.sh ./Benchmark_dwf --mpi 2.2.2.1 --dslash-unroll --threads 8 --grid 16.16.16.16 --cacheblocking 4.4.4.4 

Where omp_bind.sh does the following:
```
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

Performance:

#### Expected AMD EPYC 7601 dual socket (single prec, single node 32+32 cores) performance using NUMA MPI mapping): 

mpirun  -np 8 ./omp_bind.sh ./Benchmark_dwf --threads 8 --mpi 2.2.2.1 --dslash-unroll --grid 16.16.16.16 --cacheblocking 4.4.4.4

TBA

### Build setup for BlueGene/Q

To be written...

### Build setup for ARM Neon

To be written...

### Build setup for laptops, other compilers, non-cluster builds

Many versions of g++ and clang++ work with Grid, and involve merely replacing CXX (and MPICXX),
and omit the enable-mkl flag. 

Single node builds are enabled with 
```
            --enable-comms=none
```

FFTW support that is not in the default search path may then enabled with
```
    --with-fftw=<installpath>
```

BLAS will not be compiled in by default, and Lanczos will default to Eigen diagonalisation.

