# Grid
<table>
<tr>
    <td>Last stable release</td>
    <td><a href="https://travis-ci.org/paboyle/Grid">
    <img src="https://travis-ci.org/paboyle/Grid.svg?branch=master"></a>
    </td>
</tr>
<tr>
    <td>Development branch</td>
    <td><a href="https://travis-ci.org/paboyle/Grid">
    <img src="https://travis-ci.org/paboyle/Grid.svg?branch=develop"></a>
    </td>
</tr>
</table>

**Data parallel C++ mathematical object library.**

License: GPL v2.

Last update Nov 2016.

_Please do not send pull requests to the `master` branch which is reserved for releases._

### Bug report

_To help us tracking and solving more efficiently issues with Grid, please report problems using the issue system of GitHub rather than sending emails to Grid developers._

When you file an issue, please go though the following checklist:

1. Check that the code is pointing to the `HEAD` of `develop` or any commit in `master` which is tagged with a version number. 
2. Give a description of the target platform (CPU, network, compiler). Please give the full CPU part description, using for example `cat /proc/cpuinfo | grep 'model name' | uniq` (Linux) or `sysctl machdep.cpu.brand_string` (macOS) and the full output the `--version` option of your compiler.
3. Give the exact `configure` command used.
4. Attach `config.log`.
5. Attach `config.summary`.
6. Attach the output of `make V=1`.
7. Describe the issue and any previous attempt to solve it. If relevant, show how to reproduce the issue using a minimal working example.



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
Presently SSE4 (128 bit) AVX, AVX2, QPX (256 bit), IMCI, and AVX512 (512 bit) targets are supported (ARM NEON on the way).

These are presented as `vRealF`, `vRealD`, `vComplexF`, and `vComplexD` internal vector data types. These may be useful in themselves for other programmers.
The corresponding scalar types are named `RealF`, `RealD`, `ComplexF` and `ComplexD`.

MPI, OpenMP, and SIMD parallelism are present in the library.
Please see https://arxiv.org/abs/1512.03487 for more detail.

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
../configure --enable-precision=double --enable-simd=AVX --enable-comms=mpi-auto --prefix=<path>
```

where `--enable-precision=` set the default precision,
`--enable-simd=` set the SIMD type, `--enable-
comms=`, and `<path>` should be replaced by the prefix path where you want to
install Grid. Other options are detailed in the next section, you can also use `configure
--help` to display them. Like with any other program using GNU autotool, the
`CXX`, `CXXFLAGS`, `LDFLAGS`, ... environment variables can be modified to
customise the build.

Finally, you can build and install Grid:

``` bash
make; make install
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
- `--enable-precision={single|double}`: set the default precision (default: `double`).
- `--enable-precision=<comm>`: Use `<comm>` for message passing (default: `none`). A list of possible SIMD targets is detailed in a section below.
- `--enable-rng={ranlux48|mt19937}`: choose the RNG (default: `ranlux48 `).
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
| `mpi3l[-auto]` | MPI communications using MPI 3 shared memory and leader model |
| `shmem `       | Cray SHMEM communications                                     |

For the MPI interfaces the optional `-auto` suffix instructs the `configure` scripts to determine all the necessary compilation and linking flags. This is done by extracting the informations from the MPI wrapper specified in the environment variable `MPICXX` (if not specified `configure` will scan though a list of default names).

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
| `QPX`       | QPX (256 bit)                          |

Alternatively, some CPU codenames can be directly used:

| `<code>`    | Description                            |
| ----------- | -------------------------------------- |
| `KNC`       | [Intel Xeon Phi codename Knights Corner](http://ark.intel.com/products/codename/57721/Knights-Corner) |
| `KNL`       | [Intel Xeon Phi codename Knights Landing](http://ark.intel.com/products/codename/48999/Knights-Landing) |
| `BGQ`       | Blue Gene/Q                            |

#### Notes:
- We currently support AVX512 only for the Intel compiler. Support for GCC and clang will appear in future versions of Grid when the AVX512 support within GCC and clang will be more advanced.
- For BG/Q only [bgclang](http://trac.alcf.anl.gov/projects/llvm-bgq) is supported. We do not presently plan to support more compilers for this platform.
- BG/Q performances are currently rather poor. This is being investigated for future versions.
- The vector size for the `GEN` target can be specified with the `configure` script option `--enable-gen-simd-width`.

### Build setup for Intel Knights Landing platform

The following configuration is recommended for the Intel Knights Landing platform:

``` bash
../configure --enable-precision=double\
             --enable-simd=KNL        \
             --enable-comms=mpi-auto \
             --with-gmp=<path>        \
             --with-mpfr=<path>       \
             --enable-mkl             \
             CXX=icpc MPICXX=mpiicpc
```

where `<path>` is the UNIX prefix where GMP and MPFR are installed. If you are working on a Cray machine that does not use the `mpiicpc` wrapper, please use:

``` bash
../configure --enable-precision=double\
             --enable-simd=KNL        \
             --enable-comms=mpi       \
             --with-gmp=<path>        \
             --with-mpfr=<path>       \
             --enable-mkl             \
             CXX=CC CC=cc
```