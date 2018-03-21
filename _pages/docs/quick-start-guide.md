---
title : "Documentation"
author_profile: false
excerpt: "Quick start guide"
header:
  overlay_color: "#5DADE2"
permalink: /docs/quick-start-guide/
sidebar:
  nav : docs
---
{% include toc icon="gears" title="Quick-start" %}

Please send all pull requests to the `develop` branch.
{: .notice--warning}

## Requirements

### Required libraries

* [GMP](https://gmplib.org/) is the GNU Multiple Precision Library (RHMC support).
* [MPFR](http://www.mpfr.org/) is a C library for multiple-precision floating-point computations with correct rounding (RHMC support). 
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): bootstrapping GRID downloads and uses for internal dense matrix (non-QCD operations) the Eigen library.

Grid optionally uses:

* [HDF5](https://support.hdfgroup.org/HDF5/) for structured data I/O

* [LIME](http://usqcd-software.github.io/c-lime/) for ILDG and SciDAC file format support.

* [FFTW](http://www.fftw.org/) either generic version or via the Intel MKL library.

* [LAPACK](http://www.netlib.org/lapack/) either generic version or Intel MKL library.


### Compilers

* Intel ICPC v17 and later

* Clang v3.5 and later (need 3.8 and later for OpenMP)

* GCC v4.9.x 

* GCC v6.3 and later (recommended)

**Important:**

Some versions of GCC appear to have a bug under high optimisation (-O2, -O3).

The safety of these compiler versions cannot be guaranteed at this time. Follow [Issue 100](https://github.com/paboyle/Grid/issues/100) for details and updates.

GCC v5.x, v6.1, v6.2




## Quick start
First, start by cloning the repository:

``` bash
git clone https://github.com/paboyle/Grid.git
```

Then enter the cloned directory and set up the build system:

```bash
cd Grid
./bootstrap.sh
```

Now you can execute the `configure` script to generate makefiles as in this example (here from a build directory):

``` bash
mkdir build
cd build
../configure --enable-precision=double --enable-simd=AVX --enable-comms=mpi-auto --prefix=<path>
```

where:

``` bash
  --enable-precision=single|double
```

sets the **default precision**. Since this is largely a benchmarking convenience, it is anticipated that the default precision may be removed in future implementations, and that explicit type selection be made at all points. Naturally, most code will be type templated in any case.::

``` bash
   --enable-simd=GEN|SSE4|AVX|AVXFMA|AVXFMA4|AVX2|AVX512|NEONv8|QPX
```

sets the **SIMD architecture**, 

``` bash
   --enable-comms=mpi|none
```

selects whether to use MPI communication (mpi) or no communication (none). 

``` bash
   --prefix=<path>
```

should be passed the prefix path where you want to install Grid. 

Other options are detailed in the next section, you can also use 

```bash
   configure --help
```

to display them. 

Like with any other program using GNU autotool, the 

```bash
   CXX, CXXFLAGS, LDFLAGS, ... 
```

environment variables can be modified to customise the build.

Finally, you can build and install Grid:

``` bash
make
make install   #this is optional
```

To minimise the build time, only the tests at the root of the `tests` directory are built by default. If you want to build tests in the sub-directory `<subdir>` you can execute:

``` bash
make -C tests/<subdir> tests
```
If you want to build all the tests at once just use `make tests`.

## Build configuration options

A full list of configurations options is available with the `./configure --help` command: 

* `--prefix=<path>`: installation prefix for Grid.

* `--with-gmp=<path>`: look for GMP in the UNIX prefix `<path>`

* `--with-mpfr=<path>`: look for MPFR in the UNIX prefix `<path>`

* `--with-fftw=<path>`: look for FFTW in the UNIX prefix `<path>`

* `--with-hdf5=<path>`: look for HDF5 in the UNIX prefix `<path>`

* `--with-lime=<path>`: look for the C-LIME library in the UNIX prefix `<path>`

* `--enable-sfw-fp16=<yes|no`: Enable software FP16 communications support (default `yes`)

* `--enable-lapack[=<path>]`: enable LAPACK support in Lanczos eigensolver. A UNIX prefix containing the library can be specified (optional).

* `--enable-mkl[=<path>]`: use Intel MKL for FFT (and LAPACK if enabled) routines. A UNIX prefix containing the library can be specified (optional).

* `--enable-numa`: enable [numa first touch policy](http://queue.acm.org/detail.cfm?id=2513149) optimization (default `no`)

* `--enable-simd=<code>`: setup Grid for the SIMD target `<code>` (default: `GEN`). [List of possible SIMD targets](/Grid/docs/simd_targets/).

* `--enable-gen-simd-width=<size>`: select the size (in bytes) of the generic SIMD vector type (default: 32 bytes).

* `--enable-precision={single|double}`: set the default precision (default: `double`).

* `--enable-comms=<comm>`: Use `<comm>` for message passing (default: `none`). [List of possible comm targets](/Grid/docs/comm_interfaces/). 

* `--enable-shm=<shm>`: Use `<shm>` for shared memory behaviour (default: `shmopen`). [List of possible shm targets](/Grid/docs/comm_interfaces/). 

* `--enable-shmpath=<path>`: Select `<path>` for the shared memory mmap base path for libhugetlbfs. 

* `--enable-rng={sitmo|ranlux48|mt19937}` choose the RNG (default: `sitmo`).

* `--disable-timers`: disable system dependent high-resolution timers.

* `--enable-chroma`: enable Chroma regression tests. A compiled version of Chroma is assumed to be present. 

* `--enable-doxygen-doc`: enable the Doxygen documentation generation (build with `make doxygen-doc`)


{% if site.option=="web" %}
More details on the *Getting started* menu entries on the left. 
{% endif %}

This document was updated on March 2018. 
{: .notice}

