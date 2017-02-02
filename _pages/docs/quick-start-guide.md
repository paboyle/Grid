---
layout: single
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

## Quick start
First, start by cloning the repository:

``` bash
git clone https://github.com/paboyle/Grid.git
```

Then enter the cloned directory and set up the build system:

``` bash
cd Grid
./bootstrap.sh
```

Now you can execute the `configure` script to generate makefiles as in this example (here from a build directory):

``` bash
mkdir build
cd build
../configure --enable-precision=double --enable-simd=AVX --enable-comms=mpi-auto --prefix=<path>
```

where `--enable-precision=` sets the default precision,
`--enable-simd=` sets the SIMD type, `--enable-
comms=`, and `<path>` should be replaced by the prefix path where you want to
install Grid (optional). Other options are detailed in the next section, you can also use `configure
--help` to display them. Like with any other program using GNU autotool, the
`CXX`, `CXXFLAGS`, `LDFLAGS`, ... environment variables can be modified to
customise the build.

Finally, you can build and install Grid:

``` bash
make
make install   #this is optional
```

To minimise the build time, only the tests at the root of the `tests` directory are built by default. If you want to build tests in the sub-directory `<subdir>` you can execute:

``` bash
make -C tests/<subdir> tests
```

## Build configuration options

A full list of configurations options is available with the `./configure --help` command: 

Here we report the more common ones. 

- `--prefix=<path>`: installation prefix for Grid.
- `--with-gmp=<path>`: look for GMP in the UNIX prefix `<path>`
- `--with-mpfr=<path>`: look for MPFR in the UNIX prefix `<path>`
- `--with-fftw=<path>`: look for FFTW in the UNIX prefix `<path>`
- `--enable-lapack[=<path>]`: enable LAPACK support in Lanczos eigensolver. A UNIX prefix containing the library can be specified (optional).
- `--enable-mkl[=<path>]`: use Intel MKL for FFT (and LAPACK if enabled) routines. A UNIX prefix containing the library can be specified (optional).
- `--enable-numa`: enable [numa first touch policy](http://queue.acm.org/detail.cfm?id=2513149) optimization (default `no`)
- `--enable-simd=<code>`: setup Grid for the SIMD target `<code>` (default: `GEN`). [List of possible SIMD targets](/Grid/docs/simd_targets/).
- `--enable-precision={single|double}`: set the default precision (default: `double`).
- `--enable-precision=<comm>`: Use `<comm>` for message passing (default: `none`). [List of possible comm targets](/Grid/docs/comm_interfaces/). 
- `--enable-rng={ranlux48|mt19937|sitmo}`: choose the RNG (default: `ranlux48`).
- `--disable-timers`: disable system dependent high-resolution timers.
- `--enable-chroma`: enable Chroma regression tests. A compiled version of Chroma is assumed to be present. 


More details on the *Getting started* menu entries on the left. 


This document was updated on November 2016. 
{: .notice}

