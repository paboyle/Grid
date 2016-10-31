---
layout: single
title : "Quick start guide"
author_profile: false
excerpt: "How to install"
header:
  overlay_color: "#333"
permalink: /docs/quick-start-guide/
---
### Installation
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
mkdir build; cd build
../configure --enable-precision=double --enable-simd=AVX --enable-comms=mpi-auto --prefix=<path>
```

The list of possible options from the `./configure --help` command is 

``` text
Optional Features:
  --disable-option-checking  ignore unrecognized --enable/--with options
  --disable-FEATURE       do not include FEATURE (same as --enable-FEATURE=no)
  --enable-FEATURE[=ARG]  include FEATURE [ARG=yes]
  --enable-silent-rules   less verbose build output (undo: "make V=1")
  --disable-silent-rules  verbose build output (undo: "make V=0")
  --enable-dependency-tracking
                          do not reject slow dependency extractors
  --disable-dependency-tracking
                          speeds up one-time build
  --disable-openmp        do not use OpenMP
  --enable-lapack=yes|no|prefix
                          enable LAPACK
  --enable-numa=yes|no|prefix
                          enable first touch numa opt
  --enable-simd=SSE4|AVX|AVXFMA4|AVXFMA|AVX2|AVX512|AVX512MIC|IMCI|KNL|KNC
                          Select instructions to be SSE4.0, AVX 1.0, AVX
                          2.0+FMA, AVX 512, IMCI
  --enable-precision=single|double
                          Select default word size of Real
  --enable-comms=none|mpi|mpi-auto|shmem
                          Select communications
  --enable-rng=ranlux48|mt19937
                          Select Random Number Generator to be used
  --enable-timers         Enable system dependent high res timers
  --enable-chroma         Expect chroma compiled under c++11
  --enable-doxygen        enable documentation generation with doxygen (auto)
  --enable-dot            use 'dot' to generate graphs in doxygen (auto)
  --enable-html-docs      enable HTML generation with doxygen (yes)
  --enable-latex-docs     enable LaTeX documentation generation with doxygen
                          (no)
```


and `<path>` should be replaced by the prefix path where you want to
install Grid. The `mpi-auto` communication option let `configure` to determine
automatically how to link to MPI. 
Like with any other program using GNU autotool, the `CXX`, `CXXFLAGS`, `LDFLAGS`, ... environment variables can be modified to
customise the build.

Finally, you can build and install Grid:

``` bash
make
make install
```

To minimise the build time, only the tests at the root of the `tests` directory are built by default. If you want to build tests in the sub-directory `<subdir>` you can execute:

``` bash
make -C tests/<subdir> tests
```

### Possible SIMD types

The following options can be use with the `--enable-simd=` option to target different SIMD instruction sets:

| String      | Description                            |
| ----------- | -------------------------------------- |
| `GEN`       | generic portable vector code           |
| `SSE4`      | SSE 4.2 (128 bit)                      |
| `AVX`       | AVX (256 bit)                          |
| `AVXFMA4`   | AVX (256 bit) + FMA                    |
| `AVX2`      | AVX 2 (256 bit)                        |
| `AVX512`    | AVX 512 bit                            |
| `AVX512MIC` | AVX 512 bit for Intel MIC architecture |
| `ICMI`      | Intel ICMI instructions (512 bit)      |

Alternatively, some CPU codenames can be directly used:

| String      | Description                            |
| ----------- | -------------------------------------- |
| `KNC`       | [Intel Knights Corner](http://ark.intel.com/products/codename/57721/Knights-Corner) |
| `KNL`       | [Intel Knights Landing](http://ark.intel.com/products/codename/48999/Knights-Landing) |

