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

Please send all pull requests to the `develop` branch.

License: GPL v2.

Last update 2016/08/03.

### Description
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
Presently SSE4 (128 bit) AVX, AVX2 (256 bit) and IMCI and AVX512 (512 bit) targets are supported (ARM NEON and BG/Q QPX on the way).

These are presented as `vRealF`, `vRealD`, `vComplexF`, and `vComplexD` internal vector data types. These may be useful in themselves for other programmers.
The corresponding scalar types are named `RealF`, `RealD`, `ComplexF` and `ComplexD`.

MPI, OpenMP, and SIMD parallelism are present in the library.
Please see https://arxiv.org/abs/1512.03487 for more detail.

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

Now you can execute the `configure` script to generate makefiles (here from a build directory):

``` bash
mkdir build; cd build
../configure --enable-precision=double --enable-simd=AVX --enable-comms=mpi --prefix=<path>
```

where `--enable-precision=` set the default precision (`single` or `double`), `--enable-simd=` set the SIMD type (see possible values below), `--enable-comms=` set the protocol used for communications (`none`, `mpi` or `shmem`), and `<path>` should be replaced by the prefix path where you want to install Grid. Other options are available, use `configure --help` to display them. Like with any other program using GNU autotool, the `CXX`, `CXXFLAGS`, `LDFLAGS`, ... environment variables can be modified to customise the build.

Finally, you can build and install Grid:

``` bash
make; make install
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