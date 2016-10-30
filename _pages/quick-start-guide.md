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
../configure --enable-precision=double --enable-simd=AVX --enable-comms=mpi-auto --prefix=<path>
```

where `--enable-precision=` set the default precision (`single` or `double`),
`--enable-simd=` set the SIMD type (see possible values below), `--enable-
comms=` set the protocol used for communications (`none`, `mpi`, `mpi-auto` or
`shmem`), and `<path>` should be replaced by the prefix path where you want to
install Grid. The `mpi-auto` communication option set `configure` to determine
automatically how to link to MPI. Other options are available, use `configure
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

