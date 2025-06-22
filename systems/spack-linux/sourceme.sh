source $HOME/spack/share/spack/setup-env.sh
spack load llvm@12
spack load autoconf%clang@12.0.1
spack load automake%clang@12.0.1
spack load c-lime%clang@12.0.1
spack load fftw%clang@12.0.1
spack load gmp%clang@12.0.1
spack load mpfr%clang@12.0.1
spack load openmpi%clang@12.0.1
spack load openssl%clang@12.0.1
spack load hdf5+cxx%clang@12.0.1
spack load cmake%clang@12.0.1
export FFTW=`spack find --paths fftw%clang@12.0.1    | grep ^fftw   | awk '{print $2}' `
export HDF5=`spack find --paths hdf5+cxx%clang@12.0.1   | grep ^hdf5   | awk '{print $2}' `
export CLIME=`spack find --paths c-lime%clang@12.0.1 | grep ^c-lime | awk '{print $2}' `
export MPFR=`spack find --paths mpfr%clang@12.0.1    | grep ^mpfr  | awk '{print $2}' `
export LLVM=`spack find --paths llvm@12    | grep ^llvm  | awk '{print $2}' `
export OPENSSL=`spack find --paths openssl%clang@12.0.1 | grep openssl | awk  '{print $2}' `
export GMP=`spack find --paths gmp%clang@12.0.1      | grep ^gmp | awk '{print $2}' `
export TCLAP=`spack find --paths tclap%clang@12.0.1  | grep ^tclap | awk '{print $2}' `
export LD_LIBRARY_PATH=${TCLAP}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$MPFR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$GMP/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$FFTW/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LLVM/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LLVM/lib/x86_64-unknown-linux-gnu/:$LD_LIBRARY_PATH

ulimit -s 81920
