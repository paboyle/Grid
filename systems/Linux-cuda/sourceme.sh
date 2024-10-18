. /home/paboyle/spack/share/spack/setup-env.sh
spack load cuda@12.0.0
spack load c-lime
spack load gmp
spack load mpfr
spack load hdf5
spack load fftw
spack load openmpi
export FFTW=`spack find --paths fftw | grep fftw | cut -c 14-`
export HDF5=`spack find --paths hdf5 | grep hdf5 | cut -c 14-`
export CUDA=`spack find --paths cuda@11.8.0 | grep cuda | cut -c 14-`
export CLIME=`spack find --paths c-lime | grep c-lime| cut -c 15-`
export GMP=`spack find --paths gmp | grep gmp | cut -c 12-`
export MPFR=`spack find --paths mpfr | grep mpfr | cut -c 13-`
export NVIDIALIB=$CUDA/targets/x86_64-linux/lib/
export LD_LIBRARY_PATH=$NVIDIALIB:$LD_LIBRARY_PATH:$HDF5/lib:$FFTW/lib:$CLIME/lib/:$MPFR/lib
