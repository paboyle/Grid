#!/bin/bash

if [ $1 = "install" ]
then
    dir=`pwd`
    cd $HOME
    git clone -c feature.manyFiles=true https://github.com/spack/spack.git
    source $HOME/spack/share/spack/setup-env.sh

    spack install autoconf
    spack install automake
    spack install c-lime cppflags=-fPIE
    spack install fftw
    spack install llvm
    spack install gmp
    spack install mpfr
    spack install cuda@11.8
    spack install openmpi
    spack install openssl
    spack install hdf5
else
    source $HOME/spack/share/spack/setup-env.sh
fi

spack load autoconf
spack load automake
spack load c-lime
spack load fftw
spack load llvm
spack load gmp
spack load mpfr
spack load cuda@11.8
spack load openmpi
spack load openssl
spack load hdf5

export FFTW=`spack find --paths fftw    | grep ^fftw   | awk '{print $2}' `
export HDF5=`spack find --paths hdf5    | grep ^hdf5   | awk '{print $2}' `
export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
export MPFR=`spack find --paths mpfr    | grep ^mpfr  | awk '{print $2}' `
export GMP=`spack find --paths gmp      | grep ^gmp | awk '{print $2}' `
export NVIDIA=$CUDA_HOME
export NVIDIALIB=$NVIDIA/targets/x86_64-linux/lib/
export LD_LIBRARY_PATH=$NVIDIALIB:$FFTW/lib/:$MPFR/lib:$LD_LIBRARY_PATH
