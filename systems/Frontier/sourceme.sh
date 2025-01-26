. /autofs/nccs-svm1_home1/paboyle/Crusher/Grid/spack/share/spack/setup-env.sh
spack load c-lime
module load emacs 
module load PrgEnv-gnu
module load rocm/6.0.0
module load cray-mpich
module load gmp
module load cray-fftw
module load craype-accel-amd-gfx90a
export LD_LIBRARY_PATH=/opt/gcc/mpfr/3.1.4/lib:$LD_LIBRARY_PATH
#Hack for lib
#export LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH
