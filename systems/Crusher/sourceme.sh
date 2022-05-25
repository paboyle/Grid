module load PrgEnv-gnu
module load rocm/5.1.0
module load cray-mpich/8.1.15
module load gmp
#module load cray-fftw
module load craype-accel-amd-gfx90a
export LD_LIBRARY_PATH=/opt/gcc/mpfr/3.1.4/lib:$LD_LIBRARY_PATH
#Hack for lib
export LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH
