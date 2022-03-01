module load emacs
module load PrgEnv-gnu
module load rocm/4.5.0
module load gmp
module load cray-fftw
module load craype-accel-amd-gfx908
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export LD_LIBRARY_PATH=/opt/cray/pe/gcc/mpfr/3.1.4/lib/:$LD_LIBRARY_PATH
