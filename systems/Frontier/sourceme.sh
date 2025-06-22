
echo spack
. /autofs/nccs-svm1_home1/paboyle/Crusher/Grid/spack/share/spack/setup-env.sh

module load cce/15.0.1
module load rocm/5.3.0
module load cray-fftw
module load craype-accel-amd-gfx90a

#Ugly hacks to get down level software working on current system
export LD_LIBRARY_PATH=/opt/cray/libfabric/1.20.1/lib64/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/gcc/mpfr/3.1.4/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=`pwd`/:$LD_LIBRARY_PATH
ln -s /opt/rocm-6.0.0/lib/libamdhip64.so.6 .

#echo spack load c-lime
#spack load c-lime
#module load emacs 
##module load PrgEnv-gnu
##module load cray-mpich
##module load cray-fftw
##module load craype-accel-amd-gfx90a
##export LD_LIBRARY_PATH=/opt/gcc/mpfr/3.1.4/lib:$LD_LIBRARY_PATH
#Hack for lib
##export LD_LIBRARY_PATH=`pwd`/:$LD_LIBRARY_PATH
