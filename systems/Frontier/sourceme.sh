
echo spack
. /autofs/nccs-svm1_home1/paboyle/Crusher/Grid/spack/share/spack/setup-env.sh

module load amd/7.0.2
module load cray-fftw
module load craype-accel-amd-gfx90a
mkdir $HOME/LD_PATH
ln -s /opt/rocm-6.4.2/lib/libamdhip* $HOME/LD_PATH

#Ugly hacks to get down level software working on current system
export LD_LIBRARY_PATH=/opt/cray/libfabric/1.20.1/lib64/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/gcc/mpfr/3.1.4/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=`pwd`/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/LD_PATH/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/rocm-7.0.2/lib
