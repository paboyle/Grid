source $HOME/spack/share/spack/setup-env.sh
spack load llvm@17.0.4
export LD_LIBRARY_PATH=/direct/sdcc+u/paboyle/spack/opt/spack/linux-almalinux8-icelake/gcc-8.5.0/llvm-17.0.4-laufdrcip63ivkadmtgoepwmj3dtztdu/lib:$LD_LIBRARY_PATH
module load openmpi/4.1.8
spack load c-lime
