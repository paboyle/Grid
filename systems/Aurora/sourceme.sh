
### After Update (Hopefully Working) Setup
module reset
module use /opt/aurora/24.347.0/spack/unified/0.9.2/install/modulefiles/Core
module use /opt/aurora/24.347.0/spack/unified/0.9.2/install/modulefiles/oneapi/2025.0.5
#module load py-numpy
module unload mpich
module unload oneapi
module use /soft/compilers/oneapi/2025.1.0/modulefiles
module load oneapi/public/2025.1.0
module use /home/bertoni/mpich_module/
module load aurora_test_2025.1

### Before Update Setup
#module restore
#module load oneapi/eng-compiler/2024.07.30.002
#module load mpich/opt/4.3.0rc3
#module load cmake
#module unload gcc/12.2.0
#module use /opt/aurora/24.347.0/spack/unified/0.9.2/install/modulefiles/Core
#module load gcc/13.3.0

##########
##module load oneapi/release/2023.12.15.001
##module load mpich/icc-all-debug-pmix-gpu/52.2
##module load mpich-config/mode/deterministic
##module load intel_compute_runtime/release/821.35

source ~/spack/share/spack/setup-env.sh 
spack load c-lime
spack load openssl
spack load gmp
spack load mpfr
export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
export GMP=`spack find --paths gmp | grep ^gmp | awk '{print $2}' `
export MPFR=`spack find --paths mpfr | grep ^mpfr | awk '{print $2}' `

git config --global http.proxy http://proxy.alcf.anl.gov:3128

#spack load openssl@3.3.1%gcc@12.2.0
#spack load unwind
#export UNWIND=`spack find --paths libunwind  | grep ^libunwind  | awk '{print $2}' `
#export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file"
