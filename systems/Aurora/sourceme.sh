#module load oneapi/release/2023.12.15.001
#module load mpich/icc-all-debug-pmix-gpu/52.2
#module load mpich-config/mode/deterministic
#module load intel_compute_runtime/release/821.35

#module load pti-gpu #used to be uncommented

source ~paboyle/GPT/Grid/build/sourceme.sh #added
source ~/spack/share/spack/setup-env.sh 
spack load c-lime
spack load openssl
export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
git config --global http.proxy http://proxy.alcf.anl.gov:3128
export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file"
