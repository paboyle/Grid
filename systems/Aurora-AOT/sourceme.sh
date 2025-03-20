#module load oneapi/release/2023.12.15.001
#module load mpich/icc-all-debug-pmix-gpu/52.2
#module load mpich-config/mode/deterministic
#module load intel_compute_runtime/release/821.35

source ~/spack/share/spack/setup-env.sh 
spack load c-lime
spack load openssl
export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
git config --global http.proxy http://proxy.alcf.anl.gov:3128
export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file"
