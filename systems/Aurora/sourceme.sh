export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
git config --global http.proxy http://proxy.alcf.anl.gov:3128

source ~/spack/share/spack/setup-env.sh 
spack load c-lime
spack load openssl@3.3.1%gcc@12.2.0
spack load unwind
export UNWIND=`spack find --paths libunwind  | grep ^libunwind  | awk '{print $2}' `
export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file"
