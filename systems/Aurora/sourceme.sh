#export ONEAPI_DEVICE_SELECTOR=level_zero:0.0

module use /soft/modulefiles
module load intel_compute_runtime/release/agama-devel-682.22

export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
#export MPIR_CVAR_CH4_OFI_ENABLE_HMEM=1
git config --global http.proxy http://proxy.alcf.anl.gov:3128

export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file"
