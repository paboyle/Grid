export https_proxy=http://proxy-chain.intel.com:911
module load intel-release
module load intel/mpich
export MPIR_CVAR_CH4_OFI_ENABLE_GPU_PIPELINE=1
export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file"
