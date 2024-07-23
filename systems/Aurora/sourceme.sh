source ~/spack/share/spack/setup-env.sh 
spack load c-lime

export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
#export LD_LIBRARY_PATH=${TCMALLOC}/lib:$LD_LIBRARY_PATH

export INTELGT_AUTO_ATTACH_DISABLE=1

#export ONEAPI_DEVICE_SELECTOR=level_zero:0.0

# -ftarget-register-alloc-mode=pvc:default 
# -ftarget-register-alloc-mode=pvc:small
# -ftarget-register-alloc-mode=pvc:large
# -ftarget-register-alloc-mode=pvc:auto
#export MPIR_CVAR_CH4_OFI_ENABLE_HMEM=1

export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
git config --global http.proxy http://proxy.alcf.anl.gov:3128

export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file"
