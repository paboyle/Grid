
echo spack
. /autofs/nccs-svm1_home1/paboyle/spack/share/spack/setup-env.sh

export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
export MPFR=`spack find --paths mpfr    | grep ^mpfr  | awk '{print $2}' `
export OPENSSL=`spack find --paths openssl | grep openssl | awk  '{print $2}' `
export GMP=`spack find --paths gmp      | grep ^gmp | awk '{print $2}' `

module load cce/21.0.0
module load cpe/26.03
module load rocm/7.0.2
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/rocm-7.0.2/lib/llvm/lib/:$LD_LIBRARY_PATH
