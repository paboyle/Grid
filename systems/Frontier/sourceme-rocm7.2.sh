module unload slate 2>/dev/null

echo spack
. /autofs/nccs-svm1_home1/paboyle/spack-frontier/share/spack/setup-env.sh

export SLATE=`spack find --paths slate   | grep ^slate  | awk '{print $2}' `
export BLASPP=`spack find --paths blaspp | grep ^blaspp | awk '{print $2}' `
export LAPACKPP=`spack find --paths lapackpp | grep ^lapackpp | awk '{print $2}' `

echo SLATE    $SLATE
echo LAPACKPP $LAPACKPP
echo BLASPP   $BLASPP

ls $SLATE/lib64/libslate.so $LAPACKPP/lib64/liblapackpp.so 2>/dev/null || ls $SLATE/lib $LAPACKPP/lib
ldd $SLATE/lib64/libslate.so | grep -E "blaspp|lapackpp|rocblas|amdhip|mpi_cray|libsci"
ldd $LAPACKPP/lib64/liblapackpp.so | grep -E "blaspp|libsci"

export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
export MPFR=`spack find --paths mpfr    | grep ^mpfr  | awk '{print $2}' `
export OPENSSL=`spack find --paths openssl | grep openssl | awk  '{print $2}' `
export GMP=`spack find --paths gmp      | grep ^gmp | awk '{print $2}' `

module load cce/21.0.0
module load cpe/26.03
module load rocm/7.2.0
export LD_LIBRARY_PATH=/opt/rocm-7.2.0/lib/llvm/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LAPACKPP/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$BLASPP/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SLATE/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_LD_LIBRARY_PATH
module load emacs
