export CRAY_ACCEL_TARGET=nvidia80
source /global/homes/p/pboyle/spack/share/spack/setup-env.sh
export MPFR=`spack find --paths mpfr | grep mpfr | cut -c 13-`
export GMP=`spack find --paths gmp   | grep gmp  | cut -c 12-`

module load PrgEnv-gnu cpe-cuda cudatoolkit/12.0
