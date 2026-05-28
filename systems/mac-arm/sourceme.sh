source /Users/peterboyle/QCD//Spack/spack//share/spack/setup-env.sh
export FFTW=`spack find --paths fftw    | grep ^fftw   | awk '{print $2}' `
export CLIME=`spack find --paths c-lime | grep ^c-lime | awk '{print $2}' `
export MPFR=`spack find --paths mpfr    | grep ^mpfr  | awk '{print $2}' `
export OPENSSL=`spack find --paths openssl | grep openssl | awk  '{print $2}' `
export GMP=`spack find --paths gmp      | grep ^gmp | awk '{print $2}' `

