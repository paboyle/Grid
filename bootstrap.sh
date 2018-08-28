#!/usr/bin/env bash

EIGEN_URL='http://bitbucket.org/eigen/eigen/get/3.3.5.tar.bz2'

echo "-- deploying Eigen source..."
ARC=`basename ${EIGEN_URL}`
wget ${EIGEN_URL} --no-check-certificate && ./scripts/update_eigen.sh ${ARC} && rm ${ARC}
# patch for non-portable includes in Eigen 3.3.5
# apparently already fixed in Eigen HEAD so it should not be 
# a problem in the future (A.P.)
patch Grid/Eigen/unsupported/CXX11/Tensor scripts/eigen-3.3.5.Tensor.patch

echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi
