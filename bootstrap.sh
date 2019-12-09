#!/usr/bin/env bash
set -e

EIGEN_URL='https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2'
EIGEN_SHA256SUM='685adf14bd8e9c015b78097c1dc22f2f01343756f196acdc76a678e1ae352e11'


echo "-- deploying Eigen source..."
ARC=`basename ${EIGEN_URL}`
wget ${EIGEN_URL} --no-check-certificate
if command -v sha256sum; then
   echo "$EIGEN_SHA256SUM  $(basename "$EIGEN_URL")" \
      | sha256sum --check || exit 1
else
   echo "WARNING: could not verify checksum, please install sha256sum" >&2
fi
./scripts/update_eigen.sh ${ARC}
rm ${ARC}
# patch for non-portable includes in Eigen 3.3.5
# apparently already fixed in Eigen HEAD so it should not be 
# a problem in the future (A.P.)
patch Eigen/unsupported/Eigen/CXX11/Tensor scripts/eigen-3.3.5.Tensor.patch

echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi
