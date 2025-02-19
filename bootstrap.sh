#!/usr/bin/env bash
set -e

EIGEN_URL='https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2'
EIGEN_SHA256SUM='b4c198460eba6f28d34894e3a5710998818515104d6e74e5cc331ce31e46e626'


echo "-- deploying Eigen source..."
ARC=$(basename ${EIGEN_URL})
wget ${EIGEN_URL} --no-check-certificate
if command -v sha256sum; then
   echo "$EIGEN_SHA256SUM  $(basename "$EIGEN_URL")" \
      | sha256sum --check || exit 1
else
   echo "WARNING: could not verify checksum, please install sha256sum" >&2
fi
./scripts/update_eigen.sh "${ARC}"
rm "${ARC}"
echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi
