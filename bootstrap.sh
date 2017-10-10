#!/usr/bin/env bash

EIGEN_URL='http://bitbucket.org/eigen/eigen/get/3.3.3.tar.bz2'

echo "-- deploying Eigen source..."
wget ${EIGEN_URL} --no-check-certificate
./scripts/update_eigen.sh `basename ${EIGEN_URL}`
#rm `basename ${EIGEN_URL}`

echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi
