#!/usr/bin/env bash

EIGEN_URL='http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2'
FFTW_URL=http://www.fftw.org/fftw-3.3.4.tar.gz

echo "-- deploying Eigen source..."
wget ${EIGEN_URL} --no-check-certificate
./scripts/update_eigen.sh `basename ${EIGEN_URL}`
rm `basename ${EIGEN_URL}`

echo "-- copying fftw prototypes..."
wget ${FFTW_URL}
./scripts/update_fftw.sh `basename ${FFTW_URL}`
rm `basename ${FFTW_URL}`

echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi
