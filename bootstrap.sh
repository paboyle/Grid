#!/usr/bin/env bash

EIGEN_URL='http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2'

echo "-- deploying Eigen source..."
cd lib
git clone https://github.com/eigenteam/eigen-git-mirror.git
cd ..
#wget ${EIGEN_URL} --no-check-certificate && ./scripts/update_eigen.sh `basename ${EIGEN_URL}` && rm `basename ${EIGEN_URL}`

echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi
