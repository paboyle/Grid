#!/usr/bin/env bash

EIGEN_URL='http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2'

echo "-- deploying Eigen source..."
cd lib
rm -rf Eigen
git clone https://github.com/eigenteam/eigen-git-mirror.git
mv eigen-git-mirror/Eigen .
echo 'eigen_files =\' > Eigen.inc
find Eigen -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> Eigen.inc
cd ..

echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi
