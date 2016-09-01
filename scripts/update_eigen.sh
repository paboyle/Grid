#!/usr/bin/env bash

if (( $# != 1 )); then
    echo "usage: `basename $0` <archive>" 1>&2
    exit 1
fi
ARC=$1

INITDIR=`pwd`
rm -rf lib/Eigen
ARCDIR=`tar -tf ${ARC} | head -n1 | sed -e 's@/.*@@'`
tar -xf ${ARC}
cd ${ARCDIR}
(tar -cf - Eigen --exclude='*.txt' 2>/dev/null) | tar -xf - -C ../lib/
cd ../lib
echo 'eigen_files =\' > Eigen.inc
find Eigen -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> Eigen.inc
cd ${INITDIR}
rm -rf ${ARCDIR}
