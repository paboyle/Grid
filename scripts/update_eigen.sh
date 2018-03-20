#!/usr/bin/env bash

if (( $# != 1 )); then
    echo "usage: `basename $0` <archive>" 1>&2
    exit 1
fi
ARC=$1

INITDIR=`pwd`
cd ../lib
echo 'eigen_files =\' > Eigen.inc
find Eigen -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> Eigen.inc
cd ${INITDIR}

