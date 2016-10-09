#!/usr/bin/env bash

if (( $# != 1 )); then
    echo "usage: `basename $0` <archive>" 1>&2
    exit 1
fi
ARC=$1

INITDIR=`pwd`
rm -rf lib/fftw
mkdir lib/fftw

ARCDIR=`tar -tf ${ARC} | head -n1 | sed -e 's@/.*@@'`
tar -xf ${ARC}
cp ${ARCDIR}/api/fftw3.h lib/fftw/

cd ${INITDIR}
rm -rf ${ARCDIR}
