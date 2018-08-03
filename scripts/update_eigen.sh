#!/usr/bin/env bash

if (( $# != 1 )); then
    echo "usage: `basename $0` <archive>" 1>&2
    exit 1
fi
ARC=$1

INITDIR=`pwd`
rm -f lib/Eigen
rm -rf Eigen

##################
#untar
##################
tar -xf ${ARC}
ARCDIR=`tar -tf ${ARC} | head -n1 | sed -e 's@/.*@@'`

###############################
# Link to a deterministic name
###############################

mv ${ARCDIR} Eigen
ln -s ${INITDIR}/Eigen/Eigen ${INITDIR}/lib/Eigen
ln -s ${INITDIR}/Eigen/unsupported/Eigen ${INITDIR}/lib/Eigen/unsupported

# Eigen source headers
cd ${INITDIR}/lib
echo 'eigen_files =\' > ${INITDIR}/lib/Eigen.inc
find -L Eigen -print | sed 's/^/  /;$q;s/$/ \\/' >> ${INITDIR}/lib/Eigen.inc

###################################
# back to home
###################################
cd ${INITDIR}
