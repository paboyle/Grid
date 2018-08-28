#!/usr/bin/env bash

if (( $# != 1 )); then
    echo "usage: `basename $0` <archive>" 1>&2
    exit 1
fi
ARC=$1

INITDIR=`pwd`
rm -f Grid/Eigen
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
ln -s ${INITDIR}/Eigen/Eigen ${INITDIR}/Grid/Eigen
ln -s ${INITDIR}/Eigen/unsupported/Eigen ${INITDIR}/Grid/Eigen/unsupported

# Eigen source headers
cd ${INITDIR}/Grid
echo 'eigen_files =\' > ${INITDIR}/Grid/Eigen.inc
find -L Eigen -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> ${INITDIR}/Grid/Eigen.inc

###################################
# back to home
###################################
cd ${INITDIR}
