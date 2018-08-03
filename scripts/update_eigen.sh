#!/usr/bin/env bash

if (( $# != 1 )); then
    echo "usage: `basename $0` <archive>" 1>&2
    exit 1
fi
ARC=$1

INITDIR=`pwd`

##################
#untar
##################

tar -xf ${ARC}
ARCDIR=`tar -tf ${ARC} | head -n1 | sed -e 's@/.*@@'`
rm -f ${ARC}

###############################
# Link to a deterministic name
###############################

mv ${ARCDIR}  Eigen

# Eigen source headers
cd ${INITDIR}/Eigen

echo 'eigen_files =\' > ${INITDIR}/lib/Eigen.inc
find Eigen -name "*.h" -print | sed 's/^/  /;$q;s/$/ \\/' >> ${INITDIR}/lib/Eigen.inc

cd ${INITDIR}
echo 'eigen_unsupp_files =\' >> ${INITDIR}/lib/Eigen.inc
find  Eigen/unsupported/Eigen -name "*.h" -print | sed 's/^/  /;$q;s/$/ \\/' >> ${INITDIR}/lib/Eigen.inc



###################################
# back to home
###################################
cd ${INITDIR}

#########################################
# Make grid includes happy
#########################################
mkdir ${INITDIR}/lib/Eigen/

ln -s ${INITDIR}/Eigen/Eigen/* ${INITDIR}/lib/Eigen/
ln -s ${INITDIR}/Eigen/unsupported ${INITDIR}/lib/Eigen/
