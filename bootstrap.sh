#!/usr/bin/env bash
set -e

EIGEN_URL='http://bitbucket.org/eigen/eigen/get/3.3.5.tar.bz2'

verify_helper() {
   local cmd="$1" opt="$2" file="$3"
   command -v "$cmd" || return 0  #no error if command does not exist
   "$cmd" "$opt" "$file" && return 0
   echo "ERROR: CHECKSUM VERIFICATION FAILED!!!" >&2
   exit 1
}

verify() {
   local file="scripts/checksum-$1"
   verify_helper sha256sum -c "$file".sha256sum && return 0
   verify_helper md5sum -c "$file".md5sum && return 0
   echo "WARNING: could not verify checksum" >&2
   echo "please install sha256sum or md5sum to increase security" >&2
}

echo "-- deploying Eigen source..."
ARC=`basename ${EIGEN_URL}`
wget "${EIGEN_URL}" --no-check-certificate \
   && verify "${ARC}" \
   && ./scripts/update_eigen.sh "${ARC}" \
   && rm "${ARC}"
# patch for non-portable includes in Eigen 3.3.5
# apparently already fixed in Eigen HEAD so it should not be 
# a problem in the future (A.P.)
patch Eigen/unsupported/Eigen/CXX11/Tensor scripts/eigen-3.3.5.Tensor.patch

echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi
echo "---bootstrap successfully completed---"
