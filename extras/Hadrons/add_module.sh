#!/usr/bin/env bash

if (( $# != 2 )); then
    echo "usage: `basename $0` <module name> <namespace>" 1>&2
    exit 1
fi
NAME=$1
NS=$2

mkdir -p Modules/${NS}
if [ -e "Modules/${NS}/${NAME}.cc" ] || [ -e "Modules/${NS}/${NAME}.hpp" ]; then
	echo "error: files Modules/${NS}/${NAME}.* already exists" 1>&2
	exit 1
fi
TMPCC=".${NS}.${NAME}.tmp.cc"
TMPHPP=".${NS}.${NAME}.tmp.hpp"
sed "s/___FILEBASENAME___/${NAME}/g" Modules/templates/Module_in_NS.cc.template  > ${TMPCC}
sed "s/___FILEBASENAME___/${NAME}/g" Modules/templates/Module_in_NS.hpp.template > ${TMPHPP}
sed "s/___NAMESPACE___/${NS}/g" ${TMPCC}  > Modules/${NS}/${NAME}.cc
sed "s/___NAMESPACE___/${NS}/g" ${TMPHPP} > Modules/${NS}/${NAME}.hpp
rm -f ${TMPCC} ${TMPHPP}
./make_module_list.sh
