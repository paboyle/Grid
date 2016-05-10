#!/usr/bin/env bash

if (( $# != 1 )); then
    echo "usage: `basename $0` <module name>" 1>&2
    exit 1
fi
NAME=$1

if [ -e "Modules/${NAME}.cc" ] || [ -e "Modules/${NAME}.hpp" ]; then
    echo "error: files Modules/${NAME}.* already exists" 1>&2
    exit 1
fi
sed "s/___FILEBASENAME___/${NAME}/g" Module.cc.template > Modules/${NAME}.cc
sed "s/___FILEBASENAME___/${NAME}/g" Module.hpp.template > Modules/${NAME}.hpp
./make_module_list.sh
