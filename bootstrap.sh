#!/usr/bin/env bash

./scripts/update_eigen.sh eigen-3.2.9.tar.bz2
./scripts/filelist
autoreconf -fvi
