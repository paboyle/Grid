#!/usr/bin/env bash

autoreconf -fvi
./scripts/update_eigen.sh eigen-3.2.9.tar.bz2
./scripts/filelist
