#!/usr/bin/env bash
# Environment for building Grid on Apple Silicon Mac with Homebrew dependencies.
# Usage: source systems/mac-arm/sourceme-homebrew.sh

HOMEBREW=/opt/homebrew

export GMP=${HOMEBREW}/opt/gmp
export MPFR=${HOMEBREW}/opt/mpfr
export FFTW=${HOMEBREW}/opt/fftw
export OPENSSL=${HOMEBREW}/opt/openssl@3
export CLIME=/usr/local

export PATH="${HOMEBREW}/opt/open-mpi/bin:${HOMEBREW}/bin:${PATH}"
export LDFLAGS="-L${OPENSSL}/lib"
export CPPFLAGS="-I${OPENSSL}/include"
