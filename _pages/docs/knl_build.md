---
layout: single
title : "Documentation"
author_profile: false
excerpt: "Building on a Intel Knights Landing"
header:
  overlay_color: "#5DADE2"
permalink: /docs/knl_build/
sidebar:
  nav : docs
---
{% include base_path %}
The information included in this page has been updated on *November 2016* and it is valid for the release version 0.6.0.


The following configuration is recommended for the [Intel Knights Landing](http://ark.intel.com/products/codename/48999/Knights-Landing) platform:

``` text
../configure --enable-precision=double\
             --enable-simd=KNL        \
             --enable-comms=mpi3-auto \
             --with-gmp=<path>        \
             --with-mpfr=<path>       \
             --enable-mkl             \
             CXX=icpc MPICXX=mpiicpc
```

where `<path>` is the UNIX prefix where GMP and MPFR are installed. If you are working on a Cray machine that does not use the `mpiicpc` wrapper, please use:

``` text
../configure --enable-precision=double\
             --enable-simd=KNL        \
             --enable-comms=mpi3      \
             --with-gmp=<path>        \
             --with-mpfr=<path>       \
             --enable-mkl             \
             CXX=CC CC=cc
```


#### Notes
- [GMP](https://gmplib.org/) is the GNU Multiple Precision Library.
- [MPFR](http://www.mpfr.org/) is a C library for multiple-precision floating-point computations with correct rounding.
- Both libaries are necessary for the RHMC support. 




{% include paginator.html %}