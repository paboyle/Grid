---
layout: single
title : "Documentation"
author_profile: false
excerpt: "Welcome to the Grid documentation pages"
header:
  overlay_color: "#5DADE2"
  #cta_label: "Download documentation"
  #cta_url: "https://www.google.com"
sidebar:
  nav : docs
permalink: /docs/
---

{% include base_path %}

We are currently working on the full documentation.

Use the sidebar on the left to navigate. 

_May 2017 : The API description and Lattice Theories sections in the sidebar are work in progress_. 

### Version history

* May 2017 [version 0.7.0](https://github.com/paboyle/Grid/tree/release/v0.7.0)
* November 2016 [version 0.6.0](https://github.com/paboyle/Grid/tree/release/v0.6.0)

### Description 

This library provides data parallel C++ container classes with internal memory layout
that is transformed to map efficiently to SIMD architectures. CSHIFT facilities
are provided, similar to HPF and cmfortran, and user control is given over the mapping of
array indices to both MPI tasks and SIMD processing elements.

* Identically shaped arrays then be processed with perfect data parallelisation.
* Such identically shaped arrays are called conformable arrays.

The transformation is based on the observation that Cartesian array processing involves
identical processing to be performed on different regions of the Cartesian array.

The library will both geometrically decompose into MPI tasks and across SIMD lanes.
Local vector loops are parallelised with OpenMP pragmas.

Data parallel array operations can then be specified with a SINGLE data parallel paradigm, but
optimally use MPI, OpenMP and SIMD parallelism under the hood. This is a significant simplification
for most programmers.

The layout transformations are parametrised by the SIMD vector length. This adapts according to the architecture.
Presently SSE4 (128 bit) AVX, AVX2, QPX (256 bit), IMCI, and AVX512 (512 bit) targets are supported (ARM NEON on the way).

These are presented as `vRealF`, `vRealD`, `vComplexF`, and `vComplexD` internal vector data types. These may be useful in themselves for other programmers.
The corresponding scalar types are named `RealF`, `RealD`, `ComplexF` and `ComplexD`.

MPI, OpenMP, and SIMD parallelism are present in the library.
Please see [this paper](https://arxiv.org/abs/1512.03487) for more detail.




{% include paginator.html %}
