---
title : "Documentation"
author_profile: false
excerpt: "Execution model"
header:
  overlay_color: "#5DADE2"
permalink: /docs/execution-model/
sidebar:
  nav : docs
---

Grid is intended to support performance portability across a many of platforms ranging from single processors
to message passing CPU clusters and accelerated computing nodes.

The library provides data parallel C++ container classes with internal memory layout that is transformed to map efficiently to SIMD architectures. CSHIFT facilities are provided, similar to HPF and cmfortran, and user control is given over the mapping of array indices to both MPI tasks and SIMD processing elements.

Identically shaped arrays then be processed with perfect data parallelisation.
Such identically shaped arrays are called conformable arrays.
The transformation is based on the observation that Cartesian array processing involves identical processing to be performed on different regions of the Cartesian array.

The library will both geometrically decompose into MPI tasks and across SIMD lanes. Local vector loops are parallelised with OpenMP pragmas.

Data parallel array operations can then be specified with a SINGLE data parallel paradigm, but optimally use MPI, OpenMP and SIMD parallelism under the hood. This is a significant simplification for most programmers.

The two broad optimisation targets are:

* MPI, OpenMP, and SIMD parallelism 

Presently SSE4, ARM NEON (128 bits) AVX, AVX2, QPX (256 bits), and AVX512 (512 bits) targets are supported
with aggressive use of architecture vectorisation intrinsic functions.

* MPI between nodes with and data parallel offload to GPU's.

For the latter generic C++ code is used both on the host and on the GPU, with a common vectorisation
granularity.

### Accelerator memory model

For accelerator targets it is assumed that heap allocations can be shared between the CPU
and the accelerator. This corresponds to lattice fields having their memory allocated with 
*cudaMallocManaged* with Nvidia GPU's. 

Grid does not assume that stack or data segments share a common address space with an accelerator.

* This constraint presently rules out porting Grid to AMD GPU's which do not support managed memory.

* At some point in the future a cacheing strategy may be implemented to enable running on AMD GPU's

This document was updated on March 2018. 
{: .notice}

