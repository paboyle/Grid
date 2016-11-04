---
layout: single
title : "Documentation"
author_profile: false
excerpt: "Supported SIMD architectures"
header:
  overlay_color: "#5DADE2"
permalink: /docs/simd_targets/
sidebar:
  nav : docs
---
{% include base_path %}

The following options can be used for `--enable-simd=` flag to target different SIMD instruction sets:

| `<code>`    | Description                            |
| ----------- | -------------------------------------- |
| `GEN`       | generic portable vector code           |
| `SSE4`      | SSE 4.2 (128 bit)                      |
| `AVX`       | AVX (256 bit)                          |
| `AVXFMA`    | AVX (256 bit) + FMA                    |
| `AVXFMA4`   | AVX (256 bit) + FMA4                   |
| `AVX2`      | AVX 2 (256 bit)                        |
| `AVX512`    | AVX 512 bit                            |
| `QPX`       | QPX (256 bit)                          |

Alternatively, some CPU codenames can be directly used:

| `<code>`    | Description                            |
| ----------- | -------------------------------------- |
| `KNC`       | [Intel Xeon Phi codename Knights Corner](http://ark.intel.com/products/codename/57721/Knights-Corner) |
| `KNL`       | [Intel Xeon Phi codename Knights Landing](http://ark.intel.com/products/codename/48999/Knights-Landing) |
| `BGQ`       | Blue Gene/Q                            |

#### Notes (November 2016):
- We are currently supporting AVX512 only for the Intel compiler. GCC and clang will appear in future versions, if the compilers provide sufficient support.
- For BG/Q only [bgclang](http://trac.alcf.anl.gov/projects/llvm-bgq) is supported. We do not presently plan to support more compilers for this platform.
- BG/Q performances are currently rather poor. This is being investigated for future versions.


{% include paginator.html %}