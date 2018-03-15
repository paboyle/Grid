---
##layout : single
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

{% if site.option=="web" %}
Use the sidebar on the left to navigate. 
{% endif %}

_May 2017 : The API description and Lattice Theories sections in the sidebar are work in progress_. 

### Version history

* May 2017 [version 0.7.0](https://github.com/paboyle/Grid/tree/release/v0.7.0)
* November 2016 [version 0.6.0](https://github.com/paboyle/Grid/tree/release/v0.6.0)


Grid is primarily an application development interface (API) for structured Cartesian grid codes and written in C++11.
In particular it is aimed at Lattice Field Theory simulations in general gauge theories, but with a particular emphasis
on supporting SU(3) and U(1) gauge theories relevant to hadronic physics.


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


### Who will use this library

As an application development interface *Grid* is primarily a programmers tool providing the
building blocks and primitives for constructing lattice gauge theory programmes. 

Grid functionality includes:

* Data parallel primitives, similar to QDP++
* gauge and fermion actions 
* solvers
* gauge and fermion force terms
* integrators and (R)HMC.
* parallel field I/O 
* object serialisation (text, XML, JSON...)

Grid is intended to enable the rapid and easy development of code with reasonably competitive performance.

It is first and foremost a *library* to which people can programme, and develop new algorithms and measurements.
As such, it is very much hoped that peoples principle point of contact with Grid will be in
the wonderfully rich C++ language. Since import and export procedures are provided for the opaque lattice types
it should be possible to call Grid from other code bases. 

Grid is most tightly coupled to the Hadrons package 
developed principally by Antonin Portelli. 
This package is entirely composed against the Grid data parallel interface.

Interfacing to other packages is also possible.

Several regression tests that combine Grid with Chroma are included in the Grid distribution.
Further, Grid has been successfully interfaced to 

* The Columbia Physics System
* The MILC code

### Data parallel interface

Most users will wish to interact with Grid above the data parallel *Lattice* interface. At this level
a programme is simply written as a series of statements, addressing entire lattice objects. 


Implementation details may be provided to explain how the code works, but are not strictly part of the API.

**Example** 

   For example, as an implementation detail, in a single programme multiple data (SPMD) message passing supercomputer the main programme is trivially replicated on each computing node. The data parallel operations are called *collectively* by all nodes. Any scalar values returned by the various reduction routines are the same on each node, resulting in (for example) the same decision being made by all nodes to terminate an iterative solver on the same iteration. 



### Internal development

Internal developers may contribute to Grid at a level below the data parallel interface.

Specifically, development of new lattice Dirac operators, for example, 
or any codes directly interacting with the 

* Communicators

* Simd 

* Tensor

* Stencil 

will make use of facilities provided by to assist the creation of high performance code. The internal data layout complexities
will be exposed to some degree and the interfaces are subject to change without notice as HPC architectures change.

Since some of the internal implementation details are needed to explain the design strategy of grid these will be 
documented, but labelled as *implementation dependent*

Reasonable endeavours will be made to preserve functionality where practical but no guarantees are made.



{% include paginator.html %}
