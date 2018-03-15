---
title : "API Documentation"
author_profile: false
excerpt: "Data parallel API"
header:
  overlay_color: "#5DADE2"
permalink: /docs/API/introduction.html
sidebar:
  nav : docs
---


Data parallel array indices are divided into two types. 

* Internal indices, such as complex, colour, spin degrees of freedom 
* spatial (space-time) indices.

The ranges of all internal degrees are determined by template parameters, 
and known at compile time. The ranges of spatial indices are dynamic, run time
values and the Cartesian structure information is contained and accessed via `Grid` objects.

Grid objects are the controlling entity for the decomposition of a distributed `Lattice`
array across MPI tasks, nodes, SIMD lanes, accelerators. Threaded loops are used
as appropriate on host code.

Data parallel operations can only be performed between Lattice objects constructed
from the same Grid pointer. These are called `conformable` operations.

We will focus initially on the internal indices as these are the building blocks assembled
in Lattice container classes. Every Lattice container class constructor requires a Grid object 
pointer. 

