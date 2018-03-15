---
title : "API Documentation"
author_profile: false
excerpt: "Random number generators"
header:
  overlay_color: "#5DADE2"
permalink: /docs/API/random_number_generators.html
sidebar:
  nav : docs
---

Grid provides three configure time options for random the number generator engine.

* `sitmo`
* `ranlux48`
* `mt19937`

The selection is controlled by the `--enable-rng=<option>` flag.

Sitmo is the default Grid RNG and is recommended. It is a hash based RNG that is cryptographically secure and has 

#. passed the BigCrush tests

#. can Skip forward an arbitrary distance (up to 2^256) in O(1) time

We use Skip functionality to place each site in an independent well separated stream.
The Skip was trivially parallelised, important in a many core node,
and gives very low overhead parallel RNG initialisation.

Our implementation of parallel RNG

* Has passed the BigCrush tests **drawing once from each site RNG** in a round robin fashion.

This test is applied in `tests/testu01/Test_smallcrush.cc`

The interface is as follows::

```c++
  class GridSerialRNG { 
    GridSerialRNG();
    void SeedFixedIntegers(const std::vector<int> &seeds);
  }

  class GridParallelRNG {
    GridParallelRNG(GridBase *grid);
    void SeedFixedIntegers(const std::vector<int> &seeds);
  }

  template <class vobj> void random(GridParallelRNG &rng,Lattice<vobj> &l)   { rng.fill(l,rng._uniform);  }
  template <class vobj> void gaussian(GridParallelRNG &rng,Lattice<vobj> &l) { rng.fill(l,rng._gaussian); }

  template <class sobj> void random(GridSerialRNG &rng,sobj &l)   { rng.fill(l,rng._uniform  ); }
  template <class sobj> void gaussian(GridSerialRNG &rng,sobj &l) { rng.fill(l,rng._gaussian ); }
```

* Serial RNG's are used to assign scalar fields. 

* Parallel RNG's are used to assign lattice fields and must subdivide the field grid (need not be conformable).

It is the API users responsibility to initialise, manage, save and restore these RNG state for their algorithm.
In particular there is no single globally managed RNG state. 

Input/Output routines are provided for saving and restoring RNG states.

`lib/parallelIO/BinaryIO.h`:

```c++

  ////////////////////////////////////////////////////////////////////////////
  // Read a RNG;  use IOobject and lexico map to an array of state 
  ////////////////////////////////////////////////////////////////////////////
  static void readRNG(GridSerialRNG &serial,
			     GridParallelRNG &parallel,
			     std::string file,
			     Integer offset,
			     uint32_t &nersc_csum,
			     uint32_t &scidac_csuma,
			     uint32_t &scidac_csumb)
  ////////////////////////////////////////////////////////////////////////////
  // Write a RNG; lexico map to an array of state and use IOobject
  ////////////////////////////////////////////////////////////////////////////
  static void writeRNG(GridSerialRNG &serial,
			      GridParallelRNG &parallel,
			      std::string file,
			      Integer offset,
			      uint32_t &nersc_csum,
			      uint32_t &scidac_csuma,
			      uint32_t &scidac_csumb)

lib/parallelIO/NerscIO.h::

  void writeRNGState(GridSerialRNG &serial,GridParallelRNG &parallel,std::string file);

  void readRNG(GridSerialRNG &serial,
	       GridParallelRNG &parallel,
	       std::string file,
   	       Integer offset,
 	       uint32_t &nersc_csum,
	       uint32_t &scidac_csuma,
	       uint32_t &scidac_csumb);
```

**Example**:

```c++
  NerscIO::writeRNGState(sRNG,pRNG,rfile);
```
