---
title : "API Documentation"
author_profile: false
excerpt: "Grids"
header:
  overlay_color: "#5DADE2"
permalink: /docs/API/grids.html
sidebar:
  nav : docs
---

A `Grid` object defines the geometry of a global cartesian array, and through inheritance
provides access to message passing decomposition, the local lattice, and the message passing primitives.

The constructor requires parameters to indicate how the spatial (and temporal) indices
are decomposed across MPI tasks and SIMD lanes of the vector length.
We use a partial vectorisation transformation, must select
which space-time dimensions participate in SIMD vectorisation.
The Lattice containers are defined to have opaque internal layout, hiding this layout transformation.

We define GridCartesian and GridRedBlackCartesian which both inherit from `GridBase`:

```c++
class GridCartesian        : public GridBase 
class GridRedBlackCartesian: public GridBase 
```

The simplest Cartesian Grid constructor distributes across `MPI_COMM_WORLD`:

```c++
/////////////////////////////////////////////////////////////////////////
// Construct from comm world
/////////////////////////////////////////////////////////////////////////
GridCartesian(const Coordinate &dimensions,
  const Coordinate &simd_layout,
  const Coordinate &processor_grid);
```

A second constructor will create a child communicator from a previously declared Grid.
This allows to subdivide the processor grid, and also to define lattices of differing dimensionalities and sizes,
useful for both Chiral fermions, lower dimensional operations, and multigrid:

```c++
/////////////////////////////////////////////////////////////////////////
// Constructor takes a parent grid and possibly subdivides communicator.
/////////////////////////////////////////////////////////////////////////
GridCartesian(const Coordinate &dimensions,
        const Coordinate &simd_layout,
        const Coordinate &processor_grid,
        const GridCartesian &parent,int &split_rank);
```

The Grid object provides much `internal` functionality to map a lattice site to 
a node and lexicographic index. These are not needed by code interfacing
to the data parallel layer.

**Example** (`tests/solver/Test_split_grid.cc`):

```c++
  const int Ls=8;

  ////////////////////////////////////////////
  // Grids distributed across full machine
  // pick up default command line args
  ////////////////////////////////////////////
  Grid_init(&argc,&argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  GridCartesian * UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
                            GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * rbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  /////////////////////////////////////////////
  // Split into N copies of 1^4 mpi communicators
  /////////////////////////////////////////////
  Coordinate mpi_split (mpi_layout.size(),1);
  GridCartesian * SGrid = new GridCartesian(GridDefaultLatt(),
                                GridDefaultSimd(Nd,vComplex::Nsimd()),
                                mpi_split,
                                *UGrid); 

  GridCartesian         * SFGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,SGrid);
  GridRedBlackCartesian * SrbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(SGrid);
  GridRedBlackCartesian * SFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,SGrid);
```
