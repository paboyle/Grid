Interfacing with external software
========================================

Grid provides a number of important modules, such as solvers and
eigensolvers, that are highly optimized for complex vector/SIMD
architectures, such as the Intel Xeon Phi KNL and Skylake processors.
This growing library, with appropriate interfacing, can be accessed
from existing code. Here we describe interfacing issues and provide
examples.

	  
MPI initialization
--------------------

Grid supports threaded MPI sends and receives and, if running with
more than one thread, requires the MPI_THREAD_MULTIPLE mode of message
passing. If the user initializes MPI before starting Grid, the
appropriate initialization call is::

  MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
  assert(MPI_THREAD_MULTIPLE == provided);

Grid Initialization
---------------------

Grid itself is initialized with a call::

  Grid_init(&argc, &argv);

Command line options include::

  --mpi n.n.n.n   : default MPI decomposition
  --threads n     : default number of OMP threads
  --grid n.n.n.n  : default Grid size
  
where `argc` and `argv` are constructed to simulate the command-line
options described above.  At a minimum one usually provides the
`--grid` and `--mpi` parameters.  The former specifies the lattice
dimensions and the latter specifies the grid of processors (MPI
ranks).  If these parameters are not specified with the `Grid_init`
call, they need to be supplied later when creating Grid fields.

The following Grid procedures are useful for verifying that Grid
"default" values are properly initialized.

=============================================================   ===========================================================================================================
  Grid procedure                                                  returns 
=============================================================   ===========================================================================================================
  std::vector<int> GridDefaultLatt();                            lattice size
  std::vector<int> GridDefaultSimd(int Nd,vComplex::Nsimd());    SIMD layout
  std::vector<int> GridDefaultMpi();                             MPI layout
  int Grid::GridThread::GetThreads();                            number of threads
=============================================================   ===========================================================================================================


MPI coordination
----------------

Grid wants to use its own numbering of MPI ranks and its own
assignment of the lattice coordinates with each rank.  Obviously, the
calling program and Grid must agree on these conventions.  One should
use Grid's Cartesian communicator class to discover the processor
assignments. For a four-dimensional processor grid one can define::

  static Grid::CartesianCommunicator *grid_cart = NULL;
  grid_cart = new Grid::CartesianCommunicator(processors);

where `processors` is of type `std::vector<int>`, with values matching
the MPI processor-layout dimensions specified with the `--mpi`
argument in the `Grid_Init` call.  Then each MPI rank can obtain its
processor coordinate using the Cartesian communicator instantiated
above.  For example, in four dimensions::

  std::vector<int> pePos(4);    
  for(int i=0; i<4; i++)
     pePos[i] = grid_cart->_processor_coor[i];

and each MPI process can get its world rank from its processor
coordinates using::

  int peRank = grid_cart->RankFromProcessorCoor(pePos)
	  
Conversely, each MPI process can get its processor coordinates from
its world rank using::

  grid_cart->ProcessorCoorFromRank(peRank, pePos);

If the calling program initialized MPI before initializing Grid, it is
then important for each MPI process in the calling program to reset
its rank number so it agrees with Grid::

   MPI_Comm comm;
   MPI_Comm_split(MPI_COMM_THISJOB,jobid,peRank,&comm);
   MPI_COMM_THISJOB = comm;

where `MPI_COMM_THISJOB` is initially a copy of `MPI_COMM_WORLD` (with
`jobid = 0`), or it is a split communicator with `jobid` equal to the
index number of the subcommunicator.  Once this is done,::

  MPI_Comm_rank(MPI_COMM_THISJOB, &myrank);

returns a rank that agrees with Grid's `peRank`.

QMP coordination
----------------

If the calling program uses the SciDAC QMP message-passing package, a
call to QMP_comm_split() instead can be used to reassign the ranks.
In the example below, `peGrid` gives the processor-grid dimensions,
usually set on the command line with `-qmp-geom`.

**Example**::
  
  int NDIM = QMP_get_allocated_number_of_dimensions();
  Grid::Grid_init(argc,argv);
  FgridBase::grid_initted=true;
  std::vector<int> processors;
  for(int i=0;i<NDIM;i++) processors.push_back(peGrid[i]);
  Grid::CartesianCommunicator grid_cart(processors);
  std::vector<int> pePos(NDIM);
  for(int i=NDIM-1;i>=0;i--)
     pePos[i] = grid_cart._processor_coor[i];
  int peRank = grid_cart->RankFromProcessorCoor(pePos);
  QMP_comm_split(QMP_comm_get_default(),0,peRank,&qmp_comm);
  QMP_comm_set_default(qmp_comm);

  
Mapping fields between Grid and user layouts
---------------------------------------------

In order to map data between calling-program and Grid layouts, it is
important to know how the lattice sites are distributed across the
processor grid.  A lattice site with coordinates `r[mu]` is assigned
to the processor with processor coordinates `pePos[mu]` according to
the rule::

  pePos[mu] = r[mu]/dim[mu]

where `dim[mu]` is the lattice dimension in the `mu` direction.  For
performance reasons, it is important that the external data layout
follow the same rule.  Then data mapping can be done without
requiring costly communication between ranks.  We assume this is the
case here.

When mapping data to and from Grid, one must choose a lattice object
defined on the appropriate grid, whether it be a full lattice (4D
`GridCartesian`), one of the checkerboards (4D
`GridRedBlackCartesian`), a five-dimensional full grid (5D
`GridCartesian`), or a five-dimensional checkerboard (5D
`GridRedBlackCartesian`).  For example, an improved staggered-fermion
color-vector field `cv` on a single checkerboard would be constructed
using

**Example**::

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian       RBGrid(&Grid);

  typename ImprovedStaggeredFermion::FermionField  cv(RBGrid);

The example above assumes that the grid default values were set in the
`Grid_init` call.  If not, they can be set at this point and passed
when `GridCartesian` is instantiated here.  To map data within an MPI
rank, the external code must iterate over the sites belonging to that
rank (full or checkerboard as appropriate).  Note that the site
coordinates are specified relative to the origin of the lattice
subvolume on that rank. To import data into Grid, the external data on
a single site with coordinates `r` is first copied into the
appropriate Grid scalar object `s`.  Then it is copied into the Grid
lattice field `l` with `pokeLocalSite`::

  pokeLocalSite(const sobj &s, Lattice<vobj> &l, Coordinate &r);

To export data from Grid, the reverse operation starts with::

  peekLocalSite(const sobj &s, Lattice<vobj> &l, Coordinate &r);

and then copies the single-site data from `s` into the corresponding
external type.

Here is an example that maps a single site's worth of data in a MILC
color-vector field to a Grid scalar ColourVector object `cVec` and from
there to the lattice colour-vector field `cv`, as defined above.

**Example**::

  indexToCoords(idx,r);
  ColourVector cVec;
  for(int col=0; col<Nc; col++)
      cVec()()(col) = 
          Complex(src[idx].c[col].real, src[idx].c[col].imag);

  pokeLocalSite(cVec, cv, r);

Here the `indexToCoords()` function is a MILC mapping of the MILC site
index `idx` to the 4D lattice coordinate `r`.

Grid provides block- and multiple-rhs conjugate-gradient solvers. For
this purpose it uses a 5D lattice. To map data to and from Grid data
types, the index for the right-hand-side vector becomes the zeroth
coordinate of a five-dimensional vector `r5`.  The remaining
components of `r5` contain the 4D space-time coordinates.  The
`pokeLocalSite/peekLocalSite` operations then accept the coordinate
`r5`, provided the destination/source lattice object is also 5D.  In
the example below data from a single site specified by `idx`,
belonging to a set of `Ls` MILC color-vector fields, are copied into a
Grid 5D fermion field `cv5`.

**Example**::

    GridCartesian * UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt();
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid)  typename ImprovedStaggeredFermion5D::FermionField  cv5(FrbGrid);

    std::vector<int> r(4);
    indexToCoords(idx,r);
    std::vector<int> r5(1,0);
    for( int d = 0; d < 4; d++ ) r5.push_back(r[d]);

    for( int j = 0; j < Ls; j++ ){
      r5[0] = j;
      ColourVector cVec;
      for(int col=0; col<Nc; col++){
	  cVec()()(col) = 
	      Complex(src[j][idx].c[col].real, src[j][idx].c[col].imag);
      }
      pokeLocalSite(cVec, *(out->cv), r5);
    }

