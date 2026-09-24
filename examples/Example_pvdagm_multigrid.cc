/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_multigrid.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */

//
// The canonical three-level mrhs PVdagM multigrid driver, on the library
// objects (Grid/algorithms/multigrid/PVdagMMultiGrid.h).
//
// Every parameter is in ONE serialisable struct, read from XML:
//
//   ./Example_pvdagm_multigrid --grid 48.48.48.96 --mpi ... \
//                              --pvdagm-params params.xml
//
// With no --pvdagm-params the built-in defaults run (hot-start gauge field
// unless Config is set in the file).  A missing file gets a template
// written next to it and the program exits: the template documents every
// parameter.  The effective parameters are always printed, so the log
// describes its own run.
//
// Compile-time: NBASIS defaults to 60, the production basis; a subspace file
// may hold more vectors, the load reads only the first NBASIS.  -DNBASIS=8
// cuts it down for laptop runs.
// examples/Makefile.am builds the fp32-coarse and fp32-dense-inversion
// variants below as their own binaries.
//
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/multigrid/PVdagMMultiGrid.h>

using namespace std;
using namespace Grid;

#ifndef NBASIS
#define NBASIS 60
#endif

// Precision of the coarse + coarse-coarse sector is a compile-time
// instantiation: -DCOARSE_SINGLE builds the fp32 coarse space (the dense
// bottom's apply slab is fp32 either way; its inversion is a configure
// option).  Both levels carry the same site type iVector<CoarseScalar,NBASIS>.
#ifdef COARSE_SINGLE
typedef sTComplexF CoarseScalar_t;
#else
typedef sTComplexD CoarseScalar_t;
#endif

struct PVdagMDriverParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(PVdagMDriverParams,
                                  int,         Ls,
                                  RealD,       Mass,
                                  RealD,       M5,
                                  RealD,       MobiusB,
                                  RealD,       MobiusC,
                                  std::string, Config,          // empty: hot start
                                  int,         Nrhs,
                                  int,         SolveSingleRHS,  // also run Nrhs=1 through the same objects
                                  PVdagMMultiGridParams, MultiGrid);
  PVdagMDriverParams()
    : Ls(24), Mass(0.00078), M5(1.8), MobiusB(1.5), MobiusC(0.5),
      Config(""), Nrhs(12), SolveSingleRHS(1) {};
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  PVdagMDriverParams P;
  {
    std::string pfile("");
    if( GridCmdOptionExists(argv,argv+argc,"--pvdagm-params") )
      pfile = GridCmdOptionPayload(argv,argv+argc,"--pvdagm-params");
    if ( pfile.length() ) {
      bool good; { std::ifstream f(pfile); good = f.good(); }
      if ( !good ) {
        if ( GlobalSharedMemory::WorldRank == 0 ) {
          XmlWriter WR(pfile+".templ");
          write(WR, "PVdagMDriver", P);
          std::cout << GridLogMessage << pfile << " does not exist; template written to "
                    << pfile << ".templ" << std::endl;
        }
        Grid_finalize(); return 0;
      }
      XmlReader RD(pfile);
      read(RD, "PVdagMDriver", P);
    }
    CheckValidity(P.MultiGrid);
    std::cout << GridLogMessage << "PVdagMDriver parameters ("
              << (pfile.length() ? pfile : std::string("defaults")) << "):" << std::endl;
    std::cout << P << std::endl;
  }

  Coordinate latt = GridDefaultLatt();
  Coordinate mpi  = GridDefaultMpi();
  Coordinate fsimd= GridDefaultSimd(Nd,vComplex::Nsimd());

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt,fsimd,mpi);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(P.Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(P.Ls,UGrid);

  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers({1,2,3,4});
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers({5,6,7,8});

  LatticeGaugeField Umu(UGrid);
  if ( P.Config.length() ) {
    std::cout << GridLogMessage << "Reading gauge field " << P.Config << std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu,header,P.Config);
  } else {
    std::cout << GridLogMessage << "Hot start gauge field" << std::endl;
    SU<Nc>::HotConfiguration(RNG4,Umu);
  }

  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,P.Mass,P.M5,P.MobiusB,P.MobiusC);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,   P.M5,P.MobiusB,P.MobiusC);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD> PVdagM_t;
  PVdagM_t PVdagM(Ddwf,Dpv);

  //////////////////////////////////////////////////////////////////////
  // Grids -> coarsening -> dense bottom.  Scope order is lifetime order.
  //
  // The fine operator applied during setup is always fp64, on the fp64
  // basis.  Everything the coarsening produces -- the Galerkin matrix
  // elements and the transfer operator's store -- takes the coarse precision
  // (CoarseScalar_t, a compile-time choice).
  //
  // There is exactly ONE fine transfer operator.  Its STORE follows the
  // coarse sector, since the coarse space is what it feeds, while its import
  // and export accept either fine precision: they are already a layout
  // transformation, and a scalar conversion inside one is free.  So the fp64
  // setup and an fp32 V-cycle share the same object and the same basis store.
  // FinePrecision selects only the fine operator and smoother; the outer
  // Krylov and its true-residual check stay fp64 throughout.
  //////////////////////////////////////////////////////////////////////
  typedef PVdagMMultiGridCoarsening<vSpinColourVector,CoarseScalar_t,NBASIS> Coarsening_t;
  std::cout << GridLogMessage << "Coarse sector precision (compiled): "
            << (sizeof(typename GridTypeMapper<CoarseScalar_t>::scalar_type)==sizeof(ComplexF) ? "fp32" : "fp64")
            << ", nbasis " << NBASIS << std::endl;

  MGCoarseGrids CGrids(FGrid, P.MultiGrid.Setup);
  MGFineGridsF  FGridsF(FGrid);

  LatticeGaugeFieldF UmuF(FGridsF.UGridF);
  precisionChange(UmuF,Umu);
  MobiusFermionF DdwfF(UmuF,*FGridsF.FGridF,*FGridsF.FrbGridF,*FGridsF.UGridF,*FGridsF.UrbGridF,P.Mass,P.M5,P.MobiusB,P.MobiusC);
  MobiusFermionF DpvF (UmuF,*FGridsF.FGridF,*FGridsF.FrbGridF,*FGridsF.UGridF,*FGridsF.UrbGridF,1.0,   P.M5,P.MobiusB,P.MobiusC);

  Coarsening_t  Coarsening(CGrids, FGridsF, P.MultiGrid.Setup);

  Coarsening.GetSubspace(RNG5, PVdagM);
  Coarsening.Coarsen(PVdagM);
  Coarsening.BuildDenseBottom();
  Coarsening.CertifyCoarsening(PVdagM);

  //////////////////////////////////////////////////////////////////////
  // Solves: mrhs then (optionally) single RHS through the SAME objects.
  //////////////////////////////////////////////////////////////////////
  auto RunSolve = [&](int nr)
  {
    PVdagMMultiGridSolver<MobiusFermionD,MobiusFermionF,Coarsening_t> Solver(Ddwf,Dpv,DdwfF,DpvF,Coarsening,P.MultiGrid,nr);
    std::vector<LatticeFermionD> src(nr,FGrid), sol(nr,FGrid);
    for(int r=0;r<nr;r++){ gaussian(RNG5,src[r]); sol[r]=Zero(); }
    Solver.Solve(src,sol);
  };

  RunSolve(P.Nrhs);
  if ( P.SolveSingleRHS && P.Nrhs != 1 ) RunSolve(1);

  Grid_finalize();
}
