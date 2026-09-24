/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_hdcg_multigrid.cc

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
// The two-level mrhs HDCG driver on the library objects
// (Grid/algorithms/multigrid/HDCGMultiGrid.h): the Schur-preconditioned
// M^dag M of a Mobius/Shamir fermion on one checkerboard, coarsened on V2.
//
//   ./Example_hdcg_multigrid --grid 48.48.48.96 --mpi ... --hdcg-params params.xml
//
// With no --hdcg-params the built-in defaults run (hot-start gauge field
// unless Config is set).  A missing file gets a template written next to
// it and the program exits: the template documents every parameter.
//
// Compile-time: -DNBASIS=8 cuts the basis down for laptop runs;
// -DCOARSE_SINGLE builds the fp32 coarse space, which examples/Makefile.am
// builds as its own binary, Example_hdcg_multigrid_fp32coarse.
//
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/multigrid/HDCGMultiGrid.h>

using namespace std;
using namespace Grid;

#ifndef NBASIS
#define NBASIS 62
#endif

#ifdef COARSE_SINGLE
typedef sTComplexF CoarseScalar_t;
#else
typedef sTComplexD CoarseScalar_t;
#endif

struct HDCGDriverParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HDCGDriverParams,
                                  int,         Ls,
                                  RealD,       Mass,
                                  RealD,       M5,
                                  RealD,       MobiusB,
                                  RealD,       MobiusC,
                                  std::string, Config,          // empty: hot start
                                  int,         Checkerboard,    // 0 Even, 1 Odd: the Schur operator's checkerboard
                                  int,         Nrhs,
                                  int,         SolveSingleRHS,  // also run Nrhs=1 through the same objects
                                  HDCGMultiGridParams, MultiGrid);
  HDCGDriverParams()
    : Ls(24), Mass(0.00078), M5(1.8), MobiusB(1.5), MobiusC(0.5),
      Config(""), Checkerboard(0), Nrhs(12), SolveSingleRHS(1) {};
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  HDCGDriverParams P;
  {
    std::string pfile("");
    if( GridCmdOptionExists(argv,argv+argc,"--hdcg-params") )
      pfile = GridCmdOptionPayload(argv,argv+argc,"--hdcg-params");
    if ( pfile.length() ) {
      bool good; { std::ifstream f(pfile); good = f.good(); }
      if ( !good ) {
        if ( GlobalSharedMemory::WorldRank == 0 ) {
          XmlWriter WR(pfile+".templ");
          write(WR, "HDCGDriver", P);
          std::cout << GridLogMessage << pfile << " does not exist; template written to "
                    << pfile << ".templ" << std::endl;
        }
        Grid_finalize(); return 0;
      }
      XmlReader RD(pfile);
      read(RD, "HDCGDriver", P);
    }
    CheckValidity(P.MultiGrid);
    std::cout << GridLogMessage << "HDCGDriver parameters ("
              << (pfile.length() ? pfile : std::string("defaults")) << "):" << std::endl;
    std::cout << P << std::endl;
  }
  const int cb = P.Checkerboard;

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
  SchurDiagMooeeOperator<MobiusFermionD,LatticeFermionD> HermOpEO(Ddwf);
  HermOpAdaptor<LatticeFermionD> FineOp(HermOpEO);     // Op = HermOp, for the coarsening

  //////////////////////////////////////////////////////////////////////
  // Grids -> coarsening.  Scope order is lifetime order.  The fp32 fine
  // grids and fermion operator serve the fp32 fine level of the
  // preconditioner (MultiGrid.Setup.FinePrecision, a run-time choice).
  //////////////////////////////////////////////////////////////////////
  typedef HDCGCoarsening<vSpinColourVector,CoarseScalar_t,NBASIS> Coarsening_t;
  std::cout << GridLogMessage << "Coarse sector precision (compiled): "
            << (sizeof(typename GridTypeMapper<CoarseScalar_t>::scalar_type)==sizeof(ComplexF) ? "fp32" : "fp64")
            << ", nbasis " << NBASIS << std::endl;

  MGCoarseGrids CGrids(FGrid, P.MultiGrid.Setup);
  MGFineGridsF  FGridsF(FGrid);

  LatticeGaugeFieldF UmuF(FGridsF.UGridF);
  precisionChange(UmuF,Umu);
  MobiusFermionF DdwfF(UmuF,*FGridsF.FGridF,*FGridsF.FrbGridF,*FGridsF.UGridF,*FGridsF.UrbGridF,P.Mass,P.M5,P.MobiusB,P.MobiusC);

  Coarsening_t Coarsening(CGrids, FGridsF, FrbGrid, cb, P.MultiGrid);

  Coarsening.GetSubspace(RNG5, FineOp);
  HDCGRefineSubspace(Coarsening, Ddwf, DdwfF, FineOp, P.Nrhs);    // no-op unless Refine is HDCG
  Coarsening.Coarsen(FineOp);
  Coarsening.CertifyCoarsening(FineOp);
  Coarsening.CoarseLanczos(P.Nrhs);

  //////////////////////////////////////////////////////////////////////
  // Solves: mrhs then (optionally) single RHS through the SAME objects.
  //////////////////////////////////////////////////////////////////////
  auto RunSolve = [&](int nr)
  {
    HDCGSolver<MobiusFermionD,MobiusFermionF,Coarsening_t> Solver(Ddwf,DdwfF,Coarsening,P.MultiGrid,nr);
    std::vector<LatticeFermionD> src(nr,FrbGrid), sol(nr,FrbGrid);
    for(int r=0;r<nr;r++){ src[r].Checkerboard()=cb; sol[r].Checkerboard()=cb; gaussian(RNG5,src[r]); sol[r]=Zero(); }
    Solver.Solve(src,sol);
    if ( nr == 1 ) {            // the same solve through the LinearFunction interface
      LatticeFermionD x(FrbGrid); x.Checkerboard()=cb; x=Zero();
      Solver.Solve(src[0],x);
    }
  };

  RunSolve(P.Nrhs);
  if ( P.SolveSingleRHS && P.Nrhs != 1 ) RunSolve(1);

  Grid_finalize();
}
