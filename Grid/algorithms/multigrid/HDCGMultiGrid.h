/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/HDCGMultiGrid.h

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
#pragma once

#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczos.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczosCoarse.h>
#include <Grid/algorithms/multigrid/PVdagMMultiGrid.h>
#include <Grid/algorithms/multigrid/HDCGMultiGridParams.h>

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////////////////
// The two-level mrhs HDCG (arXiv:1402.2585, mrhs form arXiv:2409.03904) as
// objects, the shape of PVdagMMultiGrid.h:
//
//   MGCoarseGrids / MGFineGridsF   shared with the PVdagM chain
//   HDCGCoarsening                 everything that is a function of the gauge
//                                  field: the raw near-null basis (and its
//                                  refinement), the transfer operators, the
//                                  Galerkin coarse operator (V2, Hermitian),
//                                  the coarse eigenvectors for deflation
//   HDCGSolver                     the solve chain on a borrowed coarsening:
//                                  smoother, deflated coarse solve, ADEF-2
//                                  preconditioner (fp64 or fp32 behind the
//                                  seam), outer mrhs fPcg / BlockCGrQ
//
// The Hermitian operator is the Schur-preconditioned M^dag M on one
// checkerboard; the near-null basis and the fine fields live on the
// checkerboarded grid.  The coarse space is unvectorised and ALWAYS the D+1
// multiRHS field; Nrhs 1 is the extent-1 case.  Grids outlive the
// coarsening outlives the solver.
//////////////////////////////////////////////////////////////////////////////////////
template<class Fobj,class CComplex,int nbasis>
class HDCGCoarsening {
public:
  typedef Lattice<Fobj>                                          FineField;
  typedef MultiGeneralCoarsenedOperatorV2<Fobj,CComplex,nbasis>  CoarseOperator;
  typedef typename CoarseOperator::CoarseVector                  CoarseVector;
  typedef typename GridTypeMapper<Fobj>::SinglePrecision         FobjF;
  typedef Lattice<FobjF>                                         FineFieldF;

  MGCoarseGrids                    &Grids;     // borrowed
  MGFineGridsF                     &GridsF;    // borrowed
  GridBase                         *FrbGrid;   // borrowed: the checkerboarded fine grid
  int                               cb;
  HDCGMultiGridParams               Params;
  NextToNextToNextToNearestStencilGeometry5D geom;   // the coarse stencil of the Schur operator
  CoarseOperator                    CoarseOp;
  //////////////////////////////////////////////////////////////////////
  // ONE fine transfer operator.  Its STORE follows the coarse sector --
  // the coarse space is what it feeds -- while its import and export
  // accept EITHER fine precision, because they are already a layout
  // transformation and a scalar conversion inside one costs nothing.
  //////////////////////////////////////////////////////////////////////
  static const bool CoarseIsSingle =
    ( sizeof(typename GridTypeMapper<CComplex>::scalar_type) == sizeof(ComplexF) );
  typedef typename std::conditional<CoarseIsSingle,FineFieldF,FineField>::type ProjectorField;
  typedef MultiRHSBlockProject<ProjectorField>                        ProjectorL1_t;

  ProjectorL1_t                     MrhsProjectorL1;
  std::vector<FineField>            rawNull;         // RAW fine near-null basis
  std::vector<CoarseVector>         evec;            // coarse eigenvectors, on Coarse5d
  std::vector<RealD>                eval;
  MultiRHSDeflation<CoarseVector>   Deflator;        // nev==0 until CoarseLanczos
  GridCartesian                    *CMrhs;           // transient solve grid, owned
  int                               nrhs;

  HDCGCoarsening(MGCoarseGrids &_Grids, MGFineGridsF &_GridsF, GridBase *_FrbGrid, int _cb,
                 const HDCGMultiGridParams &P)
    : Grids(_Grids), GridsF(_GridsF), FrbGrid(_FrbGrid), cb(_cb), Params(P),
      geom(_Grids.Coarse5d),
      CoarseOp(geom,_Grids.Coarse5d),
      CMrhs(nullptr),
      nrhs(-1)
  {
    Deflator.Deallocate();
    // The coarse stencil is a max-norm-1 box (Geometry.h) and the Schur
    // operator M^dag M reaches four fine hops, so a block shorter than 4
    // leaves couplings outside the stencil: the Fourier extraction then
    // folds them into the retained points and the Galerkin certificate
    // reads it.  Warn here; the certificate is the guard.
    for(int d=0;d<4;d++){
      if ( Params.Setup.Block1[d] < 4 )
        std::cout << GridLogWarning << "HDCGCoarsening: Block1[" << d << "] = " << Params.Setup.Block1[d]
                  << " < 4: the Schur operator's range exceeds the coarse stencil; "
                  << "expect a Galerkin defect" << std::endl;
    }
  };
  ~HDCGCoarsening()
  {
    CoarseOp.ReleaseGrid();       // borrowers let go before their grids die
    if ( CMrhs ) delete CMrhs;
  }

  int FileExists(const std::string &fname)
  {
    uint64_t exists=0;
    if ( fname.length() ) {
      if ( Grids.FGrid->IsBoss() ){ std::ifstream f(fname); exists=f.good()?1:0; }
      Grids.FGrid->GlobalSum(exists);
    }
    return (int)exists;
  }
  FineField FineNoise(GridParallelRNG &RNG)
  {
    FineField v(FrbGrid); v.Checkerboard()=cb; gaussian(RNG,v);
    return v;
  }
  // Chebyshev upper bound: the parameter, or 1.1 x a power-method estimate
  RealD UpperBound(RealD param, LinearOperatorBase<FineField> &HermOp, GridParallelRNG &RNG)
  {
    if ( param > 0.0 ) return param;
    FineField src = FineNoise(RNG);
    PowerMethod<FineField> PM;
    RealD lam = PM(HermOp,src);
    std::cout << GridLogMessage << "HDCGCoarsening: fine power method " << lam
              << ", Chebyshev upper bound " << 1.1*lam << std::endl;
    return 1.1*lam;
  }

  ////////////////////////////////////////////////////////////////////
  // The RAW basis: a refined basis on disk wins, then a raw basis on
  // disk, else built by the chosen method (and saved if named).  CG
  // refinement happens here; HDCG refinement needs a solver, see
  // HDCGRefineSubspace below.  The Aggregation is scaffolding for the
  // Chebyshev filters only: they run entirely on the fine grid.
  ////////////////////////////////////////////////////////////////////
  void GetSubspace(GridParallelRNG &RNG, LinearOperatorBase<FineField> &HermOp)
  {
    const HDCGSubspaceParams &S = Params.Subspace;
    rawNull.clear();
    rawNull.reserve(nbasis);
    for(int k=0;k<nbasis;k++){ rawNull.push_back(FineField(FrbGrid)); rawNull[k].Checkerboard()=cb; }

    if ( FileExists(S.RefinedSubspaceFile) ) {
      std::cout << GridLogMessage << "HDCGCoarsening: loading REFINED subspace " << S.RefinedSubspaceFile << std::endl;
      loadSubspace(rawNull, S.RefinedSubspaceFile, S.FileControl);
      return;
    }
    if ( FileExists(Params.Setup.SubspaceFile) ) {
      std::cout << GridLogMessage << "HDCGCoarsening: loading subspace " << Params.Setup.SubspaceFile << std::endl;
      loadSubspace(rawNull, Params.Setup.SubspaceFile, S.FileControl);
    } else {
      Aggregation<Fobj,CComplex,nbasis> Agg(Grids.Coarse5d,FrbGrid,cb);
      if ( S.Method == HDCGSubspaceMethod::ChebyshevNew ) {
        RealD hi = UpperBound(S.ChebyHi,HermOp,RNG);
        std::cout << GridLogMessage << "HDCGCoarsening: CreateSubspaceChebyshevNew hi " << hi << std::endl;
        Agg.CreateSubspaceChebyshevNew(RNG,HermOp,hi);
      } else if ( S.Method == HDCGSubspaceMethod::Chebyshev ) {
        RealD hi = UpperBound(S.ChebyHi,HermOp,RNG);
        std::cout << GridLogMessage << "HDCGCoarsening: CreateSubspaceChebyshev hi " << hi
                  << " lo " << S.ChebyLo << " order " << S.ChebyOrder << std::endl;
        Agg.CreateSubspaceChebyshev(RNG,HermOp,nbasis,hi,S.ChebyLo,S.ChebyOrder);
      } else {
        std::cout << GridLogMessage << "HDCGCoarsening: subspace from fine eigenvectors " << S.FineEvecFile
                  << " summed " << S.FineEvecSample << " per vector" << std::endl;
        LoadFineEvecSum(Agg.subspace,S.FineEvecFile,S.FineEvecSample,S.FileControl);
      }
      for(int k=0;k<nbasis;k++) rawNull[k]=Agg.subspace[k];
      if ( Params.Setup.SubspaceFile.length() )
        saveSubspace(rawNull, Params.Setup.SubspaceFile, S.FileControl);
    }
    if ( S.Refine == HDCGRefineMethod::CG ) {
      std::cout << GridLogMessage << "HDCGCoarsening: CG refinement shift " << S.RefineShift
                << " tol " << S.RefineTol << " maxit " << S.RefineMaxIt << std::endl;
      Aggregation<Fobj,CComplex,nbasis> Agg(Grids.Coarse5d,FrbGrid,cb);
      for(int k=0;k<nbasis;k++) Agg.subspace[k]=rawNull[k];
      Agg.RefineSubspace(HermOp,S.RefineShift,S.RefineTol,S.RefineMaxIt);
      for(int k=0;k<nbasis;k++) rawNull[k]=Agg.subspace[k];
      if ( S.RefinedSubspaceFile.length() )
        saveSubspace(rawNull, S.RefinedSubspaceFile, S.FileControl);
    }
  }
  // Basis vector b = precisionChange of the sum of `sample` consecutive fp32 records
  void LoadFineEvecSum(std::vector<FineField> &sub, const std::string &file, int sample, int control)
  {
#ifdef HAVE_LIME
    emptyUserRecord record;
    FineFieldF tmp(GridsF.FrbGridF), sum(GridsF.FrbGridF);
    tmp.Checkerboard()=cb; sum.Checkerboard()=cb;
    ScidacReader RD;
    RD.open(file);
    for(int b=0;b<nbasis;b++){
      sum = Zero();
      for(int n=0;n<sample;n++){
        RD.readScidacFieldRecord(tmp,record,control);
        sum = sum + tmp;
      }
      precisionChange(sub[b],sum);
    }
    RD.close();
#else
    GRID_ASSERT(0);
#endif
  }

  ////////////////////////////////////////////////////////////////////
  // HDCG refinement of the raw basis (the paper's loop): blocks of
  // nrhs vectors, normalised, solved with the refining solver, and
  // renormalised.  The solver was built on the SHIFTED operator and
  // the coarsening of it; the caller re-coarsens afterwards.
  ////////////////////////////////////////////////////////////////////
  void RefineSubspaceHDCG(LinearOperatorBase<FineField> &HermOp, TwoLevelCGmrhs<FineField> &Refiner, int nr)
  {
    std::vector<FineField> src(nr,FrbGrid), res(nr,FrbGrid);
    FineField tmp(FrbGrid); tmp.Checkerboard()=cb;
    for(int r=0;r<nr;r++){ src[r].Checkerboard()=cb; res[r].Checkerboard()=cb; }
    for(int b=0;b<nbasis;b+=nr){
      int nb = std::min(nbasis-b,nr);
      for(int r=0;r<nb;r++){
        RealD scale = std::pow(norm2(rawNull[b+r]),-0.5);
        src[r] = rawNull[b+r]*scale;
      }
      for(int r=nb;r<nr;r++) src[r] = src[0];   // pad the block; results discarded
      for(int r=0;r<nr;r++)  res[r] = Zero();
      HermOp.HermOp(src[0],tmp);
      std::cout << GridLogMessage << "HDCGCoarsening: refine block " << b << " before <n|A|n> " << norm2(tmp) << std::endl;
      Refiner(src,res);
      for(int r=0;r<nb;r++){
        RealD scale = std::pow(norm2(res[r]),-0.5);
        rawNull[b+r] = res[r]*scale;
      }
      HermOp.HermOp(rawNull[b],tmp);
      std::cout << GridLogMessage << "HDCGCoarsening: refine block " << b << " after  <n|A|n> " << norm2(tmp) << std::endl;
    }
  }

  ////////////////////////////////////////////////////////////////////
  // The Galerkin coarse operator and transfer operators from the RAW
  // basis.  CoarsenOperator block-orthonormalises its input in place,
  // so a working copy is taken and the RAW vectors survive.  The
  // operator handed in must apply the Hermitian operator through Op
  // (HermOpAdaptor).
  //
  // The RAW basis has done its work once this returns, so it is freed --
  // it is nbasis fine vectors, the single largest resident block of the
  // setup.  retain_basis keeps it, for the two callers that coarsen more
  // than once: the HDCG subspace refinement below, and a fixed-basis
  // rebuild on a changed gauge field (HMC).  Setup.RetainSubspace keeps
  // it for the whole life of the object.
  ////////////////////////////////////////////////////////////////////
  void Coarsen(LinearOperatorBase<FineField> &HermOp, int retain_basis=0)
  {
    GRID_ASSERT( rawNull.size() == nbasis );
    std::vector<FineField> sub(rawNull);

    CoarseOp.SetGrid(Grids.CoarseBatch);
    std::cout << GridLogMessage << "HDCGCoarsening: CoarsenOperator, batch " << Grids.batch << std::endl;
    // CoarsenOperator orthonormalises sub in place and imports it into OUR
    // transfer operator, which it then uses.  One basis store, imported once.
    // The projector's grid carries the STORE's SIMD layout: the fp32 fine
    // red-black grid when the coarse sector is fp32.
    MrhsProjectorL1.Allocate(nbasis,
                             CoarseIsSingle ? (GridBase *)GridsF.FrbGridF : (GridBase *)FrbGrid,
                             Grids.Coarse5d);
    CoarseOp.CoarsenOperator(HermOp,sub,Grids.Coarse5d,Grids.batch,MrhsProjectorL1);
    if ( !retain_basis && !Params.Setup.RetainSubspace ) DiscardBasis();
    MrhsProjectorL1.ReleaseScratch();
    nrhs = -1;                     // the batch grid is set; SetNrhs retargets
    ReportDeviceFootprint("after Coarsen");
  }

  ////////////////////////////////////////////////////////////////////
  // Retarget the coarse operator to Nrhs: the D+1 grid is owned here
  // and every borrower re-captures it after this call.
  ////////////////////////////////////////////////////////////////////
  void SetNrhs(int nr)
  {
    if ( nr == nrhs ) return;
    CoarseOp.ReleaseGrid();
    if ( CMrhs ) { delete CMrhs; CMrhs=nullptr; }
    Coordinate cml({nr,1,Grids.clatt[0],Grids.clatt[1],Grids.clatt[2],Grids.clatt[3]});
    CMrhs = new GridCartesian(cml,Grids.cmsimd,Grids.cmmpi);
    CoarseOp.SetGrid(CMrhs);
    nrhs = nr;
    std::cout << GridLogMessage << "HDCGCoarsening: coarse operator at Nrhs " << nr << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // Coarse eigenvectors for the deflated coarse solve.  Block Lanczos
  // runs at the solve's Nrhs on the D+1 operator; plain Lanczos at
  // Nrhs 1 (the extent-1 D+1 grid), its vectors then sliced to the D
  // grid so both feed the deflator the same way.  Loaded from EvecFile
  // when present, saved there when named.
  ////////////////////////////////////////////////////////////////////
  void CoarseLanczos(int nr)
  {
    const HDCGLanczosParams &L = Params.Lanczos;
    evec.clear(); eval.clear();
    Deflator.Deallocate();
    if ( L.Method == HDCGLanczosMethod::None ) {
      std::cout << GridLogMessage << "HDCGCoarsening: no coarse deflation" << std::endl;
      return;
    }
    if ( FileExists(L.EvecFile) ) {
      std::cout << GridLogMessage << "HDCGCoarsening: loading coarse eigenpairs " << L.EvecFile << std::endl;
      loadEigenpairs(evec,eval,L.EvecFile,Grids.Coarse5d,Params.Subspace.FileControl);
    } else {
      int Nconv=0;
      SetNrhs( L.Method == HDCGLanczosMethod::Block ? nr : 1 );
      HermitianLinearOperator<CoarseOperator,CoarseVector> HermOpC(CoarseOp);
      GridParallelRNG cRNG(CMrhs); cRNG.SeedFixedIntegers(std::vector<int>({31,32,33,34}));
      RealD hi = L.ChebyHi;
      if ( hi <= 0.0 ) {
        CoarseVector s(CMrhs); random(cRNG,s);
        PowerMethod<CoarseVector> PM;
        hi = 1.1*PM(HermOpC,s);
        std::cout << GridLogMessage << "HDCGCoarsening: coarse Chebyshev upper bound " << hi << std::endl;
      }
      Chebyshev<CoarseVector> Cheby(L.ChebyLo,hi,L.ChebyOrder);

      if ( L.Method == HDCGLanczosMethod::Block ) {
        std::cout << GridLogMessage << "HDCGCoarsening: block Lanczos, Nrhs " << nr
                  << " Nstop " << L.Nstop << " Nk " << L.Nk << " Nm " << L.Nm << std::endl;
        ImplicitlyRestartedBlockLanczosCoarse<CoarseVector>
          IRL(HermOpC,Grids.Coarse5d,CMrhs,nr,Cheby,L.Nstop,1,nr,L.Nk,L.Nm,L.Tol,L.MaxIt);
        evec.resize(L.Nm,Grids.Coarse5d);
        eval.resize(L.Nm);
        std::vector<CoarseVector> src(nr,Grids.Coarse5d);
        GridParallelRNG dRNG(Grids.Coarse5d); dRNG.SeedFixedIntegers(std::vector<int>({41,42,43,44}));
        for(int r=0;r<nr;r++) random(dRNG,src[r]);
        IRL.calc(eval,evec,src,Nconv,LanczosType::irbl);
        evec.resize(Nconv,Grids.Coarse5d);                 // only the converged pairs
        eval.resize(Nconv);
      } else {
        std::cout << GridLogMessage << "HDCGCoarsening: plain Lanczos at Nrhs 1, Nstop " << L.Nstop
                  << " Nk " << L.Nk << " Nm " << L.Nm << std::endl;
        PlainHermOp<CoarseVector>    Op(HermOpC);
        FunctionHermOp<CoarseVector> OpCheby(Cheby,HermOpC);
        ImplicitlyRestartedLanczos<CoarseVector> IRL(OpCheby,Op,L.Nstop,L.Nk,L.Nm,L.Tol,L.MaxIt);
        std::vector<CoarseVector> evec1(L.Nm,CMrhs);
        std::vector<RealD>        eval1(L.Nm);
        CoarseVector src(CMrhs); random(cRNG,src);
        IRL.calc(eval1,evec1,src,Nconv);
        evec.resize(Nconv,Grids.Coarse5d);
        eval.resize(Nconv);
        for(int e=0;e<Nconv;e++){ ExtractSliceFast(evec[e],evec1[e],0,0); eval[e]=eval1[e]; }
      }
      // The Lanczos routines set their outputs only on reaching Nstop;
      // a partial result carries no usable eigenpairs.
      if ( Nconv < L.Nstop ) {
        std::cout << GridLogWarning << "HDCGCoarsening: coarse Lanczos converged " << Nconv << " < Nstop " << L.Nstop
                  << " modes (Chebyshev window [" << L.ChebyLo << "," << hi << "] must hold at least Nstop eigenvalues"
                  << " and Nm must resolve them); coarse solve undeflated" << std::endl;
        evec.clear(); eval.clear();
        return;
      }
      if ( L.EvecFile.length() ) saveEigenpairs(evec,eval,L.EvecFile,Params.Subspace.FileControl);
    }
    const int nev = evec.size();
    Deflator.ImportEigenBasis(evec,eval);
    // The Lattice-layout copies have been imported into the BLAS store and
    // saved if they were going to be; the solve reads only the BLAS store.
    evec.clear(); evec.shrink_to_fit();
    std::cout << GridLogMessage << "HDCGCoarsening: " << nev << " coarse eigenvectors imported for deflation"
              << (nev ? " (lowest "+std::to_string(eval[0])+", highest "+std::to_string(eval[nev-1])+")" : "")
              << std::endl;
    ReportDeviceFootprint("after CoarseLanczos");
  }

  ////////////////////////////////////////////////////////////////////
  // Certificates at Nrhs 1 through the solve's own transfer operators:
  // Galerkin ||A_c x - P^dag A P x|| / ||A_c x||  (the rounding level of
  // the coarsening: fp64 coarse ~1e-13, fp32 ~1e-6), and Hermiticity
  // |<x,A_c y> - <A_c x,y>| / |<x,A_c y>|, which CG relies on.
  ////////////////////////////////////////////////////////////////////
  RealD CertifyCoarsening(LinearOperatorBase<FineField> &HermOp)
  {
    SetNrhs(1);
    CoarseVector xc(CMrhs), yc(CMrhs), Acx(CMrhs), Acy(CMrhs), PtAPx(CMrhs);
    GridParallelRNG cRNG(CMrhs); cRNG.SeedFixedIntegers(std::vector<int>({21,22,23,24}));
    random(cRNG,xc);
    random(cRNG,yc);

    CoarseOp.M(xc,Acx);                                     // A_c x
    CoarseOp.M(yc,Acy);                                     // A_c y

    std::vector<FineField> Px(1,FrbGrid), APx(1,FrbGrid);
    Px[0].Checkerboard()=cb; APx[0].Checkerboard()=cb;
    MrhsProjectorL1.blockPromote(Px,xc);                      // P x
    HermOp.HermOp(Px[0],APx[0]);                            // A P x
    MrhsProjectorL1.blockProject(APx,PtAPx);                  // P^dag A P x

    CoarseVector d(CMrhs); d = Acx - PtAPx;
    RealD rel = std::sqrt(norm2(d)/norm2(Acx));
    ComplexD xAy = TensorRemove(innerProduct(xc,Acy));
    ComplexD Axy = TensorRemove(innerProduct(Acx,yc));
    RealD herm = abs(xAy-Axy)/abs(xAy);
    std::cout << GridLogMessage << "HDCGCoarsening: GALERKIN CERTIFICATE ||A_c x - P^dag A P x||/||A_c x|| = " << rel << std::endl;
    std::cout << GridLogMessage << "HDCGCoarsening: HERMITICITY CERTIFICATE |<x,A_c y> - <A_c x,y>|/|<x,A_c y>| = " << herm << std::endl;
    return rel;
  }


  //////////////////////////////////////////////////////////////////////
  // The device memory this coarsening holds that the memory manager
  // CANNOT evict: BLAS stores, which are plain device allocations.
  // Lattice fields are omitted on purpose -- they are evictable, so they
  // cost eviction traffic, not failure.  What is reported is what has to
  // fit alongside the manager's cache, the comms buffers and the shm
  // segment, which is the budget an out-of-memory run actually breaks.
  //////////////////////////////////////////////////////////////////////
  void ReportDeviceFootprint(const std::string &when)
  {
    uint64_t p1 = MrhsProjectorL1.DeviceBytes();
    uint64_t o1 = CoarseOp.DeviceBytes();
    uint64_t df = Deflator.DeviceBytes();
    const double G = 1024.*1024.*1024.;
    std::cout << GridLogMessage << "HDCGCoarsening (" << when << "): resident device memory/rank "
              << (p1+o1+df)/G << " GB = transfer L1 " << p1/G
              << " + operator L1 " << o1/G
              << " + deflation " << df/G << " GB" << std::endl;
    std::cout << GridLogMessage << "HDCGCoarsening (" << when << "): not evictable; must fit alongside the "
              << MemoryManager::DeviceMaxBytes/G
              << " GB manager cache, the comms buffers and the shm segment" << std::endl;
  }

  void DiscardBasis(void)
  {
    rawNull.clear(); rawNull.shrink_to_fit();
  }
};

//////////////////////////////////////////////////////////////////////
// The solve chain on a borrowed coarsening.  Both fine precisions are
// built; the outer picks the fp32 preconditioner through the seam when
// Setup.FinePrecision is fp32.  OperatorShift is 0 in production; the
// HDCG-refinement solver runs on the shifted operator.
//////////////////////////////////////////////////////////////////////
template<class Matrix,class MatrixF,class Coarsening>
class HDCGSolver {
public:
  typedef typename Coarsening::FineField      FineField;
  typedef typename Coarsening::FineFieldF     FineFieldF;
  typedef typename Coarsening::CoarseOperator CoarseOperator;
  typedef typename Coarsening::CoarseVector   CoarseVector;

  Coarsening          &C;
  HDCGMultiGridParams  Params;
  int                  nrhs;
  int                  _regrid;     // FIRST member-like init: SetNrhs before members capture grids
  RealD                OperatorShift;

  // fp64 fine chain
  SchurDiagMooeeOperator<Matrix,FineField>    HermOpEO;
  ShiftedHermOpLinearOperator<FineField>      HermOpShifted;
  HermOpAdaptor<FineField>                    FineOp;         // Op = HermOp: the operator of the solve
  ShiftedHermOpLinearOperator<FineField>      SmootherOp;
  CGSmoother<FineField>                       SmootherCG;
  ChebyshevSmoother<FineField>                SmootherCheby;
  // coarse
  HermitianLinearOperator<CoarseOperator,CoarseVector> LinOpC;
  ConjugateGradient<CoarseVector>             CoarseCG;
  DoNothingGuesser<CoarseVector>              NoGuess;
  HPDSolver<CoarseVector>                     CoarseDeflatedCG;
  ChebyshevInverter<CoarseVector>             CoarseCheby;
  MrhsADEF2Preconditioner<FineField,CoarseVector,typename Coarsening::ProjectorL1_t> ADEF2;
  // fp32 fine chain, behind the seam; the coarse solve is shared
  SchurDiagMooeeOperator<MatrixF,FineFieldF>  HermOpEOF;
  ShiftedHermOpLinearOperator<FineFieldF>     HermOpShiftedF;
  HermOpAdaptor<FineFieldF>                   FineOpF;
  ShiftedHermOpLinearOperator<FineFieldF>     SmootherOpF;
  CGSmoother<FineFieldF>                      SmootherCGF;
  ChebyshevSmoother<FineFieldF>               SmootherChebyF;
  MrhsADEF2Preconditioner<FineFieldF,CoarseVector,typename Coarsening::ProjectorL1_t> ADEF2F;
  MrhsMixedPrecPreconditioner<FineField,FineFieldF> PreconSeam;
  TwoLevelCGmrhs<FineField>                   Outer;

  LinearFunction<FineField>  &PickSmoother(void) {
    return (Params.Smoother.Type==HDCGSmootherType::CG)
      ? static_cast<LinearFunction<FineField>&>(SmootherCG)
      : static_cast<LinearFunction<FineField>&>(SmootherCheby);
  }
  LinearFunction<FineFieldF> &PickSmootherF(void) {
    return (Params.Smoother.Type==HDCGSmootherType::CG)
      ? static_cast<LinearFunction<FineFieldF>&>(SmootherCGF)
      : static_cast<LinearFunction<FineFieldF>&>(SmootherChebyF);
  }
  LinearFunction<CoarseVector> &PickCoarseSolve(void) {
    return (Params.CoarseSolver.Method==HDCGCoarseMethod::DeflatedCG)
      ? static_cast<LinearFunction<CoarseVector>&>(CoarseDeflatedCG)
      : static_cast<LinearFunction<CoarseVector>&>(CoarseCheby);
  }
  template<class F, class Op, class ShiftedOp>
  static LinearOperatorBase<F> &Unshifted(RealD shift, Op &op, ShiftedOp &shifted) {
    return (shift==0.0) ? static_cast<LinearOperatorBase<F>&>(op) : static_cast<LinearOperatorBase<F>&>(shifted);
  }

  HDCGSolver(Matrix &Ddwf, MatrixF &DdwfF, Coarsening &_C, const HDCGMultiGridParams &P, int nr, RealD shift=0.0)
    : C(_C), Params(P), nrhs(nr),
      _regrid((C.SetNrhs(nr),0)),
      OperatorShift(shift),
      HermOpEO(Ddwf),
      HermOpShifted(HermOpEO,shift),
      FineOp(Unshifted<FineField>(shift,HermOpEO,HermOpShifted)),
      SmootherOp(FineOp,P.Smoother.Shift),
      SmootherCG(P.Smoother.Order,SmootherOp),
      SmootherCheby(P.Smoother.ChebyLo,P.Smoother.ChebyHi,P.Smoother.Order,FineOp),
      LinOpC(C.CoarseOp),
      CoarseCG(P.CoarseSolver.Tol,P.CoarseSolver.MaxIt,false),
      CoarseDeflatedCG(LinOpC,CoarseCG,NoGuess),
      CoarseCheby(P.CoarseSolver.ChebyLo,P.CoarseSolver.ChebyHi,P.CoarseSolver.ChebyOrder,LinOpC),
      ADEF2(FineOp,PickSmoother(),PickCoarseSolve(),PickCoarseSolve(),C.MrhsProjectorL1,C.Deflator,C.CMrhs,C.FrbGrid),
      HermOpEOF(DdwfF),
      HermOpShiftedF(HermOpEOF,shift),
      FineOpF(Unshifted<FineFieldF>(shift,HermOpEOF,HermOpShiftedF)),
      SmootherOpF(FineOpF,P.Smoother.Shift),
      SmootherCGF(P.Smoother.Order,SmootherOpF),
      SmootherChebyF(P.Smoother.ChebyLo,P.Smoother.ChebyHi,P.Smoother.Order,FineOpF),
      ADEF2F(FineOpF,PickSmootherF(),PickCoarseSolve(),PickCoarseSolve(),C.MrhsProjectorL1,C.Deflator,C.CMrhs,C.GridsF.FrbGridF),
      PreconSeam(ADEF2F,C.FrbGrid,C.GridsF.FrbGridF,nr),
      Outer(P.Outer.Tol,P.Outer.MaxIt,FineOp,
            (P.Setup.FinePrecision==MGPrecision::fp32)
              ? static_cast<MrhsPreconditioner<FineField>&>(PreconSeam)
              : static_cast<MrhsPreconditioner<FineField>&>(ADEF2),
            C.FrbGrid,
            (P.Outer.Algorithm==HDCGOuterAlgorithm::PrecBlockCGrQ) ? MrhsCGAlgorithm::PrecBlockCGrQ : MrhsCGAlgorithm::fPcg)
  {
    // BlockCGrQ assumes a stationary preconditioner: the fixed-iteration CG
    // smoother is fixed WORK but a nonlinear map, and a coarse CG solved to a
    // tolerance is not stationary either.  Chebyshev on both is.
    if ( Params.Outer.Algorithm==HDCGOuterAlgorithm::PrecBlockCGrQ &&
         ( Params.CoarseSolver.Method!=HDCGCoarseMethod::Chebyshev ||
           Params.Smoother.Type     !=HDCGSmootherType::Chebyshev ) )
      std::cout << GridLogWarning << "HDCGSolver: PrecBlockCGrQ with a non-stationary preconditioner "
                << "(needs Smoother.Type Chebyshev AND CoarseSolver.Method Chebyshev); expect stagnation"
                << std::endl;
    std::cout << GridLogMessage << "HDCGSolver: Nrhs " << nr
              << ", operator shift " << shift
              << ", smoother " << (Params.Smoother.Type==HDCGSmootherType::CG ? "CG" : "Chebyshev")
              << " order " << Params.Smoother.Order
              << ", coarse " << (Params.CoarseSolver.Method==HDCGCoarseMethod::DeflatedCG ? "deflated CG" : "Chebyshev")
              << " (deflation vectors " << C.Deflator.nev << ")"
              << ", preconditioner fine level "
              << (Params.Setup.FinePrecision==MGPrecision::fp32 ? "fp32 (seam at the outer CG)" : "fp64")
              << ", outer " << (Params.Outer.Algorithm==HDCGOuterAlgorithm::PrecBlockCGrQ ? "PrecBlockCGrQ" : "fPcg")
              << " EXACT fp64" << std::endl;
  }

  // One right-hand side through the LinearFunction interface (nrhs must be 1)
  void Solve(const FineField &src, FineField &sol)
  {
    GRID_ASSERT( nrhs == 1 );
    sol.Checkerboard() = src.Checkerboard();
    GridStopWatch w; w.Start();
    Outer(src,sol);
    w.Stop();
    FineField Ax(src.Grid()); Ax.Checkerboard()=src.Checkerboard();
    FineOp.HermOp(sol,Ax);
    Ax=Ax-src;
    std::cout << GridLogMessage << "HDCGSolver: single rhs total " << w.Elapsed()
              << "; FINAL Nrhs 1 (LinearFunction): true residual = " << std::sqrt(norm2(Ax)/norm2(src)) << std::endl;
  }
  void Solve(std::vector<FineField> &src, std::vector<FineField> &sol)
  {
    GRID_ASSERT( (int)src.size() == nrhs );
    for(int r=0;r<nrhs;r++) sol[r].Checkerboard() = src[r].Checkerboard();
    MemoryManager::EvictAll();
    MemoryManager::DropCache();

    GridStopWatch w; w.Start();
    Outer(src,sol);
    w.Stop();
    std::cout << GridLogMessage << "HDCGSolver: Nrhs "<<nrhs<<" total " << w.Elapsed()
              << "  (per RHS: " << w.useconds()/1.0e6/nrhs << " s)" << std::endl;
    {
      FineField Ax(src[0].Grid()); Ax.Checkerboard()=src[0].Checkerboard();
      RealD worst=0.0;
      for(int r=0;r<nrhs;r++){
        FineOp.HermOp(sol[r],Ax);
        Ax=Ax-src[r];
        RealD rn=std::sqrt(norm2(Ax)/norm2(src[r]));
        std::cout << GridLogMessage << "FINAL Nrhs "<<nrhs<<": rhs["<<r<<"] true residual = " << rn << std::endl;
        worst=std::max(worst,rn);
      }
      std::cout << GridLogMessage << "FINAL Nrhs "<<nrhs<<": worst-case residual = " << worst << std::endl;
    }
  }
};

//////////////////////////////////////////////////////////////////////
// HDCG refinement of the basis (the paper's loop): coarsen the SHIFTED
// operator, deflate it, solve each raw basis vector on it with a
// short-tolerance HDCG, keep the normalised solutions.  A no-op unless
// Subspace.Refine is HDCG and no refined basis was loaded.  The caller
// then coarsens the unshifted operator as usual.
//////////////////////////////////////////////////////////////////////
template<class Matrix,class MatrixF,class Coarsening>
void HDCGRefineSubspace(Coarsening &C, Matrix &Ddwf, MatrixF &DdwfF,
                        LinearOperatorBase<typename Coarsening::FineField> &HermOp, int nr)
{
  typedef typename Coarsening::FineField FineField;
  const HDCGSubspaceParams &S = C.Params.Subspace;
  if ( S.Refine != HDCGRefineMethod::HDCG ) return;
  if ( C.FileExists(S.RefinedSubspaceFile) ) return;

  std::cout << GridLogMessage << "HDCGRefineSubspace: shift " << S.RefineShift << " tol " << S.RefineTol
            << " maxit " << S.RefineMaxIt << " smoother order " << S.RefineSmootherOrder << std::endl;
  ShiftedHermOpLinearOperator<FineField> Shifted(HermOp,S.RefineShift);
  HermOpAdaptor<FineField>               ShiftedOp(Shifted);
  C.Coarsen(ShiftedOp,1);   // refinement below still needs the RAW basis
  C.CoarseLanczos(nr);

  HDCGMultiGridParams PR = C.Params;
  PR.Outer.Tol      = S.RefineTol;
  PR.Outer.MaxIt    = S.RefineMaxIt;
  PR.Smoother.Order = S.RefineSmootherOrder;
  HDCGSolver<Matrix,MatrixF,Coarsening> Refiner(Ddwf,DdwfF,C,PR,nr,S.RefineShift);
  C.RefineSubspaceHDCG(HermOp,Refiner.Outer,nr);

  if ( S.RefinedSubspaceFile.length() )
    saveSubspace(C.rawNull, S.RefinedSubspaceFile, S.FileControl);
}

NAMESPACE_END(Grid);
