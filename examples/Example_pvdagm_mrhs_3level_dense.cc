/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_mrhs_3level_dense.cc

    Copyright (C) 2026

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */

// MultiRHS (valence) THREE-level multigrid for PVdagM with a DENSE, EXACT,
// non-iterative coarse-coarse bottom  ==  Example_pvdagm_mrhs_3level.cc with
// the L3 PGCR replaced by DistributedDenseInverse (the machinery validated in
// Example_pvdagm_3level_dense.cc: VERIFY ~6e-5, apply ~6ms single-RHS).
//
// The two headline results COMPOSED: mrhs GEMM batching at the coarse level
// (the ~12x valence win) x supercoarse-dense exact bottom (the coarsest level
// reduced from dominant cost to noise).  Under mrhs the dense apply becomes a
// GEMM against the resident row-slab: ONE GlobalSum of nrhs packed vectors +
// ONE slab pass for all rhs -- per-rhs bottom cost BELOW the single-RHS 6ms.
//
// Level structure:
//   L1 (fine)         : std::vector<LatticeFermionD>, MrhsPGCRNonHermitian on PVdagM,
//                       preconditioned by the L1->L2 mrhs V-cycle (MrhsTwoLevelMG).
//   L2 (coarse)       : 6D mrhs coarse field, PGCR, preconditioned by the L2->L3 mrhs
//                       V-cycle -- DENSE coarse-coarse correction + coarse smoother.
//   L3 (coarse-coarse): DENSE row-distributed A^{-1} (cgetrf_64 + blocked cgetrs_64,
//                       fp32, ILP64 -- rank N = gSites(CC) x nbasis, e.g. 69120 at
//                       BLOCK2=8.4.3.6 x nb60).  DENSE_CC=0 reverts to the L3 PGCR.
//
// SUPERCOARSE default: BLOCK2 = 8,4,3,6  (CC = [3,6,8,8], the dense floor on
// --mpi 3.6.4.4; sigma spectrum is FLAT across all 60 vectors => keep nb60).
//
// RAW-NULL DISCIPLINE unchanged (see project_block_orthogonalise_leak): guards
// print ||<psi|psi> - I||_F at both levels (~0.23 good, ~sqrt(N) = e_k leak).
//
// Build (Frontier/HIP): LIBS += -lrocsolver -lrocblas.
// Run: --device-mem 40000 (solve-phase residency; EvictAll makes the setup
//      window; boss GCD holds 38.2GB dA + 141MB dB at BLOCK2=8.4.3.6 nb60).
//
// Env: MASS SUBSPACE_FILE NRHS
//      BLOCK (dotted, default 2.2.2.2)  BLOCK2 (dotted, default 8.4.3.6)
//      FineSmootherShift FineSmootherOrder
//      CoarseSmootherShift CoarseSmootherNstep
//      CoarseSolverTol CoarseSolverOrder
//      DENSE_CC (default 1) DENSE_CC_CHECK
//      L3_TOL L3_MAXIT L3_NSTEP   (iterative branch only)
//      OuterMmax OuterNstep OuterTol

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>

#include <unordered_map>
#include <memory>

#ifdef GRID_HIP
#include <rocsolver/rocsolver.h>
#endif

using namespace std;
using namespace Grid;

RealD FineSmootherShift    = 0.1;
int   FineSmootherOrder    = 16;
RealD CoarseSmootherShift  = 0.1;
int   CoarseSmootherNstep  = 4;
RealD CoarseSolverTol      = 0.03;
int   CoarseSolverOrder    = 200;
RealD L3Tol                = 2.5e-1;
int   L3MaxIt              = 50;
int   L3Nstep              = 50;
RealD OuterTol             = 1.0e-8;
int   OuterMmax            = 8;
int   OuterNstep           = 8;
int   Nrhs                 = 12;
int   UseDenseCC           = 1;
RealD mass                 = 0.00078;

void ParseEnvironment(void)
{
  if(getenv("MASS"))               mass               = atof(getenv("MASS"));
  if(getenv("FineSmootherShift"))  FineSmootherShift  = atof(getenv("FineSmootherShift"));
  if(getenv("FineSmootherOrder"))  FineSmootherOrder  = atoi(getenv("FineSmootherOrder"));
  if(getenv("CoarseSmootherShift"))CoarseSmootherShift= atof(getenv("CoarseSmootherShift"));
  if(getenv("CoarseSmootherNstep"))CoarseSmootherNstep= atoi(getenv("CoarseSmootherNstep"));
  if(getenv("CoarseSolverTol"))    CoarseSolverTol    = atof(getenv("CoarseSolverTol"));
  if(getenv("CoarseSolverOrder"))  CoarseSolverOrder  = atoi(getenv("CoarseSolverOrder"));
  if(getenv("L3_TOL"))             L3Tol              = atof(getenv("L3_TOL"));
  if(getenv("L3_MAXIT"))           L3MaxIt            = atoi(getenv("L3_MAXIT"));
  if(getenv("L3_NSTEP"))           L3Nstep            = atoi(getenv("L3_NSTEP"));
  if(getenv("OuterTol"))           OuterTol           = atof(getenv("OuterTol"));
  if(getenv("OuterMmax"))          OuterMmax          = atoi(getenv("OuterMmax"));
  if(getenv("OuterNstep"))         OuterNstep         = atoi(getenv("OuterNstep"));
  if(getenv("NRHS"))               Nrhs               = atoi(getenv("NRHS"));
  if(getenv("DENSE_CC"))           UseDenseCC         = atoi(getenv("DENSE_CC"));

  std::cout << GridLogMessage << "PARAM: MASS               " << mass               << std::endl;
  std::cout << GridLogMessage << "PARAM: NRHS               " << Nrhs               << std::endl;
  std::cout << GridLogMessage << "PARAM: DENSE_CC           " << UseDenseCC         << std::endl;
  std::cout << GridLogMessage << "PARAM: FineSmootherShift  " << FineSmootherShift  << std::endl;
  std::cout << GridLogMessage << "PARAM: FineSmootherOrder  " << FineSmootherOrder  << std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSmootherShift" << CoarseSmootherShift<< std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSmootherNstep" << CoarseSmootherNstep<< std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverTol    " << CoarseSolverTol    << std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverOrder  " << CoarseSolverOrder  << std::endl;
  std::cout << GridLogMessage << "PARAM: L3_TOL             " << L3Tol              << std::endl;
  std::cout << GridLogMessage << "PARAM: OuterMmax          " << OuterMmax          << std::endl;
  std::cout << GridLogMessage << "PARAM: OuterNstep         " << OuterNstep         << std::endl;
}

template <class Field>
void saveSubspace(std::vector<Field> &subspace, std::string const fname){
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacWriter SW(subspace[0].Grid()->IsBoss());
  SW.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++) SW.writeScidacFieldRecord(subspace[k], record);
  SW.close();
#endif
}
template <class Field>
void loadSubspace(std::vector<Field> &subspace, std::string const fname){
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacReader SR;
  SR.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++) SR.readScidacFieldRecord(subspace[k], record);
  SR.close();
#endif
}

//////////////////////////////////////////////////////////////////////
// A = PV^dag M (non-Hermitian), and shifted variant for smoothers.
//////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat; Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV) {};
  void OpDiag (const Field &in, Field &out) { assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll  (const Field &in, std::vector<Field> &out){ assert(0); };
  void Op     (const Field &in, Field &out){ Field tmp(in.Grid()); _Mat.M(in,tmp); _PV.Mdag(tmp,out); }
  void AdjOp  (const Field &in, Field &out){ Field tmp(in.Grid()); _PV.M(in,tmp); _Mat.Mdag(tmp,out); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ HermOp(in,out); ComplexD d=innerProduct(in,out); n1=real(d); n2=norm2(out); }
  void HermOp(const Field &in, Field &out){ Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

template<class Matrix,class Field>
class ShiftedPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat; Matrix &_PV;
public:
  RealD shift;
  ShiftedPVdagMLinearOperator(RealD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV){};
  void OpDiag (const Field &in, Field &out) { assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll  (const Field &in, std::vector<Field> &out){ assert(0); };
  void Op     (const Field &in, Field &out){ Field tmp(in.Grid()); _Mat.M(in,tmp); _PV.Mdag(tmp,out); out = out + shift*in; }
  void AdjOp  (const Field &in, Field &out){ Field tmp(in.Grid()); _PV.M(tmp,out); _Mat.Mdag(in,tmp); out = out + shift*in; }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp(const Field &in, Field &out){ Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

// Generic shift wrapper (for the coarse-level smoother on the 6D mrhs coarse operator).
template<class Field>
class ShiftedLinearOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Op; RealD shift;
public:
  ShiftedLinearOperator(RealD _shift, LinearOperatorBase<Field> &Op) : _Op(Op), shift(_shift) {}
  void OpDiag  (const Field &in, Field &out) { assert(0); }
  void OpDir   (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll (const Field &in, std::vector<Field> &out) { assert(0); }
  void Op      (const Field &in, Field &out) { _Op.Op(in,out);    out = out + shift*in; }
  void AdjOp   (const Field &in, Field &out) { _Op.AdjOp(in,out); out = out + shift*in; }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp  (const Field &in, Field &out) { Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

//////////////////////////////////////////////////////////////////////////////////////
// DistributedDenseInverse: exact solve of a (small) coarse operator by explicit,
// row-distributed dense inverse.  VALIDATED single-RHS in Example_pvdagm_3level_dense.cc
// (VERIFY ~6e-5).  Setup: probe columns -> exact chunked GlobalSum gather streamed
// into the boss-GCD device buffer -> cgetrf_64 (ILP64) -> rows of A^{-1} produced
// blockwise via cgetrs_64 on identity-column blocks and broadcast; each rank keeps
// the rows of ITS OWN sites (ownership-aligned => apply needs ONE GlobalSum only).
//
// NEW HERE: ApplyBatch -- the mrhs form.  ONE GlobalSum of all nrhs packed vectors,
// ONE pass over the slab computing all rhs (GEMM shape).  Per-rhs cost < single-RHS.
//
// Tensor-depth agnostic: site scalar_object treated as nbasis contiguous ComplexD.
//////////////////////////////////////////////////////////////////////////////////////
template<class Field>
class DistributedDenseInverse : public LinearFunction<Field> {
public:
  using LinearFunction<Field>::operator();
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object  sobj;

  GridBase *grid;
  LinearOperatorBase<Field> &_Op;   // kept for env-gated defect checks (DENSE_CC_CHECK)
  int      nbasis;
  int64_t  N;               // dense rank = gSites * nbasis
  int      nd;
  int      lsites;          // my local sites
  int      nrows;           // my rows = lsites * nbasis
  std::vector<Coordinate> myLcoor;      // local coordinate of my site ss
  std::vector<int64_t>    myGsite;      // global lex site index of my site ss
  std::vector<ComplexF>   slab;         // nrows x N row-major rows of A^{-1}
  // CHUNKROWS=1024: fatter cgetrs_64 calls (trsm efficiency grows with nrhs;
  // 256 gave ~7s/chunk => 32min setup) and 4x fewer broadcasts.  Buffers 566MB.
  static const int64_t CHUNKROWS = 1024;
  static const int MRHS_MAX = 32;
#ifdef GRID_HIP
  rocblas_float_complex *dSlab = nullptr;  // device-resident copy of the slab (all ranks)
  rocblas_handle          applyHandle;     // per-rank handle for the apply GEMM
  // Persistent apply buffers, hoisted out of the per-call path: pinned host
  // staging (full PCIe rate, no page faults) + device X/Y (no per-call
  // hipMalloc/hipFree).  Sized once for MRHS_MAX.
  ComplexF              *hXapply = nullptr;
  ComplexF              *hYapply = nullptr;
  rocblas_float_complex *dXapply = nullptr;
  rocblas_float_complex *dYapply = nullptr;
  int devSum  = 0;    // DENSE_DEVICE_SUM=1: allreduce the DEVICE buffer
                      // (GPU-aware MPI, no host staging); default host sum
#endif
  // NB: the apply GEMM (Y = slab^dag X) is a tiny-output/huge-K shape that
  // under-fills the GPU (~13ms).  The fix is a software split-K via
  // GridBLAS.gemmBatched (see MultiRHSBlockCGLinalg.h / 2409.03904 Fig 11) —
  // NOT a raw rocblas strided-batched batch, which hung on Frontier and was
  // removed.  TODO: reimplement through GridBLAS when the ~1.4 s/RHS is wanted.

  DistributedDenseInverse(LinearOperatorBase<Field> &Op, GridBase *g, int nbasis_)
    : grid(g), _Op(Op), nbasis(nbasis_)
  {
    GRID_ASSERT( sizeof(sobj) == nbasis*sizeof(ComplexD) ); // site object == nbasis ComplexD
    nd     = grid->_ndimension;
    N      = grid->gSites() * nbasis;
    lsites = grid->lSites();
    nrows  = lsites * nbasis;

    Coordinate ldims = grid->LocalDimensions();
    Coordinate gdims = grid->GlobalDimensions();

    std::cout << GridLogMessage << "DistributedDenseInverse: N = " << N
              << " (" << grid->gSites() << " sites x " << nbasis << ")"
              << "  rows/rank = " << nrows
              << "  slab = " << (double)nrows*N*sizeof(ComplexF)/1024./1024. << " MB/rank"
              << std::endl;

    ////////////////////////////////////////////////////////////////////
    // Enumerate my sites: local coords and global lexicographic indices
    ////////////////////////////////////////////////////////////////////
    myLcoor.resize(lsites);
    myGsite.resize(lsites);
    for(int ss=0; ss<lsites; ss++){
      Coordinate lcoor(nd);
      Lexicographic::CoorFromIndex(lcoor, ss, ldims);
      Coordinate gcoor(nd);
      for(int d=0; d<nd; d++) gcoor[d] = grid->_lstart[d] + lcoor[d];
      int64_t gsite;
      Lexicographic::IndexFromCoor(gcoor, gsite, gdims);
      myLcoor[ss] = lcoor;
      myGsite[ss] = gsite;
    }

    slab.resize((uint64_t)nrows * N);

    double t0 = usecond();
    ////////////////////////////////////////////////////////////////////
    // 0. Slab cache: SLAB_FILE=<stem> -> per-rank raw file <stem>.<rank>.
    //    Present: load, skip probe/gather/factor/solve (setup ~free on
    //    resubmits and sweeps).  Absent: full setup, then write it.
    //    VERIFY below runs in BOTH paths, certifying loaded data.
    ////////////////////////////////////////////////////////////////////
    bool loaded = false;
    char *sfile = getenv("SLAB_FILE");
    std::string slabfile;
    if (sfile) {
      slabfile = std::string(sfile) + "." + std::to_string(grid->ThisRank());
      FILE *f = fopen(slabfile.c_str(),"rb");
      if (f) {
        int64_t hdr[4] = {0,0,0,0};
        GRID_ASSERT( fread(hdr,sizeof(int64_t),4,f) == 4 );
        GRID_ASSERT( hdr[0] == (int64_t)0x44454E5345 );  // magic "DENSE"
        GRID_ASSERT( hdr[1] == N && hdr[2] == (int64_t)nrows && hdr[3] == (int64_t)nbasis );
        uint64_t nelem = (uint64_t)nrows * N;
        GRID_ASSERT( fread(&slab[0], sizeof(ComplexF), nelem, f) == nelem );
        fclose(f);
        loaded = true;
        std::cout << GridLogMessage << "DistributedDenseInverse: slab loaded from "
                  << slabfile << " -- skipping probe/factor/solve" << std::endl;
      } else {
        std::cout << GridLogMessage << "DistributedDenseInverse: slab cache " << slabfile
                  << " absent -- full setup, will write it" << std::endl;
      }
    }
    if (!loaded) {
    ////////////////////////////////////////////////////////////////////
    // 1. PROBE assembly of my rows of A: column (jsite,b) = Op(e_{jsite,b})
    ////////////////////////////////////////////////////////////////////
    Field e(grid);
    Field Ae(grid);
    Coordinate gcoorj(nd);
    int64_t gsitesN = grid->gSites();
    for(int64_t jsite=0; jsite<gsitesN; jsite++){
      Lexicographic::CoorFromIndex(gcoorj, jsite, gdims);
      for(int b=0; b<nbasis; b++){
        int64_t col = jsite*nbasis + b;
        sobj s = Zero();
        ((ComplexD *)&s)[b] = ComplexD(1.0,0.0);
        e = Zero();
        pokeSite(s, e, gcoorj);
        _Op.Op(e, Ae);
        for(int ss=0; ss<lsites; ss++){
          sobj t;
          peekLocalSite(t, Ae, myLcoor[ss]);
          for(int a=0; a<nbasis; a++){
            slab[ (uint64_t)(ss*nbasis+a)*N + col ] = ComplexF( ((ComplexD *)&t)[a] );
          }
        }
      }
      if ( (jsite % (gsitesN/20+1)) == 0 ) {
        std::cout << GridLogMessage << "DistributedDenseInverse: probed site "
                  << jsite << " / " << gsitesN << std::endl;
      }
    }
    double t1 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: probe assembly took "
              << (t1-t0)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // 2. Gather full A (fp32) by chunked zero-fill+GlobalSum, streamed
    //    straight into the boss-GCD device buffer (host never holds it).
    ////////////////////////////////////////////////////////////////////
    int boss = grid->IsBoss();
    std::vector<ComplexF> Afull;
#ifdef GRID_HIP
    rocblas_float_complex *dA = nullptr;
    rocblas_float_complex *dB = nullptr;   // getrs_64 RHS block (CHUNKROWS identity columns)
    rocblas_handle rochandle;
    int64_t *dIpiv = nullptr;              // ILP64 pivots, live from factor to last getrs
    uint64_t Abytes = (uint64_t)N * N * sizeof(ComplexF);
    // Make HBM space for the naked hipMalloc: flush the device-copy layer.
    // (FreePool of the allocator free-list awaits the type-dispatched fix.)
    MemoryManager::EvictAll();
    // MemoryManager::FreePool();
    if (boss) {
      auto aerr = hipMalloc((void **)&dA, Abytes);
      if (aerr != hipSuccess) {
        std::cout << GridLogMessage << "DistributedDenseInverse: hipMalloc of "
                  << Abytes/1024./1024./1024. << " GB FAILED -- reduce --device-mem "
                  << "or enable the fixed MemoryManager::FreePool()" << std::endl;
        GRID_ASSERT(aerr == hipSuccess);
      }
      std::cout << GridLogMessage << "DistributedDenseInverse: device inversion buffer allocated ("
                << Abytes/1024./1024./1024. << " GB)" << std::endl;
    }
#else
    if (boss) Afull.resize((uint64_t)N * N);
#endif
    {
      std::unordered_map<int64_t,int> rowmap;   // global row -> my slab row
      for(int ss=0; ss<lsites; ss++)
        for(int a=0; a<nbasis; a++)
          rowmap[ myGsite[ss]*nbasis + a ] = ss*nbasis + a;

      std::vector<ComplexF> chunk((uint64_t)CHUNKROWS * N);
      for(int64_t row0=0; row0<N; row0+=CHUNKROWS){
        int64_t nrow = std::min(CHUNKROWS, N-row0);
        uint64_t nelem = (uint64_t)nrow * N;
        for(uint64_t i=0;i<nelem;i++) chunk[i]=ComplexF(0.0,0.0);
        for(int64_t r=row0; r<row0+nrow; r++){
          auto it = rowmap.find(r);
          if (it != rowmap.end()) {
            uint64_t src = (uint64_t)(it->second) * N;
            uint64_t dst = (uint64_t)(r-row0) * N;
            for(int64_t j=0;j<N;j++) chunk[dst+j] = slab[src+j];
          }
        }
        grid->GlobalSumVector(&chunk[0], (int)nelem);
        if (boss) {
#ifdef GRID_HIP
          GRID_ASSERT( hipMemcpy((char *)dA + (uint64_t)row0*N*sizeof(ComplexF),
                                 &chunk[0], nelem*sizeof(ComplexF),
                                 hipMemcpyHostToDevice) == hipSuccess );
#else
          uint64_t dst = (uint64_t)row0 * N;
          for(uint64_t i=0;i<nelem;i++) Afull[dst+i] = chunk[i];
#endif
        }
      }
    }
    double t2 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: gather to boss took "
              << (t2-t1)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // 3. Factor in place on the boss (fp32, ILP64).
    //    Row-major buffer handed to column-major LAPACK => LU of A^T.
    ////////////////////////////////////////////////////////////////////
    if (boss) {
#ifdef GRID_HIP
      std::cout << GridLogMessage << "DistributedDenseInverse: rocSOLVER cgetrf_64 (ILP64 LU) N=" << N
                << " in place on resident device buffer" << std::endl;
      auto hst = rocblas_create_handle(&rochandle);
      std::cout << GridLogMessage << "DistributedDenseInverse: rocblas handle status " << (int)hst << std::endl;
      int64_t *dInfo;
      GRID_ASSERT( hipMalloc((void **)&dIpiv, N*sizeof(int64_t)) == hipSuccess );
      GRID_ASSERT( hipMalloc((void **)&dInfo, sizeof(int64_t))   == hipSuccess );
      auto st1 = rocsolver_cgetrf_64(rochandle, (int64_t)N, (int64_t)N, dA, (int64_t)N, dIpiv, dInfo);
      hipDeviceSynchronize();
      int64_t info_h = -1;
      hipMemcpy(&info_h, dInfo, sizeof(int64_t), hipMemcpyDeviceToHost);
      std::cout << GridLogMessage << "DistributedDenseInverse: cgetrf_64 status " << (int)st1
                << " info = " << (int)info_h << std::endl;
      GRID_ASSERT(st1 == rocblas_status_success);
      GRID_ASSERT(info_h == 0);
      hipFree(dInfo);
      GRID_ASSERT( hipMalloc((void **)&dB, (uint64_t)CHUNKROWS*N*sizeof(ComplexF)) == hipSuccess );
      // dA holds the LU of A^T; rows of A^{-1} are produced blockwise below via
      // cgetrs_64 on identity-column blocks: A^T X = E => X columns = rows of
      // A^{-1}, in exactly the linear layout the harvest expects.
#else
      // Eigen fallback: small local CPU tests only.
      std::cout << GridLogMessage << "DistributedDenseInverse: Eigen fallback inversion N=" << N
                << (N > 10000 ? "  (WARNING: SLOW; use the HIP/rocSOLVER path)" : "")
                << std::endl;
      typedef Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatF;
      Eigen::Map<MatF> A(reinterpret_cast<std::complex<float>*>(&Afull[0]), N, N);
      MatF Ainv = A.inverse();
      A = Ainv;
#endif
    }
    double t3 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: factorisation took "
              << (t3-t2)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // 4. Blocked solve + broadcast: rows of A^{-1} chunk by chunk; each
    //    rank keeps the rows of its own sites (ownership-aligned).
    ////////////////////////////////////////////////////////////////////
    {
      std::unordered_map<int64_t,int> rowmap;
      for(int ss=0; ss<lsites; ss++)
        for(int a=0; a<nbasis; a++)
          rowmap[ myGsite[ss]*nbasis + a ] = ss*nbasis + a;

      std::vector<ComplexF> chunk((uint64_t)CHUNKROWS * N);
      for(int64_t row0=0; row0<N; row0+=CHUNKROWS){
        int64_t nrow = std::min(CHUNKROWS, N-row0);
        uint64_t nelem = (uint64_t)nrow * N;
        if (boss) {
#ifdef GRID_HIP
          // Identity block E: column j = e_{row0+j}; solve A^T X = E so X's
          // columns are rows [row0,row0+nrow) of A^{-1}.
          for(uint64_t i=0;i<nelem;i++) chunk[i] = ComplexF(0.0,0.0);
          for(int64_t j=0;j<nrow;j++) chunk[(uint64_t)j*N + (uint64_t)(row0+j)] = ComplexF(1.0,0.0);
          GRID_ASSERT( hipMemcpy(dB, &chunk[0], nelem*sizeof(ComplexF), hipMemcpyHostToDevice) == hipSuccess );
          auto strs = rocsolver_cgetrs_64(rochandle, rocblas_operation_none,
                                          (int64_t)N, (int64_t)nrow,
                                          dA, (int64_t)N, dIpiv, dB, (int64_t)N);
          GRID_ASSERT(strs == rocblas_status_success);
          hipDeviceSynchronize();
          GRID_ASSERT( hipMemcpy(&chunk[0], dB, nelem*sizeof(ComplexF), hipMemcpyDeviceToHost) == hipSuccess );
#else
          uint64_t src = (uint64_t)row0 * N;
          for(uint64_t i=0;i<nelem;i++) chunk[i] = Afull[src+i];
#endif
        }
        grid->Broadcast(0, &chunk[0], nelem*sizeof(ComplexF));
        for(int64_t r=row0; r<row0+nrow; r++){
          auto it = rowmap.find(r);
          if (it != rowmap.end()) {
            uint64_t dst = (uint64_t)(it->second) * N;
            uint64_t src = (uint64_t)(r-row0) * N;
            for(int64_t j=0;j<N;j++) slab[dst+j] = chunk[src+j];
          }
        }
      }
    }
#ifdef GRID_HIP
    if (boss) {
      if (dA)    hipFree(dA);
      if (dB)    hipFree(dB);
      if (dIpiv) hipFree(dIpiv);
      rocblas_destroy_handle(rochandle);
    }
#endif
    double t4 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: blocked getrs solve+scatter took "
              << (t4-t3)/1.0e6 << " s" << std::endl;

    if (sfile) {
      FILE *f = fopen(slabfile.c_str(),"wb");
      GRID_ASSERT(f != nullptr);
      int64_t hdr[4] = { (int64_t)0x44454E5345, N, (int64_t)nrows, (int64_t)nbasis };
      GRID_ASSERT( fwrite(hdr,sizeof(int64_t),4,f) == 4 );
      uint64_t nelem = (uint64_t)nrows * N;
      GRID_ASSERT( fwrite(&slab[0], sizeof(ComplexF), nelem, f) == nelem );
      fclose(f);
      std::cout << GridLogMessage << "DistributedDenseInverse: slab written to " << slabfile << std::endl;
    }
    } // end !loaded: full setup path

#ifdef GRID_HIP
    ////////////////////////////////////////////////////////////////////
    // Device-resident slab on EVERY rank + per-rank rocblas handle:
    // the batched apply is then one cgemm against resident data.
    ////////////////////////////////////////////////////////////////////
    {
      uint64_t sbytes = (uint64_t)nrows * N * sizeof(ComplexF);
      GRID_ASSERT( hipMalloc((void **)&dSlab, sbytes) == hipSuccess );
      GRID_ASSERT( hipMemcpy(dSlab, &slab[0], sbytes, hipMemcpyHostToDevice) == hipSuccess );
      rocblas_create_handle(&applyHandle);
      GRID_ASSERT( hipHostMalloc((void **)&hXapply, (uint64_t)N    *MRHS_MAX*sizeof(ComplexF)) == hipSuccess );
      GRID_ASSERT( hipHostMalloc((void **)&hYapply, (uint64_t)nrows*MRHS_MAX*sizeof(ComplexF)) == hipSuccess );
      GRID_ASSERT( hipMalloc    ((void **)&dXapply, (uint64_t)N    *MRHS_MAX*sizeof(ComplexF)) == hipSuccess );
      GRID_ASSERT( hipMalloc    ((void **)&dYapply, (uint64_t)nrows*MRHS_MAX*sizeof(ComplexF)) == hipSuccess );
      devSum  = getenv("DENSE_DEVICE_SUM") ? 1 : 0;
      std::cout << GridLogMessage << "DistributedDenseInverse: slab resident on device ("
                << sbytes/1024./1024. << " MB/rank), persistent apply buffers allocated; "
                << (devSum ? "DEVICE-buffer allreduce (GPU-aware MPI)" : "host allreduce")
                << std::endl;
    }
#endif

    ////////////////////////////////////////////////////////////////////
    // 5. VERIFY: || A (Ainv x) - x || / ||x|| plus one timed apply.
    ////////////////////////////////////////////////////////////////////
    {
      Field x(grid); Field y(grid); Field z(grid);
      x = ComplexD(1.0,0.0);
      double ta = usecond();
      (*this)(x, y);
      double tb = usecond();
      _Op.Op(y, z);
      z = z - x;
      RealD rel = std::sqrt(norm2(z)/norm2(x));
      std::cout << GridLogMessage << "DistributedDenseInverse: VERIFY ||A Ainv x - x||/||x|| = "
                << rel << "   (one apply took " << (tb-ta)/1000.0 << " ms)" << std::endl;
      GRID_ASSERT(rel < 1.0e-2);
    }
    std::cout << GridLogMessage << "DistributedDenseInverse: setup complete, total "
              << (usecond()-t0)/1.0e6 << " s" << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // Single-RHS apply: ONE GlobalSum + slab GEMV + poke my own sites.
  ////////////////////////////////////////////////////////////////////
  virtual void operator()(const Field &src, Field &psi)
  {
    std::vector<ComplexD> x((uint64_t)N, ComplexD(0.0,0.0));
    for(int ss=0; ss<lsites; ss++){
      sobj s;
      peekLocalSite(s, src, myLcoor[ss]);
      for(int b=0; b<nbasis; b++)
        x[ myGsite[ss]*nbasis + b ] = ((ComplexD *)&s)[b];
    }
    grid->GlobalSumVector(&x[0], (int)N);

    std::vector<ComplexD> y(nrows);
    thread_for(r, nrows, {
      ComplexD acc(0.0,0.0);
      const ComplexF *row = &slab[(uint64_t)r * N];
      for(int64_t j=0; j<N; j++) acc += ComplexD(row[j]) * x[j];
      y[r] = acc;
    });

    for(int ss=0; ss<lsites; ss++){
      sobj s;
      for(int b=0; b<nbasis; b++)
        ((ComplexD *)&s)[b] = y[ss*nbasis + b];
      pokeLocalSite(s, psi, myLcoor[ss]);
    }

    if ( getenv("DENSE_CC_CHECK") ) {
      Field tmp(grid);
      _Op.Op(psi, tmp);
      tmp = tmp - src;
      std::cout << GridLogMessage << "DistributedDenseInverse: apply defect ||A x - b||/||b|| = "
                << std::sqrt(norm2(tmp)/norm2(src)) << std::endl;
    }
  }

  ////////////////////////////////////////////////////////////////////
  // Batched (mrhs) apply: ONE GlobalSum of all nrhs packed vectors,
  // ONE slab pass computing every rhs (GEMM shape) -- per-rhs cost
  // below the single-RHS apply.  nrhs <= 32.
  ////////////////////////////////////////////////////////////////////
  void ApplyBatch(std::vector<Field> &src, std::vector<Field> &psi)
  {
    int nr = src.size();
    GRID_ASSERT(nr <= 32);
    double t0 = usecond();

    std::vector<ComplexD> X((uint64_t)N*nr, ComplexD(0.0,0.0));   // column r = packed src[r]
    for(int rr=0; rr<nr; rr++){
      for(int ss=0; ss<lsites; ss++){
        sobj s;
        peekLocalSite(s, src[rr], myLcoor[ss]);
        for(int b=0; b<nbasis; b++)
          X[ (uint64_t)rr*N + myGsite[ss]*nbasis + b ] = ((ComplexD *)&s)[b];
      }
    }
    grid->GlobalSumVector(&X[0], (int)((uint64_t)N*nr));

#ifdef GRID_HIP
    ////////////////////////////////////////////////////////////////////
    // Device GEMM path: Y(nrows x nr) = A^{-1}_rows x X.
    // The slab buffer (row-major nrows x N) IS col-major A^T (N x nrows,
    // lda=N), so Y = op(slab,T) * X in rocblas col-major convention.
    ////////////////////////////////////////////////////////////////////
    std::vector<ComplexF> Xf((uint64_t)N*nr);
    for(uint64_t i=0;i<(uint64_t)N*nr;i++) Xf[i]=ComplexF(X[i]);
    rocblas_float_complex *dX; rocblas_float_complex *dY;
    GRID_ASSERT( hipMalloc((void **)&dX, (uint64_t)N*nr*sizeof(ComplexF))     == hipSuccess );
    GRID_ASSERT( hipMalloc((void **)&dY, (uint64_t)nrows*nr*sizeof(ComplexF)) == hipSuccess );
    GRID_ASSERT( hipMemcpy(dX, &Xf[0], (uint64_t)N*nr*sizeof(ComplexF), hipMemcpyHostToDevice) == hipSuccess );
    rocblas_float_complex one (1.0f,0.0f);
    rocblas_float_complex zero(0.0f,0.0f);
    auto gst = rocblas_cgemm(applyHandle,
                             rocblas_operation_transpose, rocblas_operation_none,
                             (rocblas_int)nrows, (rocblas_int)nr, (rocblas_int)N,
                             &one,  dSlab, (rocblas_int)N,
                                    dX,    (rocblas_int)N,
                             &zero, dY,    (rocblas_int)nrows);
    GRID_ASSERT(gst == rocblas_status_success);
    hipDeviceSynchronize();
    std::vector<ComplexF> Yf((uint64_t)nrows*nr);
    GRID_ASSERT( hipMemcpy(&Yf[0], dY, (uint64_t)nrows*nr*sizeof(ComplexF), hipMemcpyDeviceToHost) == hipSuccess );
    hipFree(dX); hipFree(dY);

    for(int rr=0; rr<nr; rr++){
      for(int ss=0; ss<lsites; ss++){
        sobj s;
        for(int b=0; b<nbasis; b++)
          ((ComplexD *)&s)[b] = ComplexD(Yf[(uint64_t)rr*nrows + (ss*nbasis+b)]);   // Y col-major
        pokeLocalSite(s, psi[rr], myLcoor[ss]);
      }
    }
#else
    std::vector<ComplexD> Y((uint64_t)nrows*nr);
    thread_for(r_, nrows, {
      ComplexD acc[32];
      for(int rr=0; rr<nr; rr++) acc[rr]=ComplexD(0.0,0.0);
      const ComplexF *row = &slab[(uint64_t)r_ * N];
      for(int64_t j=0; j<N; j++){
        ComplexD s(row[j]);
        for(int rr=0; rr<nr; rr++) acc[rr] += s * X[(uint64_t)rr*N + j];
      }
      for(int rr=0; rr<nr; rr++) Y[(uint64_t)r_*nr + rr] = acc[rr];
    });

    for(int rr=0; rr<nr; rr++){
      for(int ss=0; ss<lsites; ss++){
        sobj s;
        for(int b=0; b<nbasis; b++)
          ((ComplexD *)&s)[b] = Y[(uint64_t)(ss*nbasis+b)*nr + rr];
        pokeLocalSite(s, psi[rr], myLcoor[ss]);
      }
    }
#endif
    double t1 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: batched apply " << nr << " rhs took "
              << (t1-t0)/1000.0 << " ms  (" << (t1-t0)/1000.0/nr << " ms/rhs)" << std::endl;

    if ( getenv("DENSE_CC_CHECK") ) {
      Field tmp(grid);
      for(int rr=0; rr<nr; rr++){
        _Op.Op(psi[rr], tmp);
        tmp = tmp - src[rr];
        std::cout << GridLogMessage << "DistributedDenseInverse: batch defect["<<rr<<"] = "
                  << std::sqrt(norm2(tmp)/norm2(src[rr])) << std::endl;
      }
    }
  }

  ////////////////////////////////////////////////////////////////////
  // 6D mrhs apply: operates DIRECTLY on the packed 6D field (rhs =
  // dim 0), replacing the per-rhs slice/split layer of the wrapper.
  //   - ONE CPU view open per field, not one per site per rhs
  //   - fp32 end-to-end; the allreduce is EXACT despite fp32 because
  //     zero-fill assembly gives every element exactly one nonzero
  //     contributing rank (sum of zeros + v = v in any precision)
  //   - persistent pinned + device buffers: no per-call hipMalloc,
  //     no page faults, full PCIe rate on the H2D
  // Per-stage timers print under DENSE_APPLY_PROFILE.
  ////////////////////////////////////////////////////////////////////
  void ApplyBatch6D(const Field &in6, Field &out6, int nr)
  {
    GRID_ASSERT(nr <= MRHS_MAX);
    GRID_ASSERT(in6.Grid()->_ndimension == nd+1);   // {rhs, s, x,y,z,t}
    double t0 = usecond();
    Field &in = const_cast<Field &>(in6);
    uint64_t nX = (uint64_t)N * nr;
#ifdef GRID_HIP
    memset(hXapply, 0, nX*sizeof(ComplexF));
    {
      autoView(iv, in, CpuRead);
      Coordinate c6(nd+1);
      for(int ss=0; ss<lsites; ss++){
        for(int d=0; d<nd; d++) c6[d+1] = myLcoor[ss][d];
        for(int rr=0; rr<nr; rr++){
          c6[0] = rr;
          sobj s;
          peekLocalSite(s, iv, c6);
          for(int b=0; b<nbasis; b++)
            hXapply[(uint64_t)rr*N + myGsite[ss]*nbasis + b] = ComplexF(((ComplexD *)&s)[b]);
        }
      }
    }
    double t1 = usecond();
    double t2, t3;
    if (devSum) {
      // H2D first, then allreduce the DEVICE buffer (GPU-aware MPI, NIC-direct).
      // fp32 sum remains exact: zero-fill = one contributing rank per element.
      GRID_ASSERT( hipMemcpy(dXapply, hXapply, nX*sizeof(ComplexF), hipMemcpyHostToDevice) == hipSuccess );
      t2 = usecond();
      grid->GlobalSumVector((ComplexF *)dXapply, (int)nX);
      t3 = usecond();
    } else {
      grid->GlobalSumVector(hXapply, (int)nX);   // fp32 allreduce: exact for zero-fill assembly
      t2 = usecond();
      GRID_ASSERT( hipMemcpy(dXapply, hXapply, nX*sizeof(ComplexF), hipMemcpyHostToDevice) == hipSuccess );
      t3 = usecond();
    }
    rocblas_float_complex one (1.0f,0.0f);
    rocblas_float_complex zero(0.0f,0.0f);
    // Y = op(slab,T) . X : row-major slab (nrows x N) == col-major A^T (N x nrows,
    // lda=N), so transpose gives the nrows x N operator.  Tiny-output/huge-K =>
    // GPU under-filled (~13ms); the software-split-K fix is a TODO via
    // GridBLAS.gemmBatched (see class-header note).
    auto gst = rocblas_cgemm(applyHandle,
                             rocblas_operation_transpose, rocblas_operation_none,
                             (rocblas_int)nrows, (rocblas_int)nr, (rocblas_int)N,
                             &one,  dSlab,   (rocblas_int)N,
                                    dXapply, (rocblas_int)N,
                             &zero, dYapply, (rocblas_int)nrows);
    GRID_ASSERT(gst == rocblas_status_success);
    hipDeviceSynchronize();
    double t4 = usecond();
    GRID_ASSERT( hipMemcpy(hYapply, dYapply, (uint64_t)nrows*nr*sizeof(ComplexF), hipMemcpyDeviceToHost) == hipSuccess );
    double t5 = usecond();
    {
      autoView(ov, out6, CpuWrite);
      Coordinate c6(nd+1);
      for(int ss=0; ss<lsites; ss++){
        for(int d=0; d<nd; d++) c6[d+1] = myLcoor[ss][d];
        for(int rr=0; rr<nr; rr++){
          c6[0] = rr;
          sobj s;
          for(int b=0; b<nbasis; b++)
            ((ComplexD *)&s)[b] = ComplexD(hYapply[(uint64_t)rr*nrows + (ss*nbasis+b)]);  // Y col-major
          pokeLocalSite(s, ov, c6);
        }
      }
    }
    double t6 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: apply6D " << nr << " rhs took "
              << (t6-t0)/1000.0 << " ms" << std::endl;
    if ( getenv("DENSE_APPLY_PROFILE") ) {
      double t_sum = devSum ? (t3-t2) : (t2-t1);   // devSum runs H2D first
      double t_h2d = devSum ? (t2-t1) : (t3-t2);
      std::cout << GridLogMessage << "DistributedDenseInverse: apply6D profile:"
                << " pack "      << (t1-t0)/1000.0
                << "  allreduce " << t_sum/1000.0
                << "  H2D "       << t_h2d/1000.0
                << "  gemm "      << (t4-t3)/1000.0
                << "  D2H "       << (t5-t4)/1000.0
                << "  unpack "    << (t6-t5)/1000.0 << "  ms" << std::endl;
    }
#else
    // Host fallback (laptop loop-closer): same single-view pack, GEMV slab pass.
    std::vector<ComplexD> X(nX, ComplexD(0.0,0.0));
    {
      autoView(iv, in, CpuRead);
      Coordinate c6(nd+1);
      for(int ss=0; ss<lsites; ss++){
        for(int d=0; d<nd; d++) c6[d+1] = myLcoor[ss][d];
        for(int rr=0; rr<nr; rr++){
          c6[0] = rr;
          sobj s;
          peekLocalSite(s, iv, c6);
          for(int b=0; b<nbasis; b++)
            X[(uint64_t)rr*N + myGsite[ss]*nbasis + b] = ((ComplexD *)&s)[b];
        }
      }
    }
    grid->GlobalSumVector(&X[0], (int)nX);
    std::vector<ComplexD> Y((uint64_t)nrows*nr);
    thread_for(r_, nrows, {
      ComplexD acc[MRHS_MAX];
      for(int rr=0; rr<nr; rr++) acc[rr]=ComplexD(0.0,0.0);
      const ComplexF *row = &slab[(uint64_t)r_ * N];
      for(int64_t j=0; j<N; j++){
        ComplexD sD(row[j]);
        for(int rr=0; rr<nr; rr++) acc[rr] += sD * X[(uint64_t)rr*N + j];
      }
      for(int rr=0; rr<nr; rr++) Y[(uint64_t)rr*nrows + r_] = acc[rr];
    });
    {
      autoView(ov, out6, CpuWrite);
      Coordinate c6(nd+1);
      for(int ss=0; ss<lsites; ss++){
        for(int d=0; d<nd; d++) c6[d+1] = myLcoor[ss][d];
        for(int rr=0; rr<nr; rr++){
          c6[0] = rr;
          sobj s;
          for(int b=0; b<nbasis; b++)
            ((ComplexD *)&s)[b] = Y[(uint64_t)rr*nrows + (ss*nbasis+b)];
          pokeLocalSite(s, ov, c6);
        }
      }
    }
    double t6 = usecond();
    std::cout << GridLogMessage << "DistributedDenseInverse: apply6D " << nr << " rhs took "
              << (t6-t0)/1000.0 << " ms  (host GEMV path)" << std::endl;
#endif
  }
};

//////////////////////////////////////////////////////////////////////
// Dense CC solve on the PACKED 6D mrhs coarse-coarse field: unpack per
// rhs, ONE batched dense apply, repack.  Drop-in for the L3 PGCR.
//////////////////////////////////////////////////////////////////////
template<class CoarseCoarseField>
class MrhsDenseCCSolve : public LinearFunction<CoarseCoarseField> {
public:
  DistributedDenseInverse<CoarseCoarseField> &_Dense;
  GridBase *_CoarseCoarse5d;
  int _nrhs;
  MrhsDenseCCSolve(DistributedDenseInverse<CoarseCoarseField> &D, GridBase *cc5d, int nrhs)
    : _Dense(D), _CoarseCoarse5d(cc5d), _nrhs(nrhs) {}
  using LinearFunction<CoarseCoarseField>::operator();
  virtual void operator()(const CoarseCoarseField &in, CoarseCoarseField &out){
    if ( getenv("DENSE_CC_CHECK") ) {
      // Audit path: per-rhs 5D unpack so ApplyBatch can run the _Op defect
      // check per rhs.  ~50ms/call of slice/split overhead -- audit only.
      CoarseCoarseField tmp(in.Grid());
      tmp = in;
      std::vector<CoarseCoarseField> split_in (_nrhs,_CoarseCoarse5d);
      std::vector<CoarseCoarseField> split_out(_nrhs,_CoarseCoarse5d);
      for(int r=0;r<_nrhs;r++) ExtractSliceFast(split_in[r], tmp, r, 0);
      _Dense.ApplyBatch(split_in, split_out);
      for(int r=0;r<_nrhs;r++) InsertSliceFast(split_out[r], out, r, 0);
    } else {
      _Dense.ApplyBatch6D(in, out, _nrhs);
    }
  }
};

//////////////////////////////////////////////////////////////////////
// mrhs interfaces + single-polynomial mrhs PGCR (verbatim from Example_pvdagm_mrhs.cc)
//////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsLinearFunction {
public:
  virtual void operator()(std::vector<Field> &in, std::vector<Field> &out) = 0;
};

template<class Field>
class MrhsPGCRNonHermitian {
public:
  RealD Tolerance; Integer MaxIterations; int mmax,nstep,steps,level;
  int ZeroGuess = 0; int FirstCycle = 0;   // caller contract: zero guess => first-cycle r0 = src
  std::string name = "Level 1";
  LinearOperatorBase<Field> &Linop;
  MrhsLinearFunction<Field> &Preconditioner;
  void Level(int lv){ name = "Level " + std::to_string(lv); level=lv; }
  void Name(std::string n){ name = n; }
  void SetZeroGuess(int z){ ZeroGuess=z; }
  MrhsPGCRNonHermitian(RealD tol,Integer maxit,LinearOperatorBase<Field> &_Linop,MrhsLinearFunction<Field> &Prec,int _mmax,int _nstep)
    : Tolerance(tol),MaxIterations(maxit),Linop(_Linop),Preconditioner(Prec),mmax(_mmax),nstep(_nstep){ level=1; }
  static RealD vnorm2(std::vector<Field> &x){ RealD s=0; for(auto &f:x) s+=norm2(f); return s; }
  static ComplexD vinnerProduct(std::vector<Field> &x,std::vector<Field> &y){ ComplexD s(0); for(int r=0;r<(int)x.size();r++) s+=innerProduct(x[r],y[r]); return s; }
  static void vaxpy(std::vector<Field> &z,ComplexD a,std::vector<Field> &x,std::vector<Field> &y){ for(int r=0;r<(int)z.size();r++) axpy(z[r],a,x[r],y[r]); }
  void vOp(std::vector<Field> &in,std::vector<Field> &out){ for(int r=0;r<(int)in.size();r++) Linop.Op(in[r],out[r]); }
  void operator()(std::vector<Field> &src,std::vector<Field> &psi){
    RealD cp,ssq,rsq; int nrhs=src.size(); GridBase *grid=src[0].Grid();
    ssq=vnorm2(src); rsq=Tolerance*Tolerance*ssq;
    std::vector<Field> r(nrhs,grid);
    GridStopWatch T; T.Start(); steps=0; FirstCycle=1;
    for(int k=0;k<MaxIterations;k++){
      cp=GCRnStep(src,psi,rsq);
      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR("<<mmax<<","<<nstep<<") "<<steps<<" steps cp = "<<cp<<" target "<<rsq<<std::endl;
      if(cp<rsq){
        T.Stop(); vOp(psi,r); for(int rr=0;rr<nrhs;rr++) axpy(r[rr],-1.0,src[rr],r[rr]);
        RealD tr=vnorm2(r);
        std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR: Converged on iteration "<<steps
                 <<" computed residual "<<std::sqrt(cp/ssq)<<" true residual "<<std::sqrt(tr/ssq)<<" target "<<Tolerance<<std::endl;
        std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR Time elapsed: Total "<<T.Elapsed()<<std::endl;
        for(int rr=0;rr<nrhs;rr++){ RealD rn=std::sqrt(norm2(r[rr])/norm2(src[rr])); std::cout<<GridLogMessage<<"MrhsPGCR per-rhs true residual["<<rr<<"] = "<<rn<<std::endl; }
        return;
      }
    }
    std::cout<<GridLogMessage<<"MrhsPGCR: did not converge"<<std::endl;
  }
  RealD GCRnStep(std::vector<Field> &src,std::vector<Field> &psi,RealD rsq){
    RealD cp; ComplexD a,b,rq; RealD zAAz; int nrhs=src.size(); GridBase *grid=src[0].Grid();
    std::vector<Field> r(nrhs,grid),z(nrhs,grid),Az(nrhs,grid);
    std::vector< std::vector<Field> > q(mmax,std::vector<Field>(nrhs,grid));
    std::vector< std::vector<Field> > p(mmax,std::vector<Field>(nrhs,grid));
    std::vector<RealD> qq(mmax);
    std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR nStep("<<nstep<<")"<<std::endl;
    if (ZeroGuess && FirstCycle) { for(int rr=0;rr<nrhs;rr++){ psi[rr]=Zero(); r[rr]=src[rr]; } }
    else                         { vOp(psi,Az); for(int rr=0;rr<nrhs;rr++) r[rr]=src[rr]-Az[rr]; }
    FirstCycle=0;
    Preconditioner(r,z); vOp(z,Az); zAAz=vnorm2(Az);
    p[0]=z; q[0]=Az; qq[0]=zAAz; cp=vnorm2(r);
    for(int k=0;k<nstep;k++){
      steps++; int kp=k+1, peri_k=k%mmax, peri_kp=kp%mmax;
      rq=vinnerProduct(q[peri_k],r); a=rq/qq[peri_k];
      vaxpy(psi,a,p[peri_k],psi); vaxpy(r,-a,q[peri_k],r); cp=vnorm2(r);
      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR step["<<steps<<"]  resid "<<cp<<" target "<<rsq<<std::endl;
      if((k==nstep-1)||(cp<rsq)) return cp;
      Preconditioner(r,z); vOp(z,Az); zAAz=vnorm2(Az);
      q[peri_kp]=Az; p[peri_kp]=z;
      int northog=((kp)>(mmax-1))?(mmax-1):(kp);
      for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; GRID_ASSERT((k-back)>=0);
        b=-real(vinnerProduct(q[peri_back],Az))/qq[peri_back];
        vaxpy(p[peri_kp],b,p[peri_back],p[peri_kp]); vaxpy(q[peri_kp],b,q[peri_back],q[peri_kp]); }
      qq[peri_kp]=vnorm2(q[peri_kp]);
    }
    GRID_ASSERT(0); return cp;
  }
};

//////////////////////////////////////////////////////////////////////
// L2->L3 mrhs V-cycle: LinearFunction on the 6D mrhs COARSE field.
// The coarse-coarse solve slot now takes EITHER the dense mrhs solve
// (DENSE_CC=1) or the L3 PGCR (DENSE_CC=0).
//////////////////////////////////////////////////////////////////////
template<class CoarseField, class CoarseCoarseField>
class MrhsCoarseThreeLevelPrec : public LinearFunction<CoarseField> {
public:
  LinearOperatorBase<CoarseField>          &_CoarseOp;          // mrhs coarse op (6D)
  LinearFunction<CoarseField>              &_CoarseSmoother;    // shifted 6D coarse smoother
  MultiRHSBlockProject<CoarseField>        &_Projector;         // L2->L3 (vector-based)
  LinearFunction<CoarseCoarseField>        &_CoarseCoarseSolve; // L3 solve (6D cc)
  GridBase *_Coarse5d, *_CoarseCoarse5d, *_CoarseCoarseMrhs;
  int _nrhs;

  MrhsCoarseThreeLevelPrec(LinearOperatorBase<CoarseField> &CoarseOp,
                           LinearFunction<CoarseField> &CoarseSmoother,
                           MultiRHSBlockProject<CoarseField> &Projector,
                           LinearFunction<CoarseCoarseField> &CoarseCoarseSolve,
                           GridBase *Coarse5d, GridBase *CoarseCoarse5d, GridBase *CoarseCoarseMrhs, int nrhs)
    : _CoarseOp(CoarseOp), _CoarseSmoother(CoarseSmoother), _Projector(Projector),
      _CoarseCoarseSolve(CoarseCoarseSolve),
      _Coarse5d(Coarse5d), _CoarseCoarse5d(CoarseCoarse5d), _CoarseCoarseMrhs(CoarseCoarseMrhs), _nrhs(nrhs) {}

  using LinearFunction<CoarseField>::operator();
  virtual void operator()(const CoarseField &in, CoarseField &out) {
    int nrhs=_nrhs; double t;
    CoarseField vec1(in.Grid());
    CoarseField vec2(in.Grid());

    // trivial pre-smoother
    out = in;

    // residual (6D coarse)
    _CoarseOp.Op(out,vec1);  sub(vec1,in,vec1);

    // restrict: unpack 6D coarse -> vector<Coarse> -> blockProject -> vector<CoarseCoarse> -> pack 6D cc
    std::vector<CoarseField>       csplit(nrhs,_Coarse5d);
    std::vector<CoarseCoarseField> ccsplit(nrhs,_CoarseCoarse5d);
    CoarseCoarseField CCsrc(_CoarseCoarseMrhs);
    CoarseCoarseField CCsol(_CoarseCoarseMrhs);

    t=-usecond();
    for(int r=0;r<nrhs;r++) ExtractSliceFast(csplit[r],vec1,r,0);
    _Projector.blockProject(csplit,ccsplit);
    for(int r=0;r<nrhs;r++) InsertSliceFast(ccsplit[r],CCsrc,r,0);
    t+=usecond();
    std::cout<<GridLogMessage<<"L2->L3 restrict took "<<t/1000.0<<"ms"<<std::endl;

    // L3 solve (dense mrhs GEMM, or PGCR)
    t=-usecond();
    CCsol=Zero();
    _CoarseCoarseSolve(CCsrc,CCsol);
    t+=usecond();
    std::cout<<GridLogMessage<<"L3 coarse-coarse solve took "<<t/1000.0<<"ms  ("<<t/1000.0/nrhs<<"ms/rhs)"<<std::endl;

    // prolong: unpack 6D cc -> blockPromote -> pack 6D coarse; add correction
    t=-usecond();
    for(int r=0;r<nrhs;r++) ExtractSliceFast(ccsplit[r],CCsol,r,0);
    _Projector.blockPromote(csplit,ccsplit);
    for(int r=0;r<nrhs;r++) InsertSliceFast(csplit[r],vec1,r,0);
    add(out,out,vec1);
    t+=usecond();
    std::cout<<GridLogMessage<<"L2->L3 prolong took "<<t/1000.0<<"ms"<<std::endl;

    // residual + coarse smoother (6D coarse)
    _CoarseOp.Op(out,vec1);  sub(vec1,in,vec1);
    vec2=Zero();
    _CoarseSmoother(vec1,vec2);
    add(out,out,vec2);
  }
};

//////////////////////////////////////////////////////////////////////
// L1->L2 mrhs V-cycle (verbatim from Example_pvdagm_mrhs.cc)
//////////////////////////////////////////////////////////////////////
template<class FineField, class MrhsCoarseVector, class FineSmoother>
class MrhsTwoLevelMG : public MrhsLinearFunction<FineField> {
public:
  typedef MrhsCoarseVector CoarseVector;
  LinearOperatorBase<FineField>   &_FineOperator;
  FineSmoother                    &_PostSmoother;
  MultiRHSBlockProject<FineField> &_Projector;
  LinearFunction<CoarseVector>    &_CoarseSolve;
  GridBase *_CoarseGrid, *_CoarseGridMrhs;
  MrhsTwoLevelMG(LinearOperatorBase<FineField> &FineOp, FineSmoother &Post,
                 MultiRHSBlockProject<FineField> &Projector, LinearFunction<CoarseVector> &CoarseSolve,
                 GridBase *CoarseGrid, GridBase *CoarseGridMrhs)
    : _FineOperator(FineOp),_PostSmoother(Post),_Projector(Projector),_CoarseSolve(CoarseSolve),
      _CoarseGrid(CoarseGrid),_CoarseGridMrhs(CoarseGridMrhs){}
  virtual void operator()(std::vector<FineField> &in, std::vector<FineField> &out){
    int nrhs=in.size(); GridBase *fgrid=in[0].Grid(); double t;
    std::vector<FineField> vec1(nrhs,fgrid),vec2(nrhs,fgrid);
    for(int r=0;r<nrhs;r++) out[r]=in[r];
    for(int r=0;r<nrhs;r++){ _FineOperator.Op(out[r],vec1[r]); sub(vec1[r],in[r],vec1[r]); }
    std::vector<CoarseVector> Csrc_split(nrhs,_CoarseGrid), Csol_split(nrhs,_CoarseGrid);
    CoarseVector CsrcMrhs(_CoarseGridMrhs), CsolMrhs(_CoarseGridMrhs);
    t=-usecond();
    _Projector.blockProject(vec1,Csrc_split);
    for(int r=0;r<nrhs;r++) InsertSliceFast(Csrc_split[r],CsrcMrhs,r,0);
    t+=usecond(); std::cout<<GridLogMessage<<"Mrhs project+pack took "<<t/1000.0<<"ms"<<std::endl;
    t=-usecond(); CsolMrhs=Zero(); _CoarseSolve(CsrcMrhs,CsolMrhs); t+=usecond();
    std::cout<<GridLogMessage<<"Mrhs coarse solve took "<<t/1000.0<<"ms  ("<<t/1000.0/nrhs<<"ms/rhs)"<<std::endl;
    t=-usecond();
    for(int r=0;r<nrhs;r++) ExtractSliceFast(Csol_split[r],CsolMrhs,r,0);
    _Projector.blockPromote(vec1,Csol_split);
    for(int r=0;r<nrhs;r++) add(out[r],out[r],vec1[r]);
    t+=usecond(); std::cout<<GridLogMessage<<"Mrhs unpack+promote took "<<t/1000.0<<"ms"<<std::endl;
    for(int r=0;r<nrhs;r++){ _FineOperator.Op(out[r],vec1[r]); sub(vec1[r],in[r],vec1[r]); }
    t=-usecond();
    for(int r=0;r<nrhs;r++){ vec2[r]=Zero(); _PostSmoother(vec1[r],vec2[r]); add(out[r],out[r],vec2[r]); }
    t+=usecond(); std::cout<<GridLogMessage<<"Mrhs post-smooth took "<<t/1000.0<<"ms  ("<<t/1000.0/nrhs<<"ms/rhs)"<<std::endl;
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  ParseEnvironment();

  const int Ls=24; RealD M5=1.8, b=1.5, c=0.5;
  //  const int nbasis=56; const int nrhs=Nrhs;
  const int nbasis=60; const int nrhs=Nrhs;
  GRID_ASSERT(nrhs % vComplex::Nsimd() == 0);

  std::vector<int> lat_size {48,48,48,96};
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Level 1 blocking (default 2^4)
  Coordinate clatt = lat_size;
  Coordinate Block({2,2,2,2});
  if ( getenv("BLOCK") ){ GridCmdOptionIntVector(std::string(getenv("BLOCK")),Block); GRID_ASSERT(Block.size()==4); }
  for(int d=0;d<4;d++){ GRID_ASSERT(lat_size[d]%Block[d]==0); clatt[d]=lat_size[d]/Block[d]; }
  std::cout << GridLogMessage << "Block  " << Block  << "  coarse lattice       " << clatt << std::endl;

  // Level 2 blocking: SUPERCOARSE default 8,4,3,6 -> CC [3,6,8,8], the dense floor.
  Coordinate cclatt = clatt;
  Coordinate Block2({8,4,3,6});
  if ( getenv("BLOCK2") ){ GridCmdOptionIntVector(std::string(getenv("BLOCK2")),Block2); GRID_ASSERT(Block2.size()==4); }
  for(int d=0;d<4;d++){ GRID_ASSERT(clatt[d]%Block2[d]==0); cclatt[d]=clatt[d]/Block2[d]; }
  std::cout << GridLogMessage << "Block2 " << Block2 << "  coarse-coarse lattice " << cclatt << std::endl;

  GridCartesian *Coarse4d = SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *Coarse5d = SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);
  GridCartesian *CoarseCoarse4d = SpaceTimeGrid::makeFourDimGrid(cclatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *CoarseCoarse5d = SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarse4d);

  // 6D mrhs grids: rhs is dim 0, SIMD across rhs
  Coordinate mpi=GridDefaultMpi();
  Coordinate rhMpi ({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  Coordinate rhSimd({vComplex::Nsimd(),1,1,1,1,1});
  Coordinate rhLatt ({nrhs,1,clatt[0], clatt[1], clatt[2], clatt[3]});
  Coordinate rhLatt2({nrhs,1,cclatt[0],cclatt[1],cclatt[2],cclatt[3]});
  GridCartesian *CoarseMrhs       = new GridCartesian(rhLatt, rhSimd,rhMpi);
  GridCartesian *CoarseCoarseMrhs = new GridCartesian(rhLatt2,rhSimd,rhMpi);

  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers({5,6,7,8});

  LatticeGaugeField Umu(UGrid);
  std::cout << GridLogMessage << "Reading gauge field" << std::endl;
  FieldMetaData header;
  std::string file("/ccs/home/poare/ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0, M5,b,c);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD>        PVdagM_t;
  typedef ShiftedPVdagMLinearOperator<MobiusFermionD,LatticeFermionD> ShiftedPVdagM_t;

  // Level 1 tensor types
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>       LittleDiracOperator;
  typedef MultiGeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>  MrhsLittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector                                CoarseVector;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>                  Subspace;

  // Level 2 tensor types (coarsening deepens the nest by one iScalar)
  typedef CoarseVector::vector_object                                      CoarseSiteObj;
  typedef iScalar<vTComplex>                                              vTTComplex;
  typedef GeneralCoarsenedMatrix<CoarseSiteObj,vTTComplex,nbasis>          LittleDiracOperatorL2;
  typedef MultiGeneralCoarsenedMatrix<CoarseSiteObj,vTTComplex,nbasis>     MrhsLittleDiracOperatorL2;
  typedef LittleDiracOperatorL2::CoarseVector                              CoarseCoarseVector;
  typedef Aggregation<CoarseSiteObj,vTTComplex,nbasis>                     SubspaceL2;

  PVdagM_t        PVdagM(Ddwf,Dpv);
  ShiftedPVdagM_t ShiftedPVdagM(FineSmootherShift,Ddwf,Dpv);

  NextToNearestStencilGeometry5D geom (Coarse5d);
  NextToNearestStencilGeometry5D geom2(CoarseCoarse5d);

  //////////////////////////////////////////////////////////////////////
  // Subspace: load RAW (no Orthogonalise!), or generate.
  //////////////////////////////////////////////////////////////////////
  std::string subspace_file = "/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/subspace_nb"
                            + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));
  uint64_t file_exists=0;
  if ( UGrid->IsBoss() ){ std::ifstream f(subspace_file); file_exists=f.good()?1:0; }
  UGrid->GlobalSum(file_exists);

  const int cb=0;
  Subspace AggregatesGCR(Coarse5d,FGrid,cb);
  if ( file_exists ){
    std::cout << GridLogMessage << "*** Loading subspace from disk (kept RAW) ***" << std::endl;
    loadSubspace(AggregatesGCR.subspace, subspace_file);
  } else {
    std::cout << GridLogMessage << "*** GCR subspace generation ***" << std::endl;
    AggregatesGCR.CreateSubspaceGCR(RNG5,PVdagM,nbasis);
    saveSubspace(AggregatesGCR.subspace, subspace_file);
  }

  // RAW copy of the fine null vectors BEFORE CoarsenOperator block-orthonormalises in place.
  std::vector<LatticeFermionD> rawNull(nbasis,FGrid);
  for(int k=0;k<nbasis;k++) rawNull[k]=AggregatesGCR.subspace[k];

  //////////////////////////////////////////////////////////////////////
  // Coarsen L1->L2 and L2->L3 with SINGLE-RHS machinery; import to mrhs via
  // CopyMatrix.  The L2 (coarse-coarse) single-RHS operator is HOISTED to
  // main scope: it is tiny at the supercoarse blocking and the dense setup
  // (probe/VERIFY/defect checks) needs it alive for the whole run.
  //////////////////////////////////////////////////////////////////////
  MrhsLittleDiracOperator   mrhsLittleDiracOpPV(geom,  CoarseMrhs);
  MrhsLittleDiracOperatorL2 mrhsLittleDiracOpL2(geom2, CoarseCoarseMrhs);
  MultiRHSBlockProject<LatticeFermionD> MrhsProjector;
  MultiRHSBlockProject<CoarseVector>    MrhsProjectorL2;

  LittleDiracOperatorL2 LittleDiracOpL2(geom2,Coarse5d,CoarseCoarse5d);
  NonHermitianLinearOperator<LittleDiracOperatorL2,CoarseCoarseVector> LinOpCC5d(LittleDiracOpL2);

  {
    // --- L1->L2 single-RHS coarse operator (scoped: its padded _A is the memory peak) ---
    LittleDiracOperator LittleDiracOpPV(geom,FGrid,Coarse5d);
    LittleDiracOpPV.CoarsenOperator(PVdagM, AggregatesGCR);   // orthonormalises AggregatesGCR.subspace in place
    mrhsLittleDiracOpPV.CopyMatrix(LittleDiracOpPV);
    MrhsProjector.Allocate(nbasis,FGrid,Coarse5d);
    MrhsProjector.ImportBasis(AggregatesGCR.subspace);
    NonHermitianLinearOperator<LittleDiracOperator,CoarseVector> LinOpCoarse(LittleDiracOpPV);

    // --- psi_coarse = P^dag (RAW fine null) -> Galerkin images, NOT e_k ---
    std::vector<CoarseVector> psi_coarse(nbasis,Coarse5d);
    for(int k=0;k<nbasis;k++) AggregatesGCR.ProjectToSubspace(psi_coarse[k], rawNull[k]);
    rawNull.clear(); rawNull.shrink_to_fit();
    {
      RealD s2=0.0;
      for(int i=0;i<nbasis;i++) for(int j=0;j<nbasis;j++){
        ComplexD sij=TensorRemove(innerProduct(psi_coarse[i],psi_coarse[j]));
        ComplexD d=sij-(i==j?ComplexD(1.0):ComplexD(0.0)); s2+=real(d)*real(d)+imag(d)*imag(d);
      }
      std::cout<<GridLogMessage<<"GUARD: ||<psi_coarse|psi_coarse> - I||_F = "<<std::sqrt(s2)
               <<"   (~0.23 good; ~sqrt(N_coarse)="<<std::sqrt((double)Coarse5d->gSites())<<" = e_k leak)"<<std::endl;
    }

    // --- L2->L3 single-RHS coarsening ---
    SubspaceL2 AggregatesL2(CoarseCoarse5d,Coarse5d,cb);
    for(int k=0;k<nbasis;k++) AggregatesL2.subspace[k]=psi_coarse[k];
    LittleDiracOpL2.CoarsenOperator(LinOpCoarse, AggregatesL2);
    mrhsLittleDiracOpL2.CopyMatrix(LittleDiracOpL2);
    MrhsProjectorL2.Allocate(nbasis,Coarse5d,CoarseCoarse5d);
    MrhsProjectorL2.ImportBasis(AggregatesL2.subspace);

    // --- guard psi_cc ---
    {
      std::vector<CoarseCoarseVector> psi_cc(nbasis,CoarseCoarse5d);
      for(int k=0;k<nbasis;k++) AggregatesL2.ProjectToSubspace(psi_cc[k], psi_coarse[k]);
      RealD s2=0.0;
      for(int i=0;i<nbasis;i++) for(int j=0;j<nbasis;j++){
        ComplexD sij=TensorRemove(innerProduct(psi_cc[i],psi_cc[j]));
        ComplexD d=sij-(i==j?ComplexD(1.0):ComplexD(0.0)); s2+=real(d)*real(d)+imag(d)*imag(d);
      }
      std::cout<<GridLogMessage<<"GUARD: ||<psi_cc|psi_cc> - I||_F = "<<std::sqrt(s2)
               <<"   (~0.23 good; ~sqrt(N_cc)="<<std::sqrt((double)CoarseCoarse5d->gSites())<<" = e_k leak)"<<std::endl;
    }
  } // single-RHS FINE op + padded _A + AggregatesL2 + psi_coarse freed here

  NonHermitianLinearOperator<MrhsLittleDiracOperator,CoarseVector>         mrhsLinOpCoarse(mrhsLittleDiracOpPV);
  NonHermitianLinearOperator<MrhsLittleDiracOperatorL2,CoarseCoarseVector> mrhsLinOpCC(mrhsLittleDiracOpL2);

  //////////////////////////////////////////////////////////////////////
  // DENSE coarse-coarse bottom (constructed AFTER the fine coarsening frees
  // its memory peak; needs only the hoisted single-RHS LinOpCC5d).
  //////////////////////////////////////////////////////////////////////
  std::unique_ptr<DistributedDenseInverse<CoarseCoarseVector>> DenseCC;
  std::unique_ptr<MrhsDenseCCSolve<CoarseCoarseVector>>        MrhsDenseCC;
  if (UseDenseCC) {
    std::cout << GridLogMessage << "**********************************************" << std::endl;
    std::cout << GridLogMessage << " Dense CC inverse setup (mrhs bottom)" << std::endl;
    std::cout << GridLogMessage << "**********************************************" << std::endl;
    DenseCC.reset(new DistributedDenseInverse<CoarseCoarseVector>(LinOpCC5d, CoarseCoarse5d, nbasis));
    MrhsDenseCC.reset(new MrhsDenseCCSolve<CoarseCoarseVector>(*DenseCC, CoarseCoarse5d, nrhs));
  }

  //////////////////////////////////////////////////////////////////////
  // Solvers, innermost first.
  //////////////////////////////////////////////////////////////////////
  TrivialPrecon<CoarseVector>       simpleC;
  TrivialPrecon<CoarseCoarseVector> simpleCC;
  TrivialPrecon<LatticeFermionD>    simple_fine;

  // L3 (coarse-coarse) iterative solve: PGCR on the 6D cc operator (DENSE_CC=0 branch)
  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseVector>
    L3PGCR(L3Tol,L3MaxIt,mrhsLinOpCC,simpleCC,L3Nstep,L3Nstep);
  L3PGCR.Level(3);
  L3PGCR.Name("CCouter");

  LinearFunction<CoarseCoarseVector> *ccSolve;
  if (UseDenseCC) ccSolve = MrhsDenseCC.get();
  else            ccSolve = &L3PGCR;

  // L2 coarse smoother: shifted 6D coarse op, fixed nstep
  ShiftedLinearOperator<CoarseVector> ShiftedMrhsCoarse(CoarseSmootherShift, mrhsLinOpCoarse);
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>
    CoarseSmootherGCR(0.01,1,ShiftedMrhsCoarse,simpleC,CoarseSmootherNstep,CoarseSmootherNstep);
  CoarseSmootherGCR.Level(2);
  CoarseSmootherGCR.Name("Csmoother");
  CoarseSmootherGCR.SetZeroGuess(1);   // caller zeroes vec2: skip r0 apply every L2 iteration

  // L2->L3 V-cycle preconditioner (operates on 6D coarse field)
  MrhsCoarseThreeLevelPrec<CoarseVector,CoarseCoarseVector>
    L2to3Precon(mrhsLinOpCoarse, CoarseSmootherGCR, MrhsProjectorL2, *ccSolve,
                Coarse5d, CoarseCoarse5d, CoarseCoarseMrhs, nrhs);

  // L2 coarse solve: PGCR on 6D coarse op, preconditioned by the L2->L3 V-cycle
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>
    L2PGCR(CoarseSolverTol, CoarseSolverOrder/16, mrhsLinOpCoarse, L2to3Precon, 16, 16);
  L2PGCR.Level(2);
  L2PGCR.Name("Couter");
  L2PGCR.SetZeroGuess(1);              // caller zeroes CsolMrhs; restarts still recompute r

  // Fine smoother (per-rhs, looped in the L1->L2 V-cycle)
  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD>
    SmootherGCR(0.0,1,ShiftedPVdagM,simple_fine,FineSmootherOrder,FineSmootherOrder);
  SmootherGCR.Level(1);
  SmootherGCR.Name("Fsmoother");
  SmootherGCR.SetZeroGuess(1);         // caller zeroes vec2[r]: saves 12 fine mults/outer

  // L1->L2 V-cycle (fine); its coarse solve is the three-level L2PGCR
  typedef PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> FineSmoother_t;
  MrhsTwoLevelMG<LatticeFermionD,CoarseVector,FineSmoother_t>
    ThreeLevelPrecon(PVdagM, SmootherGCR, MrhsProjector, L2PGCR, Coarse5d, CoarseMrhs);

  // Outer mrhs solve
  MrhsPGCRNonHermitian<LatticeFermionD>
    L1PGCR(OuterTol,1000,PVdagM,ThreeLevelPrecon,OuterMmax,OuterNstep);
  L1PGCR.Level(1);
  L1PGCR.Name("Fouter");
  L1PGCR.SetZeroGuess(1);              // sol[r]=Zero() below; restarts recompute r as always

  //////////////////////////////////////////////////////////////////////
  // Sources and solve
  //////////////////////////////////////////////////////////////////////
  std::vector<LatticeFermionD> src(nrhs,FGrid), sol(nrhs,FGrid);
  for(int r=0;r<nrhs;r++){ gaussian(RNG5,src[r]); sol[r]=Zero(); }

  std::cout << GridLogMessage << "**********************************************" << std::endl;
  std::cout << GridLogMessage << " MultiRHS THREE-level solve (dense CC bottom): " << nrhs << " RHS " << std::endl;
  std::cout << GridLogMessage << "**********************************************" << std::endl;

  GridStopWatch w; w.Start();
  L1PGCR(src,sol);
  w.Stop();
  std::cout << GridLogMessage << "MultiRHS 3-level dense solve total " << w.Elapsed()
            << "  (per RHS: " << w.useconds()/1.0e6/nrhs << " s)" << std::endl;

  { LatticeFermionD Ax(FGrid); RealD worst=0.0;
    for(int r=0;r<nrhs;r++){ PVdagM.Op(sol[r],Ax); Ax=Ax-src[r];
      RealD rn=std::sqrt(norm2(Ax)/norm2(src[r]));
      std::cout << GridLogMessage << "FINAL: rhs["<<r<<"] true residual = " << rn << std::endl;
      worst=std::max(worst,rn); }
    std::cout << GridLogMessage << "FINAL: worst-case residual = " << worst << std::endl;
  }

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
