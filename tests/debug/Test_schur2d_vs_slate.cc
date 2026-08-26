/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_schur2d_vs_slate.cc

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

//////////////////////////////////////////////////////////////////////////////
// LIKE-FOR-LIKE: Grid's 2D block-cyclic Schur inverse versus SLATE's
// getrf+getri, on the SAME synthetic matrix, same ranks, same job, with EVERY
// layout cost on the clock.
//
//   Grid  :  rows -> RowsToCyclic -> BlockCyclicSchurInverse -> CyclicToRows
//   SLATE :  rows -> RowsToCyclic -> D2H -> fromScaLAPACK -> getrf -> getri
//                 -> H2D -> CyclicToRows
//
// Both legs pay the identical rows<->2D redistribution (the production import
// produces rank-major rows).  SLATE's EXTRA layout cost is the host round
// trip: fromScaLAPACK wraps host memory in exactly the layout BlockCyclicLayout
// already uses (csrc=0 cyclic, column-major local, ld=mloc, mb=nb, row-major
// process grid), so no other relayout exists to charge.
//
// Both inverses are certified by the SAME instrument: Grid's SUMMA product
// against an untouched copy, max |A.Ainv - I| checked locally, one GlobalMax.
//
// Precision inventory (both): fp64 throughout.  Algorithms differ: SLATE is
// LU with partial pivoting then getri; Grid is pivot-free recursive Schur.
//
// Build: needs SLATE (spack load slate) and
//   CXXFLAGS += -DHAVE_SLATE -I$SLATE_ROOT/include
//   LDFLAGS  += -L$SLATE_ROOT/lib -lslate -lblaspp -llapackpp
// Without HAVE_SLATE the SLATE leg reports itself as not built and the
// Grid leg still runs, so the binary is always usable.
//
// MPI threading: the SLATE build initialises MPI at MPI_THREAD_MULTIPLE before
// Grid_init (see main); Grid's own SERIALIZED level gives Bus errors in
// slate::listBcast with >1 OpenMP thread.  On an oversubscribed laptop set
// OMP_WAIT_POLICY=passive and OPENBLAS_NUM_THREADS=1, or SLATE's spinning
// tasks and OpenMPI's polling livelock each other (0.03 s -> minutes).
//
// Validated on the laptop (CPU, HostTask) 2026-08-25: 1x2 and 2x2 process
// grids, nb dividing N (48|720) and ragged (N=730, nb=50); both legs certify
// max|A.Ainv-I| ~ 1e-15.
//
// Frontier (HIP, Target::Devices, Cray MPICH, 2026-08-25): run with
// OMP_NUM_THREADS=1.  With 7 threads per rank SLATE's getrf deadlocks in its
// tile broadcast (concurrent device-buffer MPI from OpenMP tasks); with one
// thread it completes (8 GCDs, N=4096: certificate 2.7e-15).  One host
// thread per rank is the like-for-like anyway: Grid's GPU build uses none.
// Also: the site `slate` module (cpu env) must NOT be loaded -- its host-only
// libblaspp shadows the ROCm one via LD_LIBRARY_PATH and throws
// "device BLAS not available" from host_malloc_pinned.
//
//   S2D_N, S2D_NB as in Test_schur2d_scale (default nb = N/P).
//   S2D_SKIP_GETRI=1 skips the getri leg (host loop; ~4 min at N=138240).
//   S2D_NOWARM=1 skips the warm-up.
// A third leg, getrf+getrs(I), is SLATE's device-resident inverse route.
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>
#include <Grid/algorithms/multigrid/BlockCyclicSchurInverse.h>
#include <Grid/algorithms/multigrid/BlockCyclicRedistribute.h>
#ifdef HAVE_SLATE
#include <slate/slate.hh>
#endif

using namespace Grid;

static double Cabs(const ComplexD &z){ double re=z.real(), im=z.imag(); return std::sqrt(re*re+im*im); }

static ComplexD Fill(int64_t i, int64_t j)
{
  double x = std::sin(0.7*i + 1.3*j);
  double y = std::cos(1.9*i - 0.4*j);
  if ( i==j ) return ComplexD(3.0*64 + x, 0.5);
  if ( std::fabs((double)(i-j)) > 64.0 ) return ComplexD(0.0,0.0);
  return ComplexD(x,y);
}

// max |A0 . Ainv - I| over my local elements, reduced once.
static double Certify(GridBase *grid, BlockCyclicMatrix &A0, BlockCyclicMatrix &Ainv, int64_t nb, int Pr, int Pc)
{
  BlockCyclicMatrix Cert(grid, A0.layout.N, nb, Pr, Pc);
  BlockCyclicSumma  SUMMA;
  int64_t N = A0.layout.N;
  SUMMA.Multiply(ComplexD(1.0,0.0),A0,Ainv,ComplexD(0.0,0.0),Cert, 0,N,0,N,0,N);
  BlockCyclicLayout &L = Cert.layout;
  std::vector<ComplexD> hc((uint64_t)std::max<int64_t>(L.mloc*L.nloc,1));
  if ( L.mloc*L.nloc )
    acceleratorCopyFromDevice(&Cert.data[0], &hc[0], (uint64_t)L.mloc*L.nloc*sizeof(ComplexD));
  double mx = 0.0;
  for(int64_t lj=0;lj<L.nloc;lj++){
    int64_t gj = BlockCyclicLayout::LocalToGlobal(lj, nb, L.pcol, Pc);
    for(int64_t li=0;li<L.mloc;li++){
      int64_t gi = BlockCyclicLayout::LocalToGlobal(li, nb, L.prow, Pr);
      ComplexD id = (gi==gj) ? ComplexD(1.0,0.0) : ComplexD(0.0,0.0);
      mx = std::max(mx, Cabs(hc[li+lj*L.mloc]-id));
    }
  }
  RealD gmx = mx; grid->GlobalMax(gmx);
  return gmx;
}

int main(int argc, char **argv)
{
#ifdef HAVE_SLATE
  // SLATE issues MPI calls from concurrent OpenMP tasks and needs
  // MPI_THREAD_MULTIPLE.  Grid (GridStd.h #undef GRID_COMMS_THREADS) only asks
  // for MPI_THREAD_SERIALIZED, which reproducibly gives Bus errors inside
  // slate::BaseMatrix::listBcast with >1 OpenMP thread.  Grid_init honours an
  // already-initialised MPI, so initialise it here at the level SLATE needs.
  {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if ( provided != MPI_THREAD_MULTIPLE ) {
      fprintf(stderr, "MPI_THREAD_MULTIPLE not provided (got %d); SLATE leg unsafe\n", provided);
      GRID_ASSERT(provided == MPI_THREAD_MULTIPLE);
    }
  }
#endif
  Grid_init(&argc, &argv);
  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                       GridDefaultSimd(Nd, vComplexD::Nsimd()),
                                                       GridDefaultMpi());
  const int P  = grid->ProcessorCount();
  const int me = grid->ThisRank();

  int64_t N  = getenv("S2D_N")  ? atol(getenv("S2D_N"))  : 720;
  int64_t nb = getenv("S2D_NB") ? atol(getenv("S2D_NB")) : ( (N%P==0) ? N/P : 48 );
  int Pr,Pc; BlockCyclicLayout::ChooseProcessGrid(P,Pr,Pc);

  std::vector<int64_t> rowStart(P+1); rowStart[0]=0;
  for(int r=0;r<P;r++) rowStart[r+1] = rowStart[r] + N/P + ( r < (int)(N%P) ? 1 : 0 );
  int64_t myrows = rowStart[me+1]-rowStart[me];
  int64_t row0   = rowStart[me];

  std::cout << GridLogMessage << "Grid-vs-SLATE: N=" << N << " nb=" << nb << " grid " << Pr << "x" << Pc
            << "  matrix " << (double)N*N*16.0/1.0e9 << " GB" << std::endl;

  // my rows, identical for both legs
  std::vector<ComplexD> h((uint64_t)std::max<int64_t>(myrows,1)*N);
  thread_for(jj, N, { for(int64_t i=0;i<myrows;i++) h[i + jj*myrows] = Fill(row0+i, jj); });
  deviceVector<ComplexD> rows1d(h.size());

  ////////////////////////////////////////////////////////////////////////
  // WARM-UP: first-call costs (hipblasCreate, rocBLAS kernel load, SLATE's
  // device queues) are seconds and would be charged to whichever leg runs
  // first.  Run a small throwaway inverse through BOTH paths so the timed
  // legs below measure hot code.  Not reported.
  ////////////////////////////////////////////////////////////////////////
  // S2D_NOWARM=1 skips it (hang localisation).  Stage markers are flushed so
  // a hang shows WHERE even through block-buffered stdout.
  auto Stage = [&](const char *s){ std::cout << GridLogMessage << "stage: " << s << std::endl << std::flush; };
  if ( !getenv("S2D_NOWARM") ) {
    // Fixed tiny size independent of P: the purpose is handle creation and
    // kernel loading, not work.  (8*P at P=288 was N=2304 -> a 122 s SLATE
    // warm-up dominated by 288-way tile broadcasts.)  Ranks beyond the first
    // Nw/nbw own nothing in the rows layout, which RowsToCyclic handles.
    int64_t Nw = 64; int64_t nbw = 8;
    // same partition rule as the main leg: Nw/P rows each, remainder to the
    // first ranks; at P > Nw most ranks contribute zero rows.
    std::vector<int64_t> rs(P+1); rs[0]=0;
    for(int r=0;r<P;r++) rs[r+1] = rs[r] + Nw/P + ( r < (int)(Nw%P) ? 1 : 0 );
    int64_t rw = rs[me+1]-rs[me];
    std::vector<ComplexD> hw(std::max<int64_t>(rw,1)*Nw);
    for(int64_t jj=0;jj<Nw;jj++) for(int64_t i=0;i<rw;i++) hw[i+jj*rw]=Fill(rs[me]+i,jj);
    deviceVector<ComplexD> dw(hw.size()); acceleratorCopyToDevice(&hw[0],&dw[0],hw.size()*sizeof(ComplexD));
    BlockCyclicMatrix W(grid,Nw,nbw,Pr,Pc);
    Stage("warm-up Grid RowsToCyclic");
    BlockCyclicRedistribute::RowsToCyclic(grid,rs,&dw[0],rw,W);
    Stage("warm-up Grid Invert");
    BlockCyclicSchurInverse RSIw; RSIw.Invert(W);
#ifdef HAVE_SLATE
    {
      typedef std::complex<double> scalar_t;
      BlockCyclicLayout &L = W.layout;
      uint64_t nloc = (uint64_t)L.mloc*L.nloc;
      std::vector<scalar_t> hA(nloc ? nloc : 1);
      if ( nloc ) acceleratorCopyFromDevice(&W.data[0],(void *)&hA[0],nloc*sizeof(ComplexD));
      Stage("warm-up SLATE fromScaLAPACK");
      auto S = slate::Matrix<scalar_t>::fromScaLAPACK(Nw,Nw,&hA[0],(int64_t)std::max<int64_t>(L.mloc,1),
                                                      nbw,nbw,slate::GridOrder::Row,Pr,Pc,grid->communicator);
#if defined(GRID_HIP) || defined(GRID_CUDA) || defined(GRID_SYCL)
      slate::Options opts = { { slate::Option::Target, slate::Target::Devices } };
#else
      slate::Options opts = { { slate::Option::Target, slate::Target::HostTask } };
#endif
      slate::Pivots piv;
      Stage("warm-up SLATE getrf");
      slate::getrf(S,piv,opts);
      Stage("warm-up SLATE getri");
      slate::getri(S,piv,opts);
    }
#endif
    std::cout << GridLogMessage << "warm-up done (both paths, N=" << Nw << ")" << std::endl << std::flush;
  } else {
    std::cout << GridLogMessage << "warm-up SKIPPED (S2D_NOWARM)" << std::endl << std::flush;
  }

  ////////////////////////////////////////////////////////////////////////
  // LEG 1: Grid
  ////////////////////////////////////////////////////////////////////////
  {
    acceleratorCopyToDevice(&h[0], &rows1d[0], h.size()*sizeof(ComplexD));
    BlockCyclicMatrix A(grid,N,nb,Pr,Pc), A0(grid,N,nb,Pr,Pc);
    BlockCyclicSchurInverse RSI2;
    double t0=usecond();
    BlockCyclicRedistribute::RowsToCyclic(grid,rowStart,&rows1d[0],myrows,A);
    double t1=usecond();
    if ( A.data.size() )
      acceleratorCopyDeviceToDevice((void *)&A.data[0],(void *)&A0.data[0],A.data.size()*sizeof(ComplexD));
    double t2=usecond();
    Stage("Grid Invert");
    { GRID_TRACE("GridInvert"); RSI2.Invert(A); }
    double t3=usecond();
    RSI2.ReportTelemetry(grid);            // SUMMA ring/gemm/alloc breakdown (outside the timed window)
    BlockCyclicRedistribute::CyclicToRows(grid,rowStart,A,&rows1d[0],myrows);
    double t4=usecond();
    double cert = Certify(grid,A0,A,nb,Pr,Pc);
    std::cout << GridLogMessage << "GRID  : redist->2D " << (t1-t0)/1e6
              << "  invert " << (t3-t2)/1e6 << "  redist->rows " << (t4-t3)/1e6
              << "  TOTAL " << ((t1-t0)+(t3-t2)+(t4-t3))/1e6 << " s"
              << "   certificate " << cert << std::endl;
  }

  ////////////////////////////////////////////////////////////////////////
  // LEG 2: SLATE, every layout step timed and charged.
  ////////////////////////////////////////////////////////////////////////
#ifdef HAVE_SLATE
  if ( !getenv("S2D_SKIP_GETRI") ) {   // S2D_SKIP_GETRI=1: getri is a host loop, minutes at N=138240
    typedef std::complex<double> scalar_t;
    acceleratorCopyToDevice(&h[0], &rows1d[0], h.size()*sizeof(ComplexD));
    BlockCyclicMatrix A(grid,N,nb,Pr,Pc), A0(grid,N,nb,Pr,Pc);
    double t0=usecond();
    BlockCyclicRedistribute::RowsToCyclic(grid,rowStart,&rows1d[0],myrows,A);   // same as Grid leg
    double t1=usecond();
    if ( A.data.size() )
      acceleratorCopyDeviceToDevice((void *)&A.data[0],(void *)&A0.data[0],A.data.size()*sizeof(ComplexD));

    // SLATE's extra layout cost: host copy in the ScaLAPACK layout we already hold.
    BlockCyclicLayout &L = A.layout;
    uint64_t nloc = (uint64_t)L.mloc*L.nloc;
    std::vector<scalar_t> hA(nloc ? nloc : 1);
    double t2=usecond();
    if ( nloc ) acceleratorCopyFromDevice(&A.data[0], (void *)&hA[0], nloc*sizeof(ComplexD));
    double t3=usecond();

    // Wrap: same block size, csrc=0 cyclic ownership, column-major local
    // storage with lld = mloc, ROW-major process grid (rank = p*Pc + q),
    // exactly BlockCyclicLayout's conventions.  No further relayout exists.
    //
    // The communicator MUST be Grid's cartesian one, not MPI_COMM_WORLD:
    // BlockCyclicLayout numbers processes by grid->ThisRank(), and Grid's
    // OptimalCommunicator permutes ranks relative to the world communicator,
    // so under MPI_COMM_WORLD SLATE and Grid disagree on tile ownership
    // (observed: Bus error inside listBcast on a 2x2 grid).
    auto S = slate::Matrix<scalar_t>::fromScaLAPACK(N, N, &hA[0], (int64_t)std::max<int64_t>(L.mloc,1),
                                                    nb, nb, slate::GridOrder::Row, Pr, Pc, grid->communicator);
    // Target follows the Grid build: devices on GPU builds, host tasks on a
    // CPU build (laptop validation of the SLATE calls at small N).
#if defined(GRID_HIP) || defined(GRID_CUDA) || defined(GRID_SYCL)
    slate::Target target = slate::Target::Devices;
#else
    slate::Target target = slate::Target::HostTask;
#endif
    slate::Options opts = {
      { slate::Option::Target,        target },
      { slate::Option::Lookahead,     1 },
      { slate::Option::InnerBlocking, 16 },
    };
    slate::Pivots pivots;
    Stage("SLATE getrf");
    double t4=usecond();
    { GRID_TRACE("SLATE_getrf"); slate::getrf(S, pivots, opts); }   // LU, partial pivoting
    double t5=usecond();
    Stage("SLATE getri");
    { GRID_TRACE("SLATE_getri"); slate::getri(S, pivots, opts); }   // in-place inverse from the factor
    double t6=usecond();

    // back onto the device, in our layout (getri applies the pivots itself)
    if ( nloc ) acceleratorCopyToDevice((void *)&hA[0], &A.data[0], nloc*sizeof(ComplexD));
    double t7=usecond();
    BlockCyclicRedistribute::CyclicToRows(grid,rowStart,A,&rows1d[0],myrows);
    double t8=usecond();
    double cert = Certify(grid,A0,A,nb,Pr,Pc);
    std::cout << GridLogMessage << "SLATE : redist->2D " << (t1-t0)/1e6
              << "  D2H " << (t3-t2)/1e6 << "  wrap " << (t4-t3)/1e6
              << "  getrf " << (t5-t4)/1e6 << "  getri " << (t6-t5)/1e6
              << "  H2D " << (t7-t6)/1e6 << "  redist->rows " << (t8-t7)/1e6
              << "  TOTAL " << ((t1-t0)+(t8-t2))/1e6 << " s"
              << "   certificate " << cert << std::endl;
  }

  ////////////////////////////////////////////////////////////////////////
  // LEG 3: SLATE getrf + getrs(I) -- SLATE's device-resident route to an
  // explicit inverse.  getri's L-side loop is hard-coded Target::HostTask
  // (src/getri.cc: copy/gemmA/trsm/permuteRows all <HostTask>), so under
  // Target::Devices it runs as a 288-step host loop on one thread per rank.
  // getrs is two trsm's that honour the target (plus host row permutes).
  // The identity RHS is built on the host in the same ScaLAPACK layout and
  // its construction is on SLATE's clock.
  ////////////////////////////////////////////////////////////////////////
  {
    typedef std::complex<double> scalar_t;
    acceleratorCopyToDevice(&h[0], &rows1d[0], h.size()*sizeof(ComplexD));
    BlockCyclicMatrix A(grid,N,nb,Pr,Pc), A0(grid,N,nb,Pr,Pc);
    double t0=usecond();
    BlockCyclicRedistribute::RowsToCyclic(grid,rowStart,&rows1d[0],myrows,A);
    double t1=usecond();
    if ( A.data.size() )
      acceleratorCopyDeviceToDevice((void *)&A.data[0],(void *)&A0.data[0],A.data.size()*sizeof(ComplexD));
    BlockCyclicLayout &L = A.layout;
    uint64_t nloc = (uint64_t)L.mloc*L.nloc;
    std::vector<scalar_t> hA(nloc ? nloc : 1), hB(nloc ? nloc : 1);
    double t2=usecond();
    if ( nloc ) acceleratorCopyFromDevice(&A.data[0], (void *)&hA[0], nloc*sizeof(ComplexD));
    for(int64_t lj=0;lj<L.nloc;lj++){
      int64_t gj = BlockCyclicLayout::LocalToGlobal(lj, nb, L.pcol, Pc);
      for(int64_t li=0;li<L.mloc;li++){
        int64_t gi = BlockCyclicLayout::LocalToGlobal(li, nb, L.prow, Pr);
        hB[li+lj*L.mloc] = (gi==gj) ? scalar_t(1.0,0.0) : scalar_t(0.0,0.0);
      }
    }
    double t3=usecond();
    auto S = slate::Matrix<scalar_t>::fromScaLAPACK(N, N, &hA[0], (int64_t)std::max<int64_t>(L.mloc,1),
                                                    nb, nb, slate::GridOrder::Row, Pr, Pc, grid->communicator);
    auto B = slate::Matrix<scalar_t>::fromScaLAPACK(N, N, &hB[0], (int64_t)std::max<int64_t>(L.mloc,1),
                                                    nb, nb, slate::GridOrder::Row, Pr, Pc, grid->communicator);
#if defined(GRID_HIP) || defined(GRID_CUDA) || defined(GRID_SYCL)
    slate::Target target = slate::Target::Devices;
#else
    slate::Target target = slate::Target::HostTask;
#endif
    slate::Options opts = {
      { slate::Option::Target,        target },
      { slate::Option::Lookahead,     1 },
      { slate::Option::InnerBlocking, 16 },
    };
    slate::Pivots pivots;
    Stage("SLATE getrf (getrs leg)");
    double t4=usecond();
    { GRID_TRACE("SLATE_getrf"); slate::getrf(S, pivots, opts); }
    double t5=usecond();
    Stage("SLATE getrs(I)");
    { GRID_TRACE("SLATE_getrs"); slate::getrs(S, pivots, B, opts); }   // B <- A^{-1} I
    double t6=usecond();
    if ( nloc ) acceleratorCopyToDevice((void *)&hB[0], &A.data[0], nloc*sizeof(ComplexD));
    double t7=usecond();
    BlockCyclicRedistribute::CyclicToRows(grid,rowStart,A,&rows1d[0],myrows);
    double t8=usecond();
    double cert = Certify(grid,A0,A,nb,Pr,Pc);
    std::cout << GridLogMessage << "SLATE-getrs : redist->2D " << (t1-t0)/1e6
              << "  D2H+I " << (t3-t2)/1e6 << "  wrap " << (t4-t3)/1e6
              << "  getrf " << (t5-t4)/1e6 << "  getrs " << (t6-t5)/1e6
              << "  H2D " << (t7-t6)/1e6 << "  redist->rows " << (t8-t7)/1e6
              << "  TOTAL " << ((t1-t0)+(t8-t2))/1e6 << " s"
              << "   certificate " << cert << std::endl;
  }
#else
  std::cout << GridLogMessage << "SLATE : leg not built (compile with -DHAVE_SLATE and link -lslate -lblaspp -llapackpp)" << std::endl;
#endif

  Grid_finalize();
  return 0;
}
