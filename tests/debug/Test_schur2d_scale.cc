/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_schur2d_scale.cc

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
// SCALE rehearsal for the 2D distributed dense inverse: the full
// DENSE_SCHUR2D pipeline -- 1D rows -> redistribute -> invert ->
// redistribute back -> certificate -- on a SYNTHETIC matrix of any size,
// with no multigrid machinery, no configuration and no subspace file.
//
// This is the missing rung between the small-N oracle tests (which build
// the whole matrix on every host, impossible at production N) and the
// production example (which needs the full setup and 100 s of job time
// before the inverse is even reached).  Everything here is O(N^2/P) per
// rank; at N=138240 on 288 ranks it is the production problem shape
// exactly, in a driver that runs in minutes.
//
//   S2D_N   : global dimension          (default 720, laptop friendly)
//   S2D_NB  : block size                (default N/P rows-per-rank if that
//                                        is exact, else 48)
//
// The matrix is diagonally dominant (the recursion does not pivot); its
// conditioning is BENIGN, so this rehearses scale and speed, not the real
// operator's numerics -- the production VERIFY does that.
//
// Certificate: Cert = A0 . Ainv by the (independently validated) SUMMA,
// then every rank checks ITS OWN local elements against the identity.
// One GlobalMax at the end to report; the pipeline itself is pure P2P.
//
//   T1 : max|A.Ainv - I|  <  1e-8
//   T2 : round-trip redistribution of the INVERSE bitwise consistent
//        (CyclicToRows then RowsToCyclic reproduces the device data).
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>
#include <Grid/algorithms/multigrid/BlockCyclicSchurInverse.h>
#include <Grid/algorithms/multigrid/BlockCyclicRedistribute.h>

using namespace Grid;

// Portable |z|: ComplexD is std::complex on CPU builds and thrust::complex
// under HIP, where std::abs does not resolve (same trap RecursiveSchurInverse
// documents at FrobNorm2Local).  Member real()/imag() work on both.
static double Cabs(const ComplexD &z)
{
  double re = z.real(), im = z.imag();
  return std::sqrt(re*re + im*im);
}

static ComplexD Fill(int64_t i, int64_t j, int64_t N)
{
  double x = std::sin(0.7*i + 1.3*j);
  double y = std::cos(1.9*i - 0.4*j);
  if ( i==j ) return ComplexD(3.0*64 + x, 0.5);   // dominance independent of N
  // band-limit the off-diagonal so row sums stay bounded as N grows:
  // only |i-j| <= 64 entries are non-zero => sum |offdiag| <= 128*1.42 < 3*64
  if ( std::fabs((double)(i-j)) > 64.0 ) return ComplexD(0.0,0.0);
  return ComplexD(x,y);
}

int main(int argc, char **argv)
{
  // Environment walk (2026-08-27): the SAME inverse runs 20.8 s in the SLATE
  // harness and 26.5 s in the example, with the local GPU work (GEMM, leaf
  // inverse) 25-40% slower in the example.  Knobs to reproduce the example's
  // environment here, one at a time (systems/Frontier/schur2d_env.job):
  //   GRID_MPI_THREAD_MULTIPLE=1  MPI_Init_thread(MULTIPLE) before Grid_init (SLATE harness does this)
  //   S2D_BALLAST_GB=x            x GB of Lattice fields made device-resident before the invert
  //                               (the example carries ~10.6 GB of fine-grid state in the MemoryManager)
  //   OMP_NUM_THREADS             set in the job, read by nothing here but the runtime
  if ( getenv("GRID_MPI_THREAD_MULTIPLE") ) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    std::cout << "GRID_MPI_THREAD_MULTIPLE: requested MPI_THREAD_MULTIPLE, provided " << provided
              << (provided==MPI_THREAD_MULTIPLE ? " (MULTIPLE)" : " (NOT multiple)") << std::endl;
  }
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                       GridDefaultSimd(Nd, vComplexD::Nsimd()),
                                                       GridDefaultMpi());
  const int P  = grid->ProcessorCount();
  const int me = grid->ThisRank();

  int64_t N  = getenv("S2D_N")  ? atol(getenv("S2D_N"))  : 720;
  int64_t nb;
  if      ( getenv("S2D_NB") ) nb = atol(getenv("S2D_NB"));
  else if ( N % P == 0 )       nb = N/P;
  else                         nb = 48;
  GRID_ASSERT( N >= 1 ); GRID_ASSERT( nb >= 1 );

  int Pr,Pc;
  BlockCyclicLayout::ChooseProcessGrid(P,Pr,Pc);

  // uniform-as-possible 1D ownership, as the production import produces
  std::vector<int64_t> rowStart(P+1); rowStart[0]=0;
  for(int r=0;r<P;r++) rowStart[r+1] = rowStart[r] + N/P + ( r < (int)(N%P) ? 1 : 0 );
  int64_t myrows = rowStart[me+1]-rowStart[me];
  int64_t row0   = rowStart[me];

  std::cout << GridLogMessage << "Test_schur2d_scale: N=" << N << " nb=" << nb
            << "  grid " << Pr << "x" << Pc << "  rows/rank ~" << myrows
            << "  matrix " << (double)N*N*16.0/1.0e9 << " GB global, "
            << (double)myrows*N*16.0/1.0e9 << " GB/rank rows" << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // Fill MY rows only: O(N^2/P) host work, no global matrix anywhere.
  ////////////////////////////////////////////////////////////////////////
  double t0 = usecond();
  std::vector<ComplexD> h((uint64_t)std::max<int64_t>(myrows,1)*N);
  thread_for(jj, N, {
    for(int64_t i=0;i<myrows;i++) h[i + jj*myrows] = Fill(row0+i, jj, N);
  });
  deviceVector<ComplexD> rows1d(h.size());
  acceleratorCopyToDevice(&h[0], &rows1d[0], h.size()*sizeof(ComplexD));
  double t1 = usecond();

  ////////////////////////////////////////////////////////////////////////
  // The DENSE_SCHUR2D pipeline, phase-timed.  A0 keeps the original for
  // the certificate.
  ////////////////////////////////////////////////////////////////////////
  BlockCyclicMatrix A (grid,N,nb,Pr,Pc);
  BlockCyclicMatrix A0(grid,N,nb,Pr,Pc);
  BlockCyclicSchurInverse RSI2;

  // Pre-heat: drive the GCD with back-to-back zgemm for S2D_PREHEAT_S seconds
  // before the invert.  The example calls the inverse after ~100 s of full
  // load on all 288 GCDs and its LOCAL kernels run 20-40% slower than the
  // idle-start harness (GEMM 3.1 vs 2.5 s, leaf 0.76 vs 0.24 s) with the
  // wires unchanged; thread level / OMP / residency (E1-E5) did not reproduce
  // that.  If sustained load does, it is clock/power management, not code.
  if ( getenv("S2D_PREHEAT_S") ) {
    double secs = atof(getenv("S2D_PREHEAT_S"));
    const int64_t W = 4320;
    deviceVector<ComplexD> M((uint64_t)W*W), C((uint64_t)W*W);
    { ComplexD *m = &M[0]; accelerator_for(idx,(uint64_t)W*W,1,{ m[idx] = ComplexD(1.0e-3*(idx%97),1.0e-3*(idx%89)); }); accelerator_barrier(); }
    deviceVector<ComplexD*> ap(1),bp(1),cp(1); std::vector<ComplexD*> ptr(1);
    ptr[0]=&M[0]; acceleratorCopyToDevice(&ptr[0],&ap[0],sizeof(ComplexD*)); acceleratorCopyToDevice(&ptr[0],&bp[0],sizeof(ComplexD*));
    ptr[0]=&C[0]; acceleratorCopyToDevice(&ptr[0],&cp[0],sizeof(ComplexD*));
    double t0=usecond(); int n=0; double tlast=0;
    while ( (usecond()-t0)/1.0e6 < secs ) {
      double t1=usecond();
      RSI2.SUMMA.BLAS.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,(int)W,(int)W,(int)W,ComplexD(1.0,0.0),ap,(int)W,bp,(int)W,ComplexD(0.0,0.0),cp,(int)W);
      RSI2.SUMMA.BLAS.synchronise(); tlast=usecond()-t1; n++;
    }
    double tfirst = 0; (void)tfirst;
    std::cout << GridLogMessage << "Test_schur2d_scale: pre-heat " << (usecond()-t0)/1.0e6 << " s, " << n << " zgemm W=" << W
              << ", last zgemm " << tlast/1.0e6 << " s (" << 8.0*W*W*W/tlast/1.0e6 << " TF/s; idle-start rate 23.7)" << std::endl;
  }

  // Device ballast: Lattice fields written on the accelerator so they sit in
  // the MemoryManager's device LRU exactly as the example's fine-grid state does.
  typedef Lattice<iVector<iVector<vComplexD,Nc>,Ns> > BallastField;
  std::vector<BallastField> ballast;
  if ( getenv("S2D_BALLAST_GB") ) {
    double gb = atof(getenv("S2D_BALLAST_GB"));
    uint64_t fbytes = (uint64_t)grid->oSites()*sizeof(BallastField::vector_object);
    int nf = (int)(gb*1.0e9/(double)fbytes + 0.5);
    ballast.reserve(nf);
    for(int i=0;i<nf;i++){
      ballast.emplace_back(grid);
      autoView(v, ballast[i], AcceleratorWriteDiscard);
      accelerator_for(ss, grid->oSites(), 1, { v[ss] = Zero(); });
    }
    std::cout << GridLogMessage << "Test_schur2d_scale: device ballast " << nf << " fields x " << fbytes/1.0e6
              << " MB = " << nf*fbytes/1.0e9 << " GB resident (S2D_BALLAST_GB=" << gb << ")" << std::endl;
  }

  BlockCyclicRedistribute::RowsToCyclic(grid,rowStart,&rows1d[0],myrows,A);
  double t2 = usecond();
  if ( A.data.size() )
    acceleratorCopyDeviceToDevice((void *)&A.data[0],(void *)&A0.data[0],
                                  A.data.size()*sizeof(ComplexD));
  double t3 = usecond();
  RSI2.Invert(A);
  double t4 = usecond();
  BlockCyclicRedistribute::CyclicToRows(grid,rowStart,A,&rows1d[0],myrows);
  double t5 = usecond();
  RSI2.ReportTelemetry(grid);

  ////////////////////////////////////////////////////////////////////////
  // T1 : certificate by SUMMA, checked locally, reported by one GlobalMax.
  ////////////////////////////////////////////////////////////////////////
  int failures = 0;
  {
    BlockCyclicMatrix Cert(grid,N,nb,Pr,Pc);
    BlockCyclicSumma  SUMMA;
    SUMMA.Multiply(ComplexD(1.0,0.0),A0,A,ComplexD(0.0,0.0),Cert, 0,N,0,N,0,N);
    double t6 = usecond();

    BlockCyclicLayout &L = Cert.layout;
    std::vector<ComplexD> hc((uint64_t)std::max<int64_t>(L.mloc*L.nloc,1));
    if ( L.mloc*L.nloc )
      acceleratorCopyFromDevice(&Cert.data[0], &hc[0],
                                (uint64_t)L.mloc*L.nloc*sizeof(ComplexD));
    double mx = 0.0;
    for(int64_t lj=0;lj<L.nloc;lj++){
      int64_t gj = BlockCyclicLayout::LocalToGlobal(lj, nb, L.pcol, Pc);
      for(int64_t li=0;li<L.mloc;li++){
        int64_t gi = BlockCyclicLayout::LocalToGlobal(li, nb, L.prow, Pr);
        ComplexD id = (gi==gj) ? ComplexD(1.0,0.0) : ComplexD(0.0,0.0);
        mx = std::max(mx, Cabs(hc[li+lj*L.mloc]-id));
      }
    }
    RealD gmx = mx;  grid->GlobalMax(gmx);
    if ( gmx > 1.0e-8 ) failures++;

    std::cout << GridLogMessage << "Test_schur2d_scale phases (s):"
              << "  fill "         << (t1-t0)/1.0e6
              << "  redist->2D "   << (t2-t1)/1.0e6
              << "  invert "       << (t4-t3)/1.0e6
              << "  redist->rows " << (t5-t4)/1.0e6
              << "  certify "      << (t6-t5)/1.0e6 << std::endl;
    std::cout << GridLogMessage << "Test_schur2d_scale CERTIFICATE max|A.Ainv - I| = "
              << gmx << (gmx > 1.0e-8 ? "   ** FAIL **" : "   PASS") << std::endl;
  }

  ////////////////////////////////////////////////////////////////////////
  // T2 : the rows now hold the inverse; push them back out and compare
  // against A on device -- redistribution must be bitwise invertible on
  // real (non-synthetic-import) data too.
  ////////////////////////////////////////////////////////////////////////
  {
    BlockCyclicMatrix B(grid,N,nb,Pr,Pc);
    BlockCyclicRedistribute::RowsToCyclic(grid,rowStart,&rows1d[0],myrows,B);
    std::vector<ComplexD> x((uint64_t)std::max<int64_t>(A.layout.mloc*A.layout.nloc,1));
    std::vector<ComplexD> y(x.size());
    if ( A.layout.mloc*A.layout.nloc ){
      acceleratorCopyFromDevice(&A.data[0], &x[0], x.size()*sizeof(ComplexD));
      acceleratorCopyFromDevice(&B.data[0], &y[0], y.size()*sizeof(ComplexD));
    }
    int bad = 0;
    for(uint64_t i=0;i<x.size();i++) if ( !(x[i]==y[i]) ) bad++;
    uint64_t gbad = bad;  grid->GlobalSum(gbad);
    if ( gbad ) failures++;
    std::cout << GridLogMessage << "Test_schur2d_scale ROUND TRIP mismatches = "
              << gbad << (gbad ? "   ** FAIL **" : "   PASS") << std::endl;
  }

  {
    uint64_t f = failures;  grid->GlobalSum(f);
    failures = (int)f;
  }
  std::cout << GridLogMessage << (failures ? "Test_schur2d_scale: FAILURES"
                                           : "Test_schur2d_scale: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
