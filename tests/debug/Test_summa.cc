/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_summa.cc

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
// Regression gate for BlockCyclicSumma -- stage 2 of the 2D distributed
// dense inverse.  CPU build under mpirun:
//
//   mpirun -n 1 ./Test_summa --grid 16.16.16.32 --mpi 1.1.1.1
//   mpirun -n 2 ./Test_summa --grid 16.16.16.32 --mpi 1.1.1.2
//   mpirun -n 3 ./Test_summa --grid 16.16.16.48 --mpi 1.1.1.3
//   mpirun -n 4 ./Test_summa --grid 16.16.16.32 --mpi 1.1.1.4
//
// Every stage sweeps all process-grid factorisations of P (including the
// degenerate 1xP and Px1 rings) and a battery of (N,nb) with ragged
// trailing blocks, nb>N, and more processes than blocks.  Reference is a
// host triple loop on the replicated global matrix; tolerance 1e-11 on
// max element error (the distributed and reference summation orders
// differ, so bitwise equality is not expected AGAINST THE REFERENCE --
// but IS expected between repeated distributed runs, which is T5).
//
//   T1 : C = A.B, full range, all layouts x all grids.
//   T2 : C = beta C + alpha A.B, preloaded C, complex alpha/beta.
//   T3 : windowed products, block-aligned sub-ranges incl. ragged N end.
//   T4 : in-place windows: same matrix as A, B and C on disjoint windows.
//   T5 : determinism: repeated product bitwise identical.
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>
#include <Grid/algorithms/multigrid/BlockCyclicSumma.h>

using namespace Grid;

static int failures = 0;

// Portable |z|: ComplexD is std::complex on CPU builds and thrust::complex
// under HIP, where std::abs does not resolve (same trap RecursiveSchurInverse
// documents at FrobNorm2Local).  Member real()/imag() work on both.
static double Cabs(const ComplexD &z)
{
  double re = z.real(), im = z.imag();
  return std::sqrt(re*re + im*im);
}

static void Report(const std::string &name, bool pass, const std::string &detail="")
{
  std::cout << GridLogMessage << "  " << name << (pass ? "  PASS" : "  ** FAIL **");
  if ( detail.size() ) std::cout << "   " << detail;
  std::cout << std::endl;
  if ( !pass ) failures++;
}

// Deterministic pseudo-random fill: identical on every rank, no RNG state.
static ComplexD Fill(int64_t i, int64_t j, int salt)
{
  double x = std::sin(0.7*i + 1.3*j + 0.31*salt) ;
  double y = std::cos(1.9*i - 0.4*j + 0.77*salt) ;
  return ComplexD(x,y);
}

// Host reference:  C[i0:i1,j0:j1] = beta C + alpha A[i-,k-].B[k-,j-]
static void RefGemm(ComplexD alpha, const std::vector<ComplexD> &A,
                    const std::vector<ComplexD> &B,
                    ComplexD beta, std::vector<ComplexD> &C, int64_t N,
                    int64_t i0,int64_t i1,int64_t j0,int64_t j1,int64_t k0,int64_t k1)
{
  for(int64_t j=j0;j<j1;j++){
    for(int64_t i=i0;i<i1;i++){
      ComplexD acc(0.0,0.0);
      for(int64_t k=k0;k<k1;k++) acc += A[i+k*N]*B[k+j*N];
      C[i+j*N] = beta*C[i+j*N] + alpha*acc;
    }
  }
}

static double MaxDiff(const std::vector<ComplexD> &X, const std::vector<ComplexD> &Y)
{
  double m = 0.0;
  for(uint64_t i=0;i<X.size();i++) m = std::max(m, Cabs(X[i]-Y[i]));
  return m;
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                       GridDefaultSimd(Nd, vComplexD::Nsimd()),
                                                       GridDefaultMpi());
  const int P = grid->ProcessorCount();

  // All factorisations Pr*Pc == P: squarest, plus both degenerate rings.
  std::vector<std::pair<int,int>> grids;
  for(int r=1;r<=P;r++) if ( P%r==0 ) grids.push_back({r,P/r});

  struct Cfg { int64_t N; int64_t nb; };
  std::vector<Cfg> cfgs = { {24,4}, {26,4}, {17,5}, {8,8}, {30,7}, {3,4}, {33,3} };

  const double tol = 1.0e-11;
  BlockCyclicSumma SUMMA;

  std::cout << GridLogMessage << "BlockCyclicSumma regression: P=" << P
            << ", " << grids.size() << " process grids, " << cfgs.size()
            << " layouts" << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // T1 : full product, alpha=1 beta=0.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true; double worst = 0.0;
    for(auto &g : grids){
      for(auto &c : cfgs){
        int64_t N = c.N;
        std::vector<ComplexD> Ag(N*N), Bg(N*N), Cg(N*N, ComplexD(0.0,0.0)), Cd;
        for(int64_t j=0;j<N;j++) for(int64_t i=0;i<N;i++){
          Ag[i+j*N]=Fill(i,j,1); Bg[i+j*N]=Fill(i,j,2);
        }
        BlockCyclicMatrix A(grid,N,c.nb,g.first,g.second);
        BlockCyclicMatrix B(grid,N,c.nb,g.first,g.second);
        BlockCyclicMatrix C(grid,N,c.nb,g.first,g.second);
        A.ImportGlobal(Ag); B.ImportGlobal(Bg);
        SUMMA.Multiply(ComplexD(1.0,0.0),A,B,ComplexD(0.0,0.0),C, 0,N, 0,N, 0,N);
        C.ExportGlobal(Cd);
        RefGemm(ComplexD(1.0,0.0),Ag,Bg,ComplexD(0.0,0.0),Cg,N, 0,N,0,N,0,N);
        double d = MaxDiff(Cd,Cg); worst = std::max(worst,d);
        if ( d > tol ) ok = false;
      }
    }
    Report("T1  full C = A.B, all grids x layouts", ok,
           "max err "+std::to_string(worst));
  }

  ////////////////////////////////////////////////////////////////////////
  // T2 : beta C + alpha A.B with preloaded C, complex coefficients.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true; double worst = 0.0;
    ComplexD alpha(0.5,-0.25), beta(-1.0,0.75);
    for(auto &g : grids){
      for(auto &c : cfgs){
        int64_t N = c.N;
        std::vector<ComplexD> Ag(N*N), Bg(N*N), Cg(N*N), Cd;
        for(int64_t j=0;j<N;j++) for(int64_t i=0;i<N;i++){
          Ag[i+j*N]=Fill(i,j,3); Bg[i+j*N]=Fill(i,j,4); Cg[i+j*N]=Fill(i,j,5);
        }
        BlockCyclicMatrix A(grid,N,c.nb,g.first,g.second);
        BlockCyclicMatrix B(grid,N,c.nb,g.first,g.second);
        BlockCyclicMatrix C(grid,N,c.nb,g.first,g.second);
        A.ImportGlobal(Ag); B.ImportGlobal(Bg); C.ImportGlobal(Cg);
        SUMMA.Multiply(alpha,A,B,beta,C, 0,N, 0,N, 0,N);
        C.ExportGlobal(Cd);
        RefGemm(alpha,Ag,Bg,beta,Cg,N, 0,N,0,N,0,N);
        double d = MaxDiff(Cd,Cg); worst = std::max(worst,d);
        if ( d > tol ) ok = false;
      }
    }
    Report("T2  C = beta C + alpha A.B, preloaded C", ok,
           "max err "+std::to_string(worst));
  }

  ////////////////////////////////////////////////////////////////////////
  // T3 : windowed products.  Block-aligned sub-ranges, including the
  // ragged top end g1==N, on a layout with a partial trailing block.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true; double worst = 0.0;
    int64_t N = 26, nb = 4;                    // 6 full blocks + ragged 2
    struct Rng { int64_t i0,i1,j0,j1,k0,k1; };
    std::vector<Rng> rngs = {
      { 0,8,   8,16,  16,24 },                 // interior windows
      { 4,12,  0,4,   12,26 },                 // ragged k end
      { 16,26, 20,26, 0,8   },                 // ragged i and j ends
      { 0,4,   0,4,   4,8   },                 // minimal one-block windows
      { 0,26,  0,26,  8,12  },                 // full ij, thin k
    };
    for(auto &g : grids){
      std::vector<ComplexD> Ag(N*N), Bg(N*N), Cg(N*N), Cd;
      for(int64_t j=0;j<N;j++) for(int64_t i=0;i<N;i++){
        Ag[i+j*N]=Fill(i,j,6); Bg[i+j*N]=Fill(i,j,7);
      }
      for(auto &r : rngs){
        for(int64_t j=0;j<N;j++) for(int64_t i=0;i<N;i++) Cg[i+j*N]=Fill(i,j,8);
        BlockCyclicMatrix A(grid,N,nb,g.first,g.second);
        BlockCyclicMatrix B(grid,N,nb,g.first,g.second);
        BlockCyclicMatrix C(grid,N,nb,g.first,g.second);
        A.ImportGlobal(Ag); B.ImportGlobal(Bg); C.ImportGlobal(Cg);
        SUMMA.Multiply(ComplexD(1.0,0.0),A,B,ComplexD(1.0,0.0),C,
                       r.i0,r.i1, r.j0,r.j1, r.k0,r.k1);
        C.ExportGlobal(Cd);
        RefGemm(ComplexD(1.0,0.0),Ag,Bg,ComplexD(1.0,0.0),Cg,N,
                r.i0,r.i1, r.j0,r.j1, r.k0,r.k1);
        double d = MaxDiff(Cd,Cg); worst = std::max(worst,d);
        if ( d > tol ) ok = false;
      }
    }
    Report("T3  windowed products, ragged ends", ok,
           "max err "+std::to_string(worst));
  }

  ////////////////////////////////////////////////////////////////////////
  // T4 : in-place windows of ONE matrix, exactly the stage-3 usage:
  //   M[0:b, 2b:3b] = M[0:b, b:2b] . M[b:2b, 2b:3b]
  // C window disjoint from both operand windows (asserted in Multiply).
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true; double worst = 0.0;
    int64_t nb = 4, N = 4*nb;
    for(auto &g : grids){
      std::vector<ComplexD> Mg(N*N), Md, Mr;
      for(int64_t j=0;j<N;j++) for(int64_t i=0;i<N;i++) Mg[i+j*N]=Fill(i,j,9);
      BlockCyclicMatrix M(grid,N,nb,g.first,g.second);
      M.ImportGlobal(Mg);
      SUMMA.Multiply(ComplexD(1.0,0.0),M,M,ComplexD(0.0,0.0),M,
                     0,nb,  2*nb,3*nb,  nb,2*nb);
      M.ExportGlobal(Md);
      Mr = Mg;
      RefGemm(ComplexD(1.0,0.0),Mg,Mg,ComplexD(0.0,0.0),Mr,N,
              0,nb, 2*nb,3*nb, nb,2*nb);
      double d = MaxDiff(Md,Mr); worst = std::max(worst,d);
      if ( d > tol ) ok = false;
    }
    Report("T4  in-place disjoint windows (stage-3 usage)", ok,
           "max err "+std::to_string(worst));
  }

  ////////////////////////////////////////////////////////////////////////
  // T5 : determinism.  The summation order is fixed (ascending k-block),
  // so repeated distributed products must agree BITWISE -- the property
  // the P2P design buys and a collective reduce cannot promise.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok = true;
    for(auto &g : grids){
      int64_t N = 30, nb = 7;
      std::vector<ComplexD> Ag(N*N), Bg(N*N), C1, C2;
      for(int64_t j=0;j<N;j++) for(int64_t i=0;i<N;i++){
        Ag[i+j*N]=Fill(i,j,10); Bg[i+j*N]=Fill(i,j,11);
      }
      BlockCyclicMatrix A(grid,N,nb,g.first,g.second);
      BlockCyclicMatrix B(grid,N,nb,g.first,g.second);
      BlockCyclicMatrix C(grid,N,nb,g.first,g.second);
      A.ImportGlobal(Ag); B.ImportGlobal(Bg);
      SUMMA.Multiply(ComplexD(1.0,0.0),A,B,ComplexD(0.0,0.0),C, 0,N,0,N,0,N);
      C.ExportGlobal(C1);
      SUMMA.Multiply(ComplexD(1.0,0.0),A,B,ComplexD(0.0,0.0),C, 0,N,0,N,0,N);
      C.ExportGlobal(C2);
      for(uint64_t i=0;i<C1.size();i++)
        if ( !(C1[i]==C2[i]) ) ok = false;     // BITWISE, not toleranced
    }
    Report("T5  repeated product bitwise identical", ok);
  }

  {
    uint64_t f = failures;  grid->GlobalSum(f);
    if ( f && !failures )
      std::cout << GridLogMessage << "  ** failures on OTHER ranks: " << f << " **" << std::endl;
    failures = (int)f;
  }
  std::cout << GridLogMessage << (failures ? "Test_summa: FAILURES"
                                           : "Test_summa: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
