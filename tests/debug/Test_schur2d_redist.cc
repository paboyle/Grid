/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_schur2d_redist.cc

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
// Regression gate for BlockCyclicRedistribute and the full stage-4 path
//
//   1D rank-major rows -> block cyclic -> Invert -> back to 1D rows
//
// which is exactly what DENSE_SCHUR2D runs inside DenseCoarseMatrix.
// CPU build under mpirun at n = 1,2,3,4.
//
//   T1 : RowsToCyclic against a direct ImportGlobal of the same matrix --
//        BITWISE (pure data movement, no arithmetic).
//   T2 : round trip rows -> 2D -> rows -- BITWISE, uniform AND non-uniform
//        rowStart, layouts with ragged trailing blocks.
//   T3 : full pipeline inverse against a host Gauss-Jordan reference.
//   T4 : CROSS-IMPLEMENTATION: the same matrix inverted by the 1D
//        RecursiveSchurInverse and by the 2D pipeline; results compared
//        element-wise.  Two independent implementations, two independent
//        decompositions, one answer.
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>
#include <Grid/algorithms/multigrid/RecursiveSchurInverse.h>
#include <Grid/algorithms/multigrid/BlockCyclicSchurInverse.h>
#include <Grid/algorithms/multigrid/BlockCyclicRedistribute.h>

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

static ComplexD Fill(int64_t i, int64_t j, int salt)
{
  double x = std::sin(0.7*i + 1.3*j + 0.31*salt);
  double y = std::cos(1.9*i - 0.4*j + 0.77*salt);
  return ComplexD(x,y);
}
static void MakeMatrix(std::vector<ComplexD> &A, int64_t N, int salt)
{
  A.resize((uint64_t)N*N);
  for(int64_t j=0;j<N;j++)
    for(int64_t i=0;i<N;i++)
      A[i+j*N] = Fill(i,j,salt) + ((i==j) ? ComplexD(3.0*N,0.5) : ComplexD(0.0,0.0));
}
static void HostInverse(std::vector<ComplexD> A, std::vector<ComplexD> &X, int64_t N)
{
  X.assign((uint64_t)N*N, ComplexD(0.0,0.0));
  for(int64_t i=0;i<N;i++) X[i+i*N] = ComplexD(1.0,0.0);
  for(int64_t c=0;c<N;c++){
    int64_t piv=c; double mx = Cabs(A[c+c*N]);
    for(int64_t r=c+1;r<N;r++)
      if ( Cabs(A[r+c*N]) > mx ){ mx=Cabs(A[r+c*N]); piv=r; }
    GRID_ASSERT( mx > 0.0 );
    if ( piv != c )
      for(int64_t j=0;j<N;j++){ std::swap(A[c+j*N],A[piv+j*N]); std::swap(X[c+j*N],X[piv+j*N]); }
    ComplexD d = ComplexD(1.0,0.0)/A[c+c*N];
    for(int64_t j=0;j<N;j++){ A[c+j*N]*=d; X[c+j*N]*=d; }
    for(int64_t r=0;r<N;r++){
      if ( r==c ) continue;
      ComplexD f = A[r+c*N];
      if ( f == ComplexD(0.0,0.0) ) continue;
      for(int64_t j=0;j<N;j++){ A[r+j*N]-=f*A[c+j*N]; X[r+j*N]-=f*X[c+j*N]; }
    }
  }
}

// Ownership tables over the ranks: uniform-ish, and deliberately lopsided.
static std::vector<int64_t> MakeRowStart(int64_t N, int P, int lopsided)
{
  std::vector<int64_t> t(P+1); t[0]=0;
  for(int r=0;r<P;r++){
    int64_t base = N/P;
    int64_t extra = ( r < (int)(N%P) ) ? 1 : 0;
    int64_t n = base + extra;
    if ( lopsided ){                       // shift work toward low ranks
      if ( r==0 && P>1 ) n = std::min(N, base+extra+3);
      if ( r==P-1 )      n = N - t[r];     // remainder
    }
    t[r+1] = std::min(N, t[r]+n);
  }
  t[P]=N;
  return t;
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                       GridDefaultSimd(Nd, vComplexD::Nsimd()),
                                                       GridDefaultMpi());
  const int P  = grid->ProcessorCount();
  const int me = grid->ThisRank();

  std::vector<std::pair<int,int>> grids;
  for(int r=1;r<=P;r++) if ( P%r==0 ) grids.push_back({r,P/r});

  struct Cfg { int64_t N; int64_t nb; };
  std::vector<Cfg> cfgs = { {24,4}, {26,4}, {17,5}, {30,7}, {33,3} };

  std::cout << GridLogMessage << "BlockCyclicRedistribute regression: P=" << P << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // T1 + T2 : distribution correctness and round trip, bitwise.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok1 = true, ok2 = true;
    for(int lop=0; lop<2; lop++){
      for(auto &g : grids){
        for(auto &c : cfgs){
          int64_t N = c.N;
          std::vector<int64_t> rowStart = MakeRowStart(N,P,lop);
          int64_t myrows = rowStart[me+1]-rowStart[me];

          std::vector<ComplexD> Ag;  MakeMatrix(Ag, N, 21);

          // my 1D rows, column major ld = myrows
          std::vector<ComplexD> h((uint64_t)std::max<int64_t>(myrows,1)*N);
          for(int64_t j=0;j<N;j++)
            for(int64_t i=0;i<myrows;i++)
              h[i + j*myrows] = Ag[(rowStart[me]+i) + j*N];
          deviceVector<ComplexD> rows1d(h.size());
          acceleratorCopyToDevice(&h[0], &rows1d[0], h.size()*sizeof(ComplexD));

          BlockCyclicMatrix A(grid,N,c.nb,g.first,g.second);
          BlockCyclicMatrix R(grid,N,c.nb,g.first,g.second);
          BlockCyclicRedistribute::RowsToCyclic(grid,rowStart,&rows1d[0],myrows,A);
          R.ImportGlobal(Ag);

          // T1: bitwise against direct import
          {
            std::vector<ComplexD> x((uint64_t)A.layout.mloc*A.layout.nloc);
            std::vector<ComplexD> y(x.size());
            if ( x.size() ){
              acceleratorCopyFromDevice(&A.data[0], &x[0], x.size()*sizeof(ComplexD));
              acceleratorCopyFromDevice(&R.data[0], &y[0], y.size()*sizeof(ComplexD));
            }
            for(uint64_t i=0;i<x.size();i++) if ( !(x[i]==y[i]) ) ok1 = false;
          }
          // T2: round trip, bitwise
          {
            deviceVector<ComplexD> back(h.size());
            std::vector<ComplexD>  hb(h.size(), ComplexD(0.0,0.0));
            acceleratorCopyToDevice(&hb[0], &back[0], hb.size()*sizeof(ComplexD));
            BlockCyclicRedistribute::CyclicToRows(grid,rowStart,A,&back[0],myrows);
            acceleratorCopyFromDevice(&back[0], &hb[0], hb.size()*sizeof(ComplexD));
            for(int64_t j=0;j<N;j++)
              for(int64_t i=0;i<myrows;i++)
                if ( !(hb[i+j*myrows]==h[i+j*myrows]) ) ok2 = false;
          }
        }
      }
    }
    Report("T1  RowsToCyclic == ImportGlobal, bitwise", ok1);
    Report("T2  round trip rows->2D->rows, bitwise, incl. lopsided rowStart", ok2);
  }

  ////////////////////////////////////////////////////////////////////////
  // T3 + T4 : the DENSE_SCHUR2D pipeline against the host reference and
  // against the INDEPENDENT 1D RecursiveSchurInverse.
  ////////////////////////////////////////////////////////////////////////
  {
    bool ok3 = true, ok4 = true;
    double worst3 = 0.0, worst4 = 0.0;
    BlockCyclicSchurInverse RSI2;
    for(auto &g : grids){
      for(auto &c : cfgs){
        int64_t N = c.N;
        std::vector<int64_t> rowStart = MakeRowStart(N,P,0);
        int64_t myrows = rowStart[me+1]-rowStart[me];

        std::vector<ComplexD> Ag, Ref;
        MakeMatrix(Ag, N, 22);
        HostInverse(Ag, Ref, N);

        std::vector<ComplexD> h((uint64_t)std::max<int64_t>(myrows,1)*N);
        for(int64_t j=0;j<N;j++)
          for(int64_t i=0;i<myrows;i++)
            h[i + j*myrows] = Ag[(rowStart[me]+i) + j*N];

        double mxref = 0.0;
        for(auto &z : Ref) mxref = std::max(mxref, Cabs(z));

        // ---- 2D pipeline: rows -> cyclic -> invert -> rows ----
        std::vector<ComplexD> h2d(h.size());
        {
          deviceVector<ComplexD> rows1d(h.size());
          acceleratorCopyToDevice(&h[0], &rows1d[0], h.size()*sizeof(ComplexD));
          BlockCyclicMatrix A(grid,N,c.nb,g.first,g.second);
          BlockCyclicRedistribute::RowsToCyclic(grid,rowStart,&rows1d[0],myrows,A);
          RSI2.Invert(A);
          BlockCyclicRedistribute::CyclicToRows(grid,rowStart,A,&rows1d[0],myrows);
          acceleratorCopyFromDevice(&rows1d[0], &h2d[0], h2d.size()*sizeof(ComplexD));
        }
        for(int64_t j=0;j<N;j++)
          for(int64_t i=0;i<myrows;i++){
            double d = Cabs(h2d[i+j*myrows]-Ref[(rowStart[me]+i)+j*N])/mxref;
            worst3 = std::max(worst3,d);
            if ( d > 1.0e-9 ) ok3 = false;
          }

        // ---- 1D RecursiveSchurInverse on the same matrix ----
        {
          BlockRows Ar;  Ar.Resize(myrows, N);
          acceleratorCopyToDevice(&h[0], &Ar.data[0], h.size()*sizeof(ComplexD));
          std::vector<int64_t> rs = rowStart;
          RecursiveSchurInverse RSI1(grid, N, rs, 1<<20);
          RSI1.Invert(Ar);
          std::vector<ComplexD> h1d(h.size());
          acceleratorCopyFromDevice(&Ar.data[0], &h1d[0], h1d.size()*sizeof(ComplexD));
          for(int64_t j=0;j<N;j++)
            for(int64_t i=0;i<myrows;i++){
              double d = Cabs(h2d[i+j*myrows]-h1d[i+j*myrows])/mxref;
              worst4 = std::max(worst4,d);
              if ( d > 1.0e-9 ) ok4 = false;
            }
        }
      }
    }
    Report("T3  2D pipeline vs host reference", ok3, "worst "+std::to_string(worst3));
    Report("T4  2D pipeline vs 1D RecursiveSchurInverse", ok4, "worst "+std::to_string(worst4));
  }

  {
    uint64_t f = failures;  grid->GlobalSum(f);
    if ( f && !failures )
      std::cout << GridLogMessage << "  ** failures on OTHER ranks: " << f << " **" << std::endl;
    failures = (int)f;
  }
  std::cout << GridLogMessage << (failures ? "Test_schur2d_redist: FAILURES"
                                           : "Test_schur2d_redist: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
