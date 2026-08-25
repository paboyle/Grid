/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_pgcr_history.cc

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
// Regression gate for the PERSISTENT history arrays in
// PrecGeneralisedConjugateResidualNonHermitian (q, p, qq allocated once per
// grid instead of per GCRnStep call).
//
//   mpirun -n 1 ./Test_pgcr_history --grid 8.8.8.8 --mpi 1.1.1.1
//
// Persistence is only safe if no stale history is ever read.  The test:
//   1. solve  src1 -> x1              (fresh allocation)
//   2. solve  src2 -> y               (DIRTIES the persistent history)
//   3. solve  src1 -> x2              (must reuse arrays, must not see src2)
//   T1 : x1 == x2 BITWISE  -- any consumption of stale history breaks this
//   T2 : true residual of x1 below tolerance -- the solver still solves
//   T3 : a second solver object with a different mmax on the same grid
//        (exercises the (re)allocation path) also converges.
//
// Wilson at mass 0.5 on a hot configuration: the Hermitian part is positive
// definite, so unpreconditioned GCR converges without breakdown.
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>

using namespace Grid;

static int failures = 0;
static void Report(const std::string &name, bool pass, const std::string &detail="")
{
  std::cout << GridLogMessage << "  " << name << (pass ? "  PASS" : "  ** FAIL **");
  if ( detail.size() ) std::cout << "   " << detail;
  std::cout << std::endl;
  if ( !pass ) failures++;
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                       GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers(seeds);

  LatticeGaugeFieldD Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4, Umu);

  RealD mass = 0.5;
  WilsonFermionD Dw(Umu, *UGrid, *UrbGrid, mass);
  NonHermitianLinearOperator<WilsonFermionD, LatticeFermionD> Op(Dw);
  TrivialPrecon<LatticeFermionD> simple;

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> PGCR(1.0e-10, 2000, Op, simple, 4, 8);
  PGCR.SetZeroGuess(1);

  LatticeFermionD src1(UGrid), src2(UGrid), x1(UGrid), x2(UGrid), y(UGrid), r(UGrid);
  gaussian(RNG4, src1);
  gaussian(RNG4, src2);

  x1 = Zero(); PGCR(src1, x1);
  y  = Zero(); PGCR(src2, y);      // dirty the history with an unrelated solve
  x2 = Zero(); PGCR(src1, x2);

  {
    r = x1 - x2;
    RealD dd = norm2(r);
    Report("T1  repeat solve after a different source is BITWISE identical", dd == 0.0,
           "|x1-x2|^2 = "+std::to_string(dd));
  }
  {
    Op.Op(x1, r); r = r - src1;
    RealD res = std::sqrt(norm2(r)/norm2(src1));
    Report("T2  true residual of the solve", res < 1.0e-9, "rel "+std::to_string(res));
  }
  {
    PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> PGCR2(1.0e-10, 2000, Op, simple, 2, 8);
    PGCR2.SetZeroGuess(1);
    PGCR2.LogCoefficients(1);          // exercise the coefficient log format
    LatticeFermionD x3(UGrid); x3 = Zero(); PGCR2(src1, x3);
    Op.Op(x3, r); r = r - src1;
    RealD res = std::sqrt(norm2(r)/norm2(src1));
    Report("T3  second solver object, different mmax, converges", res < 1.0e-9, "rel "+std::to_string(res));
  }

  std::cout << GridLogMessage << (failures ? "Test_pgcr_history: FAILURES"
                                           : "Test_pgcr_history: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
