/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_poly_smoother.cc

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
// Fixed-polynomial smoothers (Smoothers.h) against the adaptive GCR they
// replace, on a shifted Wilson operator (non-Hermitian, spectrum in the right
// half plane):
//
//   record : GCR(mmax=2, 8 steps) with a coefficient recorder on 16 sources
//   T1     : GCRReplaySmoother on a FRESH source reduces the residual to
//            within a factor 2 of the live GCR on the same source
//   T2     : replay and GCR solutions agree to the coefficient spread (<10%)
//   T3     : ChebyshevNonHermitianSmoother with Op=HermOp reproduces the legacy
//            HermOp ChebyshevSmoother (real spectrum); on the complex-spectrum
//            Wilson op its (poor) reduction is printed for information only
//   T4     : the historical HermOp ChebyshevSmoother compiles with both
//            constructor signatures and reduces the MdagM residual
//   T5     : replay is bitwise repeatable (no reductions => no reordering)
//
//   mpirun -n 1 ./Test_poly_smoother --grid 8.8.8.8 --mpi 1.1.1.1
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>

using namespace Grid;

static int failures = 0;
static void Report(const std::string &name, bool pass, const std::string &detail="")
{
  std::cout << GridLogMessage << "  " << name << (pass ? "  PASS" : "  ** FAIL **");
  if ( detail.size() ) std::cout << "   " << detail;
  std::cout << std::endl;
  if ( !pass ) failures++;
}

template<class Field>
class ShiftedOp : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Op; RealD shift;
public:
  ShiftedOp(RealD s, LinearOperatorBase<Field> &Op) : _Op(Op), shift(s) {}
  void OpDiag  (const Field &in, Field &out) { GRID_ASSERT(0); }
  void OpDir   (const Field &in, Field &out,int dir,int disp) { GRID_ASSERT(0); }
  void OpDirAll(const Field &in, std::vector<Field> &out) { GRID_ASSERT(0); }
  void Op      (const Field &in, Field &out) { _Op.Op(in,out);    out = out + shift*in; }
  void AdjOp   (const Field &in, Field &out) { _Op.AdjOp(in,out); out = out + shift*in; }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ GRID_ASSERT(0); }
  void HermOp  (const Field &in, Field &out) { Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

template<class Field>
RealD Residual(LinearOperatorBase<Field> &Op, const Field &src, const Field &x)
{
  Field r(src.Grid()); Op.Op(x,r); r = r - src;
  return std::sqrt(norm2(r)/norm2(src));
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

  RealD mass = 0.5, shift = 0.1;
  WilsonFermionD Dw(Umu, *UGrid, *UrbGrid, mass);
  NonHermitianLinearOperator<WilsonFermionD, LatticeFermionD> Op(Dw);
  ShiftedOp<LatticeFermionD> SOp(shift, Op);
  TrivialPrecon<LatticeFermionD> simple;

  const int mmax = 2, nstep = 8, ncal = 16;

  //////////////////////////////////////////////////////////////////////
  // record
  //////////////////////////////////////////////////////////////////////
  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> GCR(0.0, 1, SOp, simple, mmax, nstep);
  GCR.SetZeroGuess(1); GCR.Name("smoother");
  GCRCoefficients rec;
  GCR.SetCoefficientRecorder(&rec);
  LatticeFermionD src(UGrid), x(UGrid);
  for(int c=0;c<ncal;c++){ gaussian(RNG4,src); x = Zero(); GCR(src,x); }
  GCR.SetCoefficientRecorder(nullptr);
  rec.Flush();
  rec.Report("smoother");
  Report("record: steps and calls", rec.Steps()==nstep && rec.Calls()==ncal,
         std::to_string(rec.Steps())+" steps, "+std::to_string(rec.Calls())+" calls");

  //////////////////////////////////////////////////////////////////////
  // fresh source: live GCR vs replay
  //////////////////////////////////////////////////////////////////////
  gaussian(RNG4,src);
  LatticeFermionD xg(UGrid), xr(UGrid), xr2(UGrid), xc(UGrid), d(UGrid);
  xg = Zero(); GCR(src,xg);
  RealD rg = Residual(SOp,src,xg);

  GCRReplaySmoother<LatticeFermionD> Replay(SOp, rec);
  Replay(src,xr);
  RealD rr = Residual(SOp,src,xr);
  Report("T1 replay residual within 2x of live GCR", rr < 2.0*rg,
         "GCR |r|/|r0| = "+std::to_string(rg)+"  replay "+std::to_string(rr));

  d = xr - xg;
  RealD rel = std::sqrt(norm2(d)/norm2(xg));
  Report("T2 replay solution vs GCR solution", rel < 0.1, "rel "+std::to_string(rel));

  //////////////////////////////////////////////////////////////////////
  // Chebyshev 1/x on [lo,hi], hi from a power iteration on SOp
  //////////////////////////////////////////////////////////////////////
  RealD hi;
  {
    LatticeFermionD v(UGrid), Av(UGrid); gaussian(RNG4,v);
    RealD n = std::sqrt(norm2(v)); v = v*(1.0/n);
    for(int i=0;i<60;i++){ SOp.Op(v,Av); hi = std::sqrt(norm2(Av)); v = Av*(1.0/hi); }
  }
  RealD lo = 0.5;
  std::cout << GridLogMessage << "power iteration |lambda_max| ~ " << hi << "  Chebyshev range [" << lo << "," << 1.05*hi << "]" << std::endl;
  ChebyshevNonHermitianSmoother<LatticeFermionD> Cheb(lo, 1.05*hi, nstep, SOp);
  Cheb(src,xc);
  RealD rc = Residual(SOp,src,xc);
  // Not gated: Wilson's spectrum has O(1) imaginary parts and a real-interval
  // Chebyshev fit to 1/x degrades exponentially off the axis (Bernstein
  // ellipse).  Printed as the reminder of what a non-real spectrum does.
  std::cout << GridLogMessage << "  (info) ChebyshevNonHermitian on the COMPLEX-spectrum Wilson op: |r|/|r0| = " << rc
            << "   (GCR " << rg << ")  -- expected poor; PVdagM smoother ops have real coefficients" << std::endl;

  //////////////////////////////////////////////////////////////////////
  // legacy HermOp Chebyshev smoother, both constructor forms, on MdagM
  //////////////////////////////////////////////////////////////////////
  {
    MdagMLinearOperator<WilsonFermionD, LatticeFermionD> HermOp(Dw);
    LatticeFermionD hsrc(UGrid), hx(UGrid), hr(UGrid); gaussian(RNG4,hsrc);
    RealD hhi = 0.0;
    { LatticeFermionD v(UGrid), Av(UGrid); gaussian(RNG4,v); RealD n=std::sqrt(norm2(v)); v=v*(1.0/n);
      for(int i=0;i<40;i++){ HermOp.HermOp(v,Av); hhi=std::sqrt(norm2(Av)); v=Av*(1.0/hhi); } }
    ChebyshevSmoother<LatticeFermionD>                 S4(0.5, 1.05*hhi, 12, HermOp);
    ChebyshevSmoother<LatticeFermionD, WilsonFermionD> S5(0.5, 1.05*hhi, 12, HermOp, Dw);   // historical 5-arg form
    S4(hsrc,hx); HermOp.HermOp(hx,hr); hr = hr - hsrc;
    RealD r4 = std::sqrt(norm2(hr)/norm2(hsrc));
    S5(hsrc,hx); HermOp.HermOp(hx,hr); hr = hr - hsrc;
    RealD r5 = std::sqrt(norm2(hr)/norm2(hsrc));
    Report("T4 legacy ChebyshevSmoother (4- and 5-arg) reduces MdagM residual", r4 < 0.5 && r5 == r4,
           "|r|/|r0| = "+std::to_string(r4)+" / "+std::to_string(r5));

    // T3: the Op()-based Clenshaw on a REAL-spectrum operator must reproduce
    // the legacy HermOp smoother: same coefficients, same recurrence.
    struct HermAsOp : public LinearOperatorBase<LatticeFermionD> {
      LinearOperatorBase<LatticeFermionD> &H;
      HermAsOp(LinearOperatorBase<LatticeFermionD> &h) : H(h) {}
      void OpDiag  (const LatticeFermionD &in, LatticeFermionD &out) { GRID_ASSERT(0); }
      void OpDir   (const LatticeFermionD &in, LatticeFermionD &out,int dir,int disp) { GRID_ASSERT(0); }
      void OpDirAll(const LatticeFermionD &in, std::vector<LatticeFermionD> &out) { GRID_ASSERT(0); }
      void Op      (const LatticeFermionD &in, LatticeFermionD &out) { H.HermOp(in,out); }
      void AdjOp   (const LatticeFermionD &in, LatticeFermionD &out) { H.HermOp(in,out); }
      void HermOpAndNorm(const LatticeFermionD &in, LatticeFermionD &out,RealD &n1,RealD &n2){ GRID_ASSERT(0); }
      void HermOp  (const LatticeFermionD &in, LatticeFermionD &out) { H.HermOp(in,out); }
    } HOp(HermOp);
    ChebyshevNonHermitianSmoother<LatticeFermionD> C3(0.5, 1.05*hhi, 12, HOp);
    LatticeFermionD hx3(UGrid), dd(UGrid);
    C3(hsrc,hx3); HermOp.HermOp(hx3,hr); hr = hr - hsrc;
    RealD r3 = std::sqrt(norm2(hr)/norm2(hsrc));
    S4(hsrc,hx); dd = hx - hx3;
    RealD reldiff = std::sqrt(norm2(dd)/norm2(hx));
    Report("T3 ChebyshevNonHermitian(Op=HermOp) == legacy ChebyshevSmoother", r3 < 0.5 && reldiff < 1.0e-12,
           "|r|/|r0| = "+std::to_string(r3)+"  rel diff of solutions "+std::to_string(reldiff));
  }

  //////////////////////////////////////////////////////////////////////
  // determinism
  //////////////////////////////////////////////////////////////////////
  Replay(src,xr2);
  d = xr - xr2;
  Report("T5 replay bitwise repeatable", norm2(d)==0.0);

  std::cout << GridLogMessage << (failures ? "Test_poly_smoother: FAILURES" : "Test_poly_smoother: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
