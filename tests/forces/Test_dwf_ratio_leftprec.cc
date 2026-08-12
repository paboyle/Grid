/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/forces/Test_dwf_ratio_leftprec.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
/*  END LEGAL */

//
// Correctness of TwoFlavourRatioLeftPrecPseudoFermionAction against the
// decades-proven TwoFlavourRatioPseudoFermionAction.  Both classes compute
// the SAME action S = phi^dag V (MdagM)^-1 Vdag phi through different solve
// chains (normal-equations vs left-preconditioned F = Vdag M), so with
// twin-seeded refreshes and 1e-12 solvers they must agree to solver
// tolerance.  Tests:
//
//   E0a/E0b : heatbath identity, S == 0.5|eta|^2 after RNG refresh, for
//             BOTH classes (E0a also validates the twin-eta capture).
//   E1      : S_classic == S_leftprec       (relative, ~1e-8)
//   E2      : deriv_classic == deriv_leftprec (pointwise field norm, ~1e-8)
//   F1      : ForceTest (Test_double_ratio.cc idiom) on the LeftPrec class.
//
// All asserts are hard: this is the regression gate for the new class.
// Run small, e.g.:  ./Test_dwf_ratio_leftprec --grid 8.8.8.8
//
#include <Grid/Grid.h>
#include <Grid/qcd/action/pseudofermion/TwoFlavourRatioLeftPrec.h>

using namespace std;
using namespace Grid;

////////////////////////////////////////////////////////////////////
// Minimal LinearOperator for the composite F = Vdag M, exposing the
// Hermitian normal operator FdagF for CG.  Stencil entries assert.
////////////////////////////////////////////////////////////////////
template<class Impl>
class VdagMNormalOperator : public LinearOperatorBase<typename Impl::FermionField> {
public:
  typedef typename Impl::FermionField Field;
  FermionOperator<Impl> &VOp;
  FermionOperator<Impl> &MOp;
  VdagMNormalOperator(FermionOperator<Impl> &V,FermionOperator<Impl> &M) : VOp(V), MOp(M) {};

  void Fapply(const Field &in, Field &out) {          // out = Vdag M in
    Field tmp(in.Grid());
    MOp.M(in,tmp);
    VOp.Mdag(tmp,out);
  }
  void FdagApply(const Field &in, Field &out) {       // out = Mdag V in
    Field tmp(in.Grid());
    VOp.M(in,tmp);
    MOp.Mdag(tmp,out);
  }
  virtual void Op     (const Field &in, Field &out) { Fapply(in,out); }
  virtual void AdjOp  (const Field &in, Field &out) { FdagApply(in,out); }
  virtual void HermOp (const Field &in, Field &out) {
    Field tmp(in.Grid());
    Fapply(in,tmp);
    FdagApply(tmp,out);
  }
  virtual void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2) {
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1 = real(dot);
    n2 = norm2(out);
  }
  virtual void OpDiag (const Field &in, Field &out)                  { GRID_ASSERT(0); }
  virtual void OpDir  (const Field &in, Field &out,int dir,int disp) { GRID_ASSERT(0); }
  virtual void OpDirAll(const Field &in, std::vector<Field> &out)    { GRID_ASSERT(0); }
};

////////////////////////////////////////////////////////////////////
// F-contract LinearFunctions for the test, both via CG on FdagF:
//   forward :  F x = b     ==>  x = (FdagF)^-1 Fdag b
//   adjoint :  Fdag z = b  ==>  z = F (FdagF)^-1 b
////////////////////////////////////////////////////////////////////
template<class Impl>
class ForwardFSolve : public LinearFunction<typename Impl::FermionField> {
public:
  typedef typename Impl::FermionField Field;
  using LinearFunction<Field>::operator();
  VdagMNormalOperator<Impl> &FdagF; RealD tol; Integer maxit;
  ForwardFSolve(VdagMNormalOperator<Impl> &Op,RealD _tol,Integer _maxit) : FdagF(Op), tol(_tol), maxit(_maxit) {};
  void operator()(const Field &in, Field &out) {
    Field src(in.Grid());
    FdagF.FdagApply(in,src);
    ConjugateGradient<Field> CG(tol,maxit);
    out = Zero();
    CG(FdagF,src,out);
  }
};
template<class Impl>
class AdjointFSolve : public LinearFunction<typename Impl::FermionField> {
public:
  typedef typename Impl::FermionField Field;
  using LinearFunction<Field>::operator();
  VdagMNormalOperator<Impl> &FdagF; RealD tol; Integer maxit;
  AdjointFSolve(VdagMNormalOperator<Impl> &Op,RealD _tol,Integer _maxit) : FdagF(Op), tol(_tol), maxit(_maxit) {};
  void operator()(const Field &in, Field &out) {
    Field y(in.Grid());
    y = Zero();
    ConjugateGradient<Field> CG(tol,maxit);
    CG(FdagF,in,y);
    FdagF.Fapply(y,out);
  }
};
template<class Matrix,class Field>
class NormalEqSolve : public LinearFunction<Field> {   // out = (MdagM)^-1 in
public:
  using LinearFunction<Field>::operator();
  Matrix &_Mat; RealD tol; Integer maxit;
  NormalEqSolve(Matrix &Mat,RealD _tol,Integer _maxit) : _Mat(Mat), tol(_tol), maxit(_maxit) {};
  void operator()(const Field &in, Field &out) {
    MdagMLinearOperator<Matrix,Field> MdagM(_Mat);
    ConjugateGradient<Field> CG(tol,maxit);
    out = Zero();
    CG(MdagM,in,out);
  }
};

////////////////////////////////////////////////////////////////////
// ForceTest idiom from Test_double_ratio.cc (midpoint derivative)
////////////////////////////////////////////////////////////////////
template<class Gimpl>
void ForceTest(Action<LatticeGaugeField> &action,LatticeGaugeField & U,MomentumFilterBase<LatticeGaugeField> &Filter)
{
  GridBase *UGrid = U.Grid();

  std::vector<int> seeds({1,2,3,5});
  GridSerialRNG            sRNG;         sRNG.SeedFixedIntegers(seeds);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds);

  LatticeColourMatrix Pmu(UGrid);
  LatticeGaugeField P(UGrid);
  LatticeGaugeField UdSdU(UGrid);

  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
  std::cout << GridLogMessage << " Force test for "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "*********************************************************"<<std::endl;

  RealD eps=0.005;

  Gimpl::generate_momenta(P,sRNG,RNG4);
  Filter.applyFilter(P);

  action.refresh(U,sRNG,RNG4);

  RealD S1 = action.S(U);

  Gimpl::update_field(P,U,eps);

  action.deriv(U,UdSdU);
  UdSdU = Ta(UdSdU);
  Filter.applyFilter(UdSdU);

  DumpSliceNorm("Force",UdSdU,Nd-1);

  Gimpl::update_field(P,U,eps);

  RealD S2 = action.S(U);

  LatticeComplex dS(UGrid); dS = Zero();
  for(int mu=0;mu<Nd;mu++){
    auto UdSdUmu = PeekIndex<LorentzIndex>(UdSdU,mu);
    Pmu= PeekIndex<LorentzIndex>(P,mu);
    dS = dS - trace(Pmu*UdSdUmu)*eps*2.0*2.0;
  }
  ComplexD dSpred    = sum(dS);
  RealD diff =  S2-S1-dSpred.real();

  std::cout<< GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<< GridLogMessage << "S1 : "<< S1    <<std::endl;
  std::cout<< GridLogMessage << "S2 : "<< S2    <<std::endl;
  std::cout<< GridLogMessage << "dS : "<< S2-S1 <<std::endl;
  std::cout<< GridLogMessage << "dSpred : "<< dSpred.real() <<std::endl;
  std::cout<< GridLogMessage << "diff : "<< diff<<std::endl;
  std::cout<< GridLogMessage << "diff/dS : "<< diff/(S2-S1)<<std::endl;
  std::cout<< GridLogMessage << "*********************************************************"<<std::endl;
  std::cout<< GridLogMessage << "Done" <<std::endl;
  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout << std::setprecision(14);

  const int Ls=8;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  GridParallelRNG RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField U(UGrid);
  SU<Nc>::HotConfiguration(RNG4,U);

  ////////////////////////////////////////////////////////////////
  // Operators: quotient pair (V = PV mass 1, M light-ish), Mobius,
  // campaign b,c.  Heavyish M so CG is quick on a hot configuration.
  ////////////////////////////////////////////////////////////////
  RealD mden = 0.2;
  RealD mnum = 1.0;
  RealD M5   = 1.8;
  RealD b    = 1.5;
  RealD c    = 0.5;

  WilsonImplParams p;
  p.boundary_phases[0] = 1.0;
  p.boundary_phases[1] = 1.0;
  p.boundary_phases[2] = 1.0;
  p.boundary_phases[3] = -1.0;

  MobiusFermionD DenOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mden,M5,b,c,p);
  MobiusFermionD NumOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mnum,M5,b,c,p);

  RealD tol     = 1.0e-12;
  Integer maxit = 30000;

  typedef WilsonImplD::FermionField FermionField;

  ////////////////////////////////////////////////////////////////
  // Solvers.  Classic: CG as OperatorFunction on the MdagM linop the
  // action supplies.  LeftPrec: F-contract solves via CG on FdagF.
  ////////////////////////////////////////////////////////////////
  ConjugateGradient<FermionField> CG(tol,maxit);

  VdagMNormalOperator<WilsonImplD> FdagF(NumOp,DenOp);
  ForwardFSolve<WilsonImplD>       Ffwd (FdagF,tol,maxit);
  AdjointFSolve<WilsonImplD>       Fadj (FdagF,tol,maxit);
  NormalEqSolve<MobiusFermionD,FermionField> VdagVinv(NumOp,tol,maxit);

  TwoFlavourRatioPseudoFermionAction<WilsonImplD>         Classic (NumOp,DenOp,CG,CG);
  TwoFlavourRatioLeftPrecPseudoFermionAction<WilsonImplD> LeftPrec(NumOp,DenOp,Ffwd,Fadj,Fadj,VdagVinv);

  ////////////////////////////////////////////////////////////////
  // Twin-seeded refreshes: identical eta into both classes.
  ////////////////////////////////////////////////////////////////
  std::vector<int> seedsR({9,11,13,17});
  GridSerialRNG   sRNGa;         sRNGa.SeedFixedIntegers(seedsR);
  GridSerialRNG   sRNGb;         sRNGb.SeedFixedIntegers(seedsR);
  GridParallelRNG RNG5a(FGrid);  RNG5a.SeedFixedIntegers(seedsR);
  GridParallelRNG RNG5b(FGrid);  RNG5b.SeedFixedIntegers(seedsR);
  GridParallelRNG RNG5c(FGrid);  RNG5c.SeedFixedIntegers(seedsR);

  FermionField etaTwin(FGrid);
  gaussian(RNG5c,etaTwin);                    // identical to both refresh draws

  Classic.refresh (U,sRNGa,RNG5a);
  LeftPrec.refresh(U,sRNGb,RNG5b);

  ////////////////////////////////////////////////////////////////
  // E0 : heatbath identity for both classes
  ////////////////////////////////////////////////////////////////
  RealD Sexpect = 0.5*norm2(etaTwin);
  RealD Sc = Classic.S(U);
  RealD Sl = LeftPrec.S(U);

  RealD e0a = std::abs(Sc-Sexpect)/Sexpect;
  RealD e0b = std::abs(Sl-Sexpect)/Sexpect;
  std::cout << GridLogMessage << "=========================================================" << std::endl;
  std::cout << GridLogMessage << " E0 heatbath identity: 0.5|eta|^2 = " << Sexpect << std::endl;
  std::cout << GridLogMessage << "   classic  S = " << Sc << "  rel defect " << e0a
	    << ( e0a < 1.0e-8 ? "   PASS" : "   FAIL" ) << std::endl;
  std::cout << GridLogMessage << "   leftprec S = " << Sl << "  rel defect " << e0b
	    << ( e0b < 1.0e-8 ? "   PASS" : "   FAIL" ) << std::endl;

  ////////////////////////////////////////////////////////////////
  // E1 : action equivalence
  ////////////////////////////////////////////////////////////////
  RealD e1 = std::abs(Sc-Sl)/std::abs(Sc);
  std::cout << GridLogMessage << " E1 action equivalence: rel diff = " << e1
	    << ( e1 < 1.0e-8 ? "   PASS" : "   FAIL" ) << std::endl;

  ////////////////////////////////////////////////////////////////
  // E2 : derivative equivalence (pointwise field comparison)
  ////////////////////////////////////////////////////////////////
  LatticeGaugeField dSdUc(UGrid);
  LatticeGaugeField dSdUl(UGrid);
  LatticeGaugeField dDiff(UGrid);

  Classic.deriv (U,dSdUc);
  LeftPrec.deriv(U,dSdUl);
  dDiff = dSdUc - dSdUl;

  RealD e2 = std::sqrt( norm2(dDiff) / norm2(dSdUc) );
  std::cout << GridLogMessage << " E2 deriv equivalence: |diff|/|classic| = " << e2
	    << ( e2 < 1.0e-8 ? "   PASS" : "   FAIL" ) << std::endl;
  std::cout << GridLogMessage << "   |dSdU classic |^2 = " << norm2(dSdUc) << std::endl;
  std::cout << GridLogMessage << "   |dSdU leftprec|^2 = " << norm2(dSdUl) << std::endl;
  std::cout << GridLogMessage << "=========================================================" << std::endl;

  ////////////////////////////////////////////////////////////////
  // F1 : standalone force test on the LeftPrec class
  ////////////////////////////////////////////////////////////////
  MomentumFilterNone<LatticeGaugeField> FilterNone;
  ForceTest<GimplTypesR>(LeftPrec,U,FilterNone);

  ////////////////////////////////////////////////////////////////
  // Summary + hard asserts (this is the regression gate)
  ////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=========================================================" << std::endl;
  std::cout << GridLogMessage << " SUMMARY" << std::endl;
  std::cout << GridLogMessage << "   E0a classic  heatbath defect : " << e0a << std::endl;
  std::cout << GridLogMessage << "   E0b leftprec heatbath defect : " << e0b << std::endl;
  std::cout << GridLogMessage << "   E1  action   equivalence     : " << e1  << std::endl;
  std::cout << GridLogMessage << "   E2  deriv    equivalence     : " << e2  << std::endl;
  std::cout << GridLogMessage << "=========================================================" << std::endl;

  GRID_ASSERT(e0a < 1.0e-8);
  GRID_ASSERT(e0b < 1.0e-8);
  GRID_ASSERT(e1  < 1.0e-8);
  GRID_ASSERT(e2  < 1.0e-8);

  std::cout << GridLogMessage << "All equivalence tests PASSED" << std::endl;

  Grid_finalize();
}
