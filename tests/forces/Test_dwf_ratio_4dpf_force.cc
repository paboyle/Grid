/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/forces/Test_dwf_ratio_4dpf_force.cc

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
// Validation of TwoFlavourRatio4DPseudoFermionAction (non-EO, LinearFunction
// solver slots).  Three tests, run for BOTH wall conventions:
//
//  T1 HeatbathIdentityTest : refresh then S; PASS iff S == 0.5*|eta4|^2 to
//     solver tolerance.  This adjudicates the 4D effective-operator
//     composition identity [P M^-1 V Pdag][P V^-1 M Pdag] = 1 for the chosen
//     (P,Pdag) wall pair.  NO PREDICTION is made about which convention
//     passes -- that is what the test decides.
//  T2 ForceTest (idiom from Test_double_ratio.cc) : midpoint-derivative
//     check of deriv against S.  Should PASS for BOTH conventions (S and
//     deriv use the same literal-adjoint pair by construction).
//  T3 Trivial-ratio control (V == M) : T1 with NumOp = DenOp.  The solve
//     cancels against the multiply, so S = 0.5|eta4|^2 requires only
//     P Pdag = 1_4d.  Should PASS for BOTH conventions; isolates plumbing
//     from the composition identity.
//
// Solvers here are plain CG on the normal equations (CGNR), tolerance 1e-12,
// so every defect above ~1e-10 is structural, not solver noise.  Run small,
// e.g.:  ./Test_dwf_ratio_4dpf_force --grid 8.8.8.8
//
#include <Grid/Grid.h>
#include <Grid/qcd/action/pseudofermion/TwoFlavourRatio4DPseudoFermion.h>

using namespace std;
using namespace Grid;

////////////////////////////////////////////////////////////////////
// LinearFunction wrappers: direct M^-1 and M^-dag via CG on the
// normal equations.  These stand in for the MG-GCR stack in this
// test; the action class sees only LinearFunction.
////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class CGNRLinearFunction : public LinearFunction<Field> {   // out = M^-1 in
public:
  using LinearFunction<Field>::operator();
  Matrix &_Mat; RealD tol; Integer maxit;
  CGNRLinearFunction(Matrix &Mat,RealD _tol,Integer _maxit) : _Mat(Mat), tol(_tol), maxit(_maxit) {};
  void operator()(const Field &in, Field &out) {
    MdagMLinearOperator<Matrix,Field> MdagM(_Mat);
    Field src(in.Grid());
    _Mat.Mdag(in,src);                 // src = Mdag in
    ConjugateGradient<Field> CG(tol,maxit);
    out = Zero();
    CG(MdagM,src,out);                 // out = (MdagM)^-1 Mdag in = M^-1 in
  }
};
template<class Matrix,class Field>
class CGNRDagLinearFunction : public LinearFunction<Field> { // out = M^-dag in
public:
  using LinearFunction<Field>::operator();
  Matrix &_Mat; RealD tol; Integer maxit;
  CGNRDagLinearFunction(Matrix &Mat,RealD _tol,Integer _maxit) : _Mat(Mat), tol(_tol), maxit(_maxit) {};
  void operator()(const Field &in, Field &out) {
    MdagMLinearOperator<Matrix,Field> MdagM(_Mat);
    Field tmp(in.Grid());
    tmp = Zero();
    ConjugateGradient<Field> CG(tol,maxit);
    CG(MdagM,in,tmp);                  // tmp = (MdagM)^-1 in
    _Mat.M(tmp,out);                   // out = M (MdagM)^-1 in = M^-dag in
  }
};

////////////////////////////////////////////////////////////////////
// T1 / T3 : heatbath composition-identity test.
// Twin-seeded RNG reproduces the eta4 drawn inside refresh.
////////////////////////////////////////////////////////////////////
template<class Impl>
RealD HeatbathIdentityTest(TwoFlavourRatio4DPseudoFermionAction<Impl> &action,
			   LatticeGaugeField &U,
			   GridCartesian *UGrid,
			   const std::string &tag)
{
  typedef typename Impl::FermionField FermionField;

  std::vector<int> seeds({9,11,13,17});
  GridSerialRNG   sRNG;                 sRNG.SeedFixedIntegers(seeds);
  GridParallelRNG RNG4(UGrid);          RNG4.SeedFixedIntegers(seeds);
  GridParallelRNG RNG4check(UGrid);     RNG4check.SeedFixedIntegers(seeds);

  FermionField eta4check(UGrid);
  gaussian(RNG4check,eta4check);        // identical to the draw inside refresh

  action.refresh(U,sRNG,RNG4);
  RealD S       = action.S(U);
  RealD Sexpect = 0.5*norm2(eta4check);
  RealD defect  = std::abs(S-Sexpect)/Sexpect;

  std::cout << GridLogMessage << "=========================================================" << std::endl;
  std::cout << GridLogMessage << " HeatbathIdentityTest ["<<tag<<"]" << std::endl;
  std::cout << GridLogMessage << "   S              = " << S << std::endl;
  std::cout << GridLogMessage << "   0.5|eta4|^2    = " << Sexpect << std::endl;
  std::cout << GridLogMessage << "   relative defect = " << defect
	    << ( defect < 1.0e-8 ? "   PASS" : "   FAIL" ) << std::endl;
  std::cout << GridLogMessage << "=========================================================" << std::endl;
  return defect;
}

////////////////////////////////////////////////////////////////////
// T2 : ForceTest idiom from Test_double_ratio.cc (midpoint derivative)
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

  // Use the derivative
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
  //  GRID_ASSERT(diff<1.0);
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
  // Operators: Mobius, campaign-like b,c; heavyish masses so CGNR
  // is fast and well-conditioned even on a hot configuration.
  ////////////////////////////////////////////////////////////////
  RealD mden = 0.2;
  RealD mnum = 0.5;
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

  RealD tol   = 1.0e-12;
  Integer maxit = 20000;

  typedef WilsonImplD::FermionField FermionField;
  CGNRLinearFunction<MobiusFermionD,FermionField>    MinvSolver   (DenOp,tol,maxit);
  CGNRDagLinearFunction<MobiusFermionD,FermionField> MdagInvSolver(DenOp,tol,maxit);
  CGNRLinearFunction<MobiusFermionD,FermionField>    VinvSolver   (NumOp,tol,maxit);

  ////////////////////////////////////////////////////////////////
  // Actions: both wall conventions, plus the V==M trivial control
  ////////////////////////////////////////////////////////////////
  TwoFlavourRatio4DPseudoFermionAction<WilsonImplD> ActSol(NumOp,DenOp,MinvSolver,MdagInvSolver,MinvSolver,VinvSolver,1);
  TwoFlavourRatio4DPseudoFermionAction<WilsonImplD> ActSrc(NumOp,DenOp,MinvSolver,MdagInvSolver,MinvSolver,VinvSolver,0);
  TwoFlavourRatio4DPseudoFermionAction<WilsonImplD> ActTrivSol(DenOp,DenOp,MinvSolver,MdagInvSolver,MinvSolver,MinvSolver,1);
  TwoFlavourRatio4DPseudoFermionAction<WilsonImplD> ActTrivSrc(DenOp,DenOp,MinvSolver,MdagInvSolver,MinvSolver,MinvSolver,0);

  ////////////////////////////////////////////////////////////////
  // T3 controls first (must both pass; isolates plumbing)
  ////////////////////////////////////////////////////////////////
  RealD d3s = HeatbathIdentityTest(ActTrivSol,U,UGrid,"T3 trivial V==M, solution walls");
  RealD d3q = HeatbathIdentityTest(ActTrivSrc,U,UGrid,"T3 trivial V==M, source walls");

  ////////////////////////////////////////////////////////////////
  // T1 : the composition-identity adjudication
  ////////////////////////////////////////////////////////////////
  RealD d1s = HeatbathIdentityTest(ActSol,U,UGrid,"T1 ratio, solution walls");
  RealD d1q = HeatbathIdentityTest(ActSrc,U,UGrid,"T1 ratio, source walls");

  ////////////////////////////////////////////////////////////////
  // T2 : force consistency (expected PASS for both conventions)
  ////////////////////////////////////////////////////////////////
  MomentumFilterNone<LatticeGaugeField> FilterNone;
  ForceTest<GimplTypesR>(ActSol,U,FilterNone);
  ForceTest<GimplTypesR>(ActSrc,U,FilterNone);

  ////////////////////////////////////////////////////////////////
  // Summary
  ////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "=========================================================" << std::endl;
  std::cout << GridLogMessage << " SUMMARY (relative heatbath defects)" << std::endl;
  std::cout << GridLogMessage << "   T3 trivial  solution walls : " << d3s << std::endl;
  std::cout << GridLogMessage << "   T3 trivial  source   walls : " << d3q << std::endl;
  std::cout << GridLogMessage << "   T1 ratio    solution walls : " << d1s << std::endl;
  std::cout << GridLogMessage << "   T1 ratio    source   walls : " << d1q << std::endl;
  std::cout << GridLogMessage << " T3 must pass for both; T1 selects the wall convention." << std::endl;
  std::cout << GridLogMessage << "=========================================================" << std::endl;

  GRID_ASSERT(d3s < 1.0e-8);
  GRID_ASSERT(d3q < 1.0e-8);

  Grid_finalize();
}
