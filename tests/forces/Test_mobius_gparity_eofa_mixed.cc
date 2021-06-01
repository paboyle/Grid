/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/forces/Test_mobius_gparity_eofa_mixed.cc

Copyright (C) 2017

Author: Christopher Kelly <ckelly@bnl.gov>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

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

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
 ;

typedef GparityWilsonImplD FermionImplPolicyD;
typedef GparityMobiusEOFAFermionD FermionActionD;
typedef typename FermionActionD::FermionField FermionFieldD;

typedef GparityWilsonImplF FermionImplPolicyF;
typedef GparityMobiusEOFAFermionF FermionActionF;
typedef typename FermionActionF::FermionField FermionFieldF;

NAMESPACE_BEGIN(Grid);

  template<class FermionOperatorD, class FermionOperatorF, class SchurOperatorD, class  SchurOperatorF> 
  class MixedPrecisionConjugateGradientOperatorFunction : public OperatorFunction<typename FermionOperatorD::FermionField> {
  public:
    typedef typename FermionOperatorD::FermionField FieldD;
    typedef typename FermionOperatorF::FermionField FieldF;

    using OperatorFunction<FieldD>::operator();

    RealD   Tolerance;
    RealD   InnerTolerance; //Initial tolerance for inner CG. Defaults to Tolerance but can be changed
    Integer MaxInnerIterations;
    Integer MaxOuterIterations;
    GridBase* SinglePrecGrid4; //Grid for single-precision fields
    GridBase* SinglePrecGrid5; //Grid for single-precision fields
    RealD OuterLoopNormMult; //Stop the outer loop and move to a final double prec solve when the residual is OuterLoopNormMult * Tolerance

    FermionOperatorF &FermOpF;
    FermionOperatorD &FermOpD;;
    SchurOperatorF &LinOpF;
    SchurOperatorD &LinOpD;

    Integer TotalInnerIterations; //Number of inner CG iterations
    Integer TotalOuterIterations; //Number of restarts
    Integer TotalFinalStepIterations; //Number of CG iterations in final patch-up step

    MixedPrecisionConjugateGradientOperatorFunction(RealD tol, 
						    Integer maxinnerit, 
						    Integer maxouterit, 
						    GridBase* _sp_grid4, 
						    GridBase* _sp_grid5, 
						    FermionOperatorF &_FermOpF,
						    FermionOperatorD &_FermOpD,
						    SchurOperatorF   &_LinOpF,
						    SchurOperatorD   &_LinOpD): 
      LinOpF(_LinOpF),
      LinOpD(_LinOpD),
      FermOpF(_FermOpF),
      FermOpD(_FermOpD),
      Tolerance(tol), 
      InnerTolerance(tol), 
      MaxInnerIterations(maxinnerit), 
      MaxOuterIterations(maxouterit), 
      SinglePrecGrid4(_sp_grid4),
      SinglePrecGrid5(_sp_grid5),
      OuterLoopNormMult(100.) 
    { 
    };

    void operator()(LinearOperatorBase<FieldD> &LinOpU, const FieldD &src, FieldD &psi) {

      std::cout << GridLogMessage << " Mixed precision CG wrapper operator() "<<std::endl;

      SchurOperatorD * SchurOpU = static_cast<SchurOperatorD *>(&LinOpU);
      assert(&(SchurOpU->_Mat)==&(LinOpD._Mat));

      ////////////////////////////////////////////////////////////////////////////////////
      // Must snarf a single precision copy of the gauge field in Linop_d argument
      ////////////////////////////////////////////////////////////////////////////////////
      //typedef typename FermionOperatorF::GaugeField GaugeFieldF;
      //typedef typename FermionOperatorF::GaugeLinkField GaugeLinkFieldF;
      //typedef typename FermionOperatorD::GaugeField GaugeFieldD;
      //typedef typename FermionOperatorD::GaugeLinkField GaugeLinkFieldD;

      //GridBase * GridPtrF = SinglePrecGrid4;
      //GridBase * GridPtrD = FermOpD.Umu.Grid();
      //GaugeFieldF     U_f  (GridPtrF);
      //GaugeLinkFieldF Umu_f(GridPtrF);

      ////////////////////////////////////////////////////////////////////////////////////
      // Moving this to a Clone method of fermion operator would allow to duplicate the 
      // physics parameters and decrease gauge field copies
      ////////////////////////////////////////////////////////////////////////////////////

      //typedef typename std::decay<decltype(PeekIndex<LorentzIndex>(FermOpD.Umu, 0))>::type DoubleS

      //GaugeLinkFieldD Umu_d(GridPtrD);
      //for(int mu=0;mu<Nd*2;mu++){ 
      //Umu_d = PeekIndex<LorentzIndex>(FermOpD.Umu, mu);
      //precisionChange(Umu_f,Umu_d);
      //PokeIndex<LorentzIndex>(FermOpF.Umu, Umu_f, mu);
      //}

      precisionChange(FermOpF.Umu, FermOpD.Umu);

      pickCheckerboard(Even,FermOpF.UmuEven,FermOpF.Umu);
      pickCheckerboard(Odd ,FermOpF.UmuOdd ,FermOpF.Umu);

      ////////////////////////////////////////////////////////////////////////////////////
      // Make a mixed precision conjugate gradient
      ////////////////////////////////////////////////////////////////////////////////////
      MixedPrecisionConjugateGradient<FieldD,FieldF> MPCG(Tolerance,MaxInnerIterations,MaxOuterIterations,SinglePrecGrid5,LinOpF,LinOpD);
      MPCG.InnerTolerance = InnerTolerance;
      std::cout << GridLogMessage << "Calling mixed precision Conjugate Gradient" <<std::endl;
      MPCG(src,psi);
    }
  };

NAMESPACE_END(Grid);



int main (int argc, char** argv)
{
  Grid_init(&argc, &argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate mpi_layout  = GridDefaultMpi();

  const int Ls = 8;

  GridCartesian         *UGridD   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGridD = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridD);
  GridCartesian         *FGridD   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGridD);
  GridRedBlackCartesian *FrbGridD = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGridD);

  GridCartesian         *UGridF   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian         *FGridF   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGridF);
  GridRedBlackCartesian *FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGridF);

  std::vector<int> seeds4({1,2,3,5});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG5(FGridD);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGridD);  RNG4.SeedFixedIntegers(seeds4);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  LatticeGaugeFieldD Ud(UGridD);
  SU<Nc>::HotConfiguration(RNG4,Ud);

  LatticeGaugeFieldF Uf(UGridF);
  precisionChange(Uf, Ud);

  RealD b  = 2.5;
  RealD c  = 1.5;
  RealD mf = 0.01;
  RealD mb = 1.0;
  RealD M5 = 1.8;
  FermionActionD::ImplParams params;
  params.twists[0] = 1; //GPBC in X
  params.twists[Nd-1] = 1; //APRD in T

  std::vector<int> gtwists(4,0);
  gtwists[0] = 1;

  ConjugateGimplD::setDirections(gtwists);

  FermionActionD LopD(Ud, *FGridD, *FrbGridD, *UGridD, *UrbGridD, mf, mf, mb, 0.0, -1, M5, b, c, params);
  FermionActionD RopD(Ud, *FGridD, *FrbGridD, *UGridD, *UrbGridD, mb, mf, mb, -1.0, 1, M5, b, c, params);

  FermionActionF LopF(Uf, *FGridF, *FrbGridF, *UGridF, *UrbGridF, mf, mf, mb, 0.0, -1, M5, b, c, params);
  FermionActionF RopF(Uf, *FGridF, *FrbGridF, *UGridF, *UrbGridF, mb, mf, mb, -1.0, 1, M5, b, c, params);


  OneFlavourRationalParams OFRp(0.95, 100.0, 5000, 1.0e-12, 12);
  ConjugateGradient<FermionFieldD> CG(1.0e-10, 10000);


  typedef SchurDiagMooeeOperator<FermionActionD,FermionFieldD> EOFAschuropD;
  typedef SchurDiagMooeeOperator<FermionActionF,FermionFieldF> EOFAschuropF;
  
  EOFAschuropD linopL_D(LopD);
  EOFAschuropD linopR_D(RopD);

  EOFAschuropF linopL_F(LopF);
  EOFAschuropF linopR_F(RopF);

  typedef MixedPrecisionConjugateGradientOperatorFunction<FermionActionD, FermionActionF, EOFAschuropD, EOFAschuropF> EOFA_mxCG;

  EOFA_mxCG MCG_L(1e-10, 10000, 1000, UGridF, FrbGridF, LopF, LopD, linopL_F, linopL_D);
  MCG_L.InnerTolerance = 1e-5;

  EOFA_mxCG MCG_R(1e-10, 10000, 1000, UGridF, FrbGridF, RopF, RopD, linopR_F, linopR_D);
  MCG_R.InnerTolerance = 1e-5;

  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicyD> MeofaD(LopD, RopD, CG, CG, CG, CG, CG, OFRp, true);
  ExactOneFlavourRatioMixedPrecHeatbathPseudoFermionAction<FermionImplPolicyD, FermionImplPolicyF> MeofaMx(LopF, RopF, LopD, RopD, MCG_L, MCG_R, MCG_L, MCG_R, MCG_L, MCG_R, OFRp, true);
  
  FermionFieldD eta(FGridD);
  gaussian(RNG5, eta);

  MeofaD.refresh(Ud, eta);
  MeofaMx.refresh(Ud, eta);

  FermionFieldD diff_phi(FGridD);
  diff_phi = MeofaD.getPhi() - MeofaMx.getPhi();

  RealD n = norm2(diff_phi);
  
  std::cout << GridLogMessage << "Phi(double)=" << norm2(MeofaD.getPhi()) << " Phi(mixed)=" << norm2(MeofaMx.getPhi()) << " diff=" << n << std::endl;

  assert(n < 1e-8);

  RealD Sd = MeofaD.S(Ud);
  RealD Smx = MeofaMx.S(Ud);

  std::cout << GridLogMessage << "Initial action double=" << Sd << " mixed=" << Smx << " diff=" << Sd-Smx << std::endl;

  assert(fabs(Sd-Smx) < 1e-6);

  SU<Nc>::HotConfiguration(RNG4,Ud);
  precisionChange(Uf, Ud);

  Sd = MeofaD.S(Ud);
  Smx = MeofaMx.S(Ud);

  std::cout << GridLogMessage << "After randomizing U, action double=" << Sd << " mixed=" << Smx << " diff=" << Sd-Smx << std::endl;

  assert(fabs(Sd-Smx) < 1e-6);

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
}
