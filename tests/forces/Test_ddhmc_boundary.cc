   /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_ddhmc_boundary.cc

    Copyright (C) 2021

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/qcd/action/momentum/DirichletFilter.h>
#include <Grid/qcd/action/momentum/DDHMCfilter.h>
#include <Grid/qcd/action/fermion/DirichletFermionOperator.h>
#include <Grid/qcd/action/fermion/SchurFactoredFermionOperator.h>
#include <Grid/qcd/action/pseudofermion/DomainDecomposedBoundaryTwoFlavourPseudoFermion.h>
#include <Grid/qcd/action/pseudofermion/DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion.h>
#include <Grid/qcd/action/pseudofermion/DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion.h>

using namespace std;
using namespace Grid;

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
  
  RealD eps=0.001;

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Refresh "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  Gimpl::generate_momenta(P,sRNG,RNG4);
  Filter.applyFilter(P);

  SU<Nc>::HotConfiguration(RNG4,U);

  action.refresh(U,sRNG,RNG4);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Action "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

  RealD S1 = action.S(U);

  Gimpl::update_field(P,U,eps);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Derivative "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  action.deriv(U,UdSdU);
  UdSdU = Ta(UdSdU);
  Filter.applyFilter(UdSdU);

  DumpSliceNorm("Force",UdSdU,Nd-1);
  
  Gimpl::update_field(P,U,eps);
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Action "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
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
  std::cout<< GridLogMessage << "*********************************************************"<<std::endl;
  assert(diff<1.0);
  std::cout<< GridLogMessage << "Done" <<std::endl;
  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  LatticeGaugeField U(UGrid);

  RealD beta=6.0;
  WilsonGaugeActionR PlaqAction(beta);
  IwasakiGaugeActionR RectAction(beta);

  MomentumFilterNone<LatticeGaugeField> FilterNone;
  ForceTest<GimplTypesR>(PlaqAction,U,FilterNone);
  ForceTest<GimplTypesR>(RectAction,U,FilterNone);

  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=0.01; 
  RealD M5=1.8; 

  typedef DirichletFermionOperator<WilsonImplR> DirichletFermion;

  Coordinate Block({16,16,16,8});

  DomainWallFermionR DdwfPeriodic(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionR Ddwf(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DirichletFermion   DdwfDirichlet(Ddwf,Block);

  DomainWallFermionR PVPeriodic(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,M5);
  DomainWallFermionR PV(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,M5);
  DirichletFermion   PVDirichlet(PV,Block);

  double StoppingCondition = 1.0e-12;
  double MaxCGIterations = 10000;
  ConjugateGradient<LatticeFermion>  CG(StoppingCondition,MaxCGIterations);

  //  DDHMCFilter<LatticeGaugeField> FilterDDHMC(Block,1);
  DirichletFilter<LatticeGaugeField> FilterDDHMC(Block);

  //////////////////// Two Flavour Determinant Ratio ///////////////////////////////
  typedef WilsonImplR FermionImplPolicy;
  TwoFlavourRatioPseudoFermionAction<FermionImplPolicy> Nf2(PVPeriodic, DdwfPeriodic,CG,CG);
  ForceTest<GimplTypesR>(Nf2,U,FilterNone);

  //////////////////// Two Flavour Determinant force test Even Odd ///////////////////////////////
  TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> Nf2eo(PVPeriodic, DdwfPeriodic,CG,CG);
  ForceTest<GimplTypesR>(Nf2eo,U,FilterNone);
  
  //////////////////// DDHMC Boundary force ///////////////////////////////
  
  SchurFactoredFermionOperator<DomainWallFermionR::Impl_t>
    SchurDwf(DdwfPeriodic,
	     DdwfDirichlet,
	     CG,Block);

  SchurFactoredFermionOperator<DomainWallFermionR::Impl_t>
    SchurPV(PVPeriodic,
	    PVDirichlet,
	    CG,Block);


  
  DomainDecomposedBoundaryTwoFlavourPseudoFermion<DomainWallFermionR::Impl_t> DBPFA(SchurDwf,CG,CG);
  ForceTest<GimplTypesR>(DBPFA,U,FilterDDHMC);

  DomainDecomposedBoundaryTwoFlavourRatioPseudoFermion<DomainWallFermionR::Impl_t> DBPFRA(SchurPV,SchurDwf,CG,CG);
  ForceTest<GimplTypesR>(DBPFRA,U,FilterDDHMC);

  DomainDecomposedBoundaryTwoFlavourBosonPseudoFermion<DomainWallFermionR::Impl_t> DBPFBA(SchurPV,CG,CG);
  ForceTest<GimplTypesR>(DBPFBA,U,FilterDDHMC);

  Grid_finalize();

}
