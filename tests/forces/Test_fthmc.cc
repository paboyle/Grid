/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_fthmc.cc

    Copyright (C) 2022

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
#include <Grid/Grid.h>
#include <Grid/qcd/smearing/GaugeConfigurationMasked.h>
#include <Grid/qcd/smearing/JacobianAction.h>

using namespace std;
using namespace Grid;

typedef MobiusFermionD FermionAction;
typedef WilsonImplD FimplD;
typedef WilsonImplD FermionImplPolicy;

template<class Gimpl>
void ForceTest(Action<LatticeGaugeField> &action,ConfigurationBase<LatticeGaugeField> & smU,MomentumFilterBase<LatticeGaugeField> &Filter)
{
  LatticeGaugeField U = smU.get_U(false); // unsmeared config
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
  
  RealD eps=0.01;

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Refresh "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  Gimpl::generate_momenta(P,sRNG,RNG4);
  //  Filter.applyFilter(P);

  action.refresh(smU,sRNG,RNG4);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Action "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

  RealD S1 = action.S(smU);

  Gimpl::update_field(P,U,eps);
  smU.set_Field(U);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Derivative "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  action.deriv(smU,UdSdU);
  UdSdU = Ta(UdSdU);
  //  Filter.applyFilter(UdSdU);

  DumpSliceNorm("Force",UdSdU,Nd-1);
  
  Gimpl::update_field(P,U,eps);
  smU.set_Field(U);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Action "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  RealD S2 = action.S(smU);

  // Use the derivative
  LatticeComplex dS(UGrid); dS = Zero();
  for(int mu=0;mu<Nd;mu++){
    auto UdSdUmu = PeekIndex<LorentzIndex>(UdSdU,mu);
    Pmu= PeekIndex<LorentzIndex>(P,mu);
    dS = dS - trace(Pmu*UdSdUmu)*eps*2.0*HMC_MOMENTUM_DENOMINATOR;
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
  //  assert(diff<1.0);
  std::cout<< GridLogMessage << "Done" <<std::endl;
  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout << std::setprecision(14);
  Coordinate latt_size   = GridDefaultLatt();
  Coordinate mpi_layout  = GridDefaultMpi();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate shm;
  GlobalSharedMemory::GetShmDims(mpi_layout,shm);

  const int Ls=12;
  const int Nt = latt_size[3];
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  
  ///////////////////// Gauge Field and Gauge Forces ////////////////////////////
  LatticeGaugeField U(UGrid);

#if  0
  FieldMetaData header;
  std::string file("./ckpoint_lat.2000");
  NerscIO::readConfiguration(U,header,file);
#else
  std::vector<int> seeds({1,2,3,4,5,6,7,8});
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds);
  SU<Nc>::HotConfiguration(RNG4,U);
#endif

  
  WilsonGaugeActionR  PlaqAction(6.0);
  IwasakiGaugeActionR RectAction(2.13);
  PlaqAction.is_smeared = true;  
  RectAction.is_smeared = true;  

  ////////////////////////////////////
  // Fermion Action
  ////////////////////////////////////
  RealD mass=0.01; 
  RealD pvmass=1.0; 
  RealD M5=1.8; 
  RealD b=1.5;
  RealD c=0.5;
  
  // Double versions
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);
  FermionAction DdwfPeriodic(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,Params);
  FermionAction PVPeriodic  (U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pvmass,M5,b,c,Params);

  double StoppingCondition = 1.0e-8;
  double MaxCGIterations = 50000;
  ConjugateGradient<LatticeFermion>  CG(StoppingCondition,MaxCGIterations);

  TwoFlavourRatioPseudoFermionAction<FimplD> Nf2(PVPeriodic, DdwfPeriodic,CG,CG);
  Nf2.is_smeared = true;  
  
  ////////////////////////////////////////////////
  // Plaquette only FTHMC smearer
  ////////////////////////////////////////////////
  double rho = 0.1;
  Smear_Stout<PeriodicGimplR> Smearer(rho);
  SmearedConfigurationMasked<PeriodicGimplR> SmartConfig(UGrid,2*Nd,Smearer);
  SmearedConfiguration<PeriodicGimplR> StoutConfig(UGrid,1,Smearer);

  JacobianAction<PeriodicGimplR> Jacobian(&SmartConfig);
  
  ////////////////////////////////////////////////
  // Run some tests
  ////////////////////////////////////////////////
  MomentumFilterNone<LatticeGaugeField> FilterNone;

  std::cout << " *********  FIELD TRANSFORM SMEARING ***** "<<std::endl;

  SmartConfig.set_Field(U);
  ForceTest<GimplTypesR>(PlaqAction,SmartConfig,FilterNone);

  SmartConfig.set_Field(U);
  ForceTest<GimplTypesR>(RectAction,SmartConfig,FilterNone);

  SmartConfig.set_Field(U);
  ForceTest<GimplTypesR>(Jacobian,SmartConfig,FilterNone);

  SmartConfig.set_Field(U);
  ForceTest<GimplTypesR>(Nf2,SmartConfig,FilterNone);

  std::cout << " *********    STOUT SMEARING ***** "<<std::endl;

  StoutConfig.set_Field(U);
  ForceTest<GimplTypesR>(PlaqAction,StoutConfig,FilterNone);

  StoutConfig.set_Field(U);
  ForceTest<GimplTypesR>(RectAction,StoutConfig,FilterNone);
  
  StoutConfig.set_Field(U);
  ForceTest<GimplTypesR>(Nf2,StoutConfig,FilterNone);
  

  Grid_finalize();
}
