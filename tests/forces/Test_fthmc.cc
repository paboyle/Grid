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


template<class Gimpl>
void ForceTest(Action<LatticeGaugeField> &action,SmearedConfigurationMasked<PeriodicGimplR> & smU,MomentumFilterBase<LatticeGaugeField> &Filter)
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
  
  RealD eps=0.005;

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Refresh "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  Gimpl::generate_momenta(P,sRNG,RNG4);
  Filter.applyFilter(P);

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
  Filter.applyFilter(UdSdU);

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
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);

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

  
  RealD beta=6.0;
  WilsonGaugeActionR  PlaqAction(beta);
  IwasakiGaugeActionR RectAction(beta);


  ////////////////////////////////////////////////
  // Plaquette only FTHMC smearer
  ////////////////////////////////////////////////
  double rho = 0.1;
  Smear_Stout<PeriodicGimplR> Smearer(rho);
  SmearedConfigurationMasked<PeriodicGimplR> SmartConfig(UGrid,2*Nd,Smearer,true);

  JacobianAction<PeriodicGimplR> Jacobian(&SmartConfig);
  
  ////////////////////////////////////////////////
  // Run some tests
  ////////////////////////////////////////////////
  MomentumFilterNone<LatticeGaugeField> FilterNone;
  SmartConfig.set_Field(U);
  ForceTest<GimplTypesR>(PlaqAction,SmartConfig,FilterNone);
  SmartConfig.set_Field(U);
  ForceTest<GimplTypesR>(RectAction,SmartConfig,FilterNone);
  SmartConfig.set_Field(U);
  ForceTest<GimplTypesR>(Jacobian,SmartConfig,FilterNone);

  Grid_finalize();
}
