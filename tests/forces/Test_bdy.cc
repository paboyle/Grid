/*

  2f Full det MdagM 10^6 force ~ 1.3e7
rid : Message : 1767.283471 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 1767.283476 s : S1 : 1.52885e+09
Grid : Message : 1767.283480 s : S2 : 1.52886e+09
Grid : Message : 1767.283482 s : dS : 8877.34
Grid : Message : 1767.283483 s : dSpred : 8877.7
Grid : Message : 1767.283484 s : diff : -0.360484
Grid : Message : 1767.283485 s : *********************************************************

  2f Full det MpcdagMpc 10^6 force ~ 1.8e6
Grid : Message : 2399.576962 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 2399.576968 s : S1 : 1.52885e+09
Grid : Message : 2399.576972 s : S2 : 1.52886e+09
Grid : Message : 2399.576974 s : dS : 9728.49
Grid : Message : 2399.576975 s : dSpred : 9726.58
Grid : Message : 2399.576976 s : diff : 1.90683
Grid : Message : 2399.576977 s : *********************************************************

  2f bdy MdagM 1500 force Force ~ 2800
Grid : Message : 4622.385061 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 4622.385067 s : S1 : 1.52885e+09
Grid : Message : 4622.385071 s : S2 : 1.52885e+09
Grid : Message : 4622.385072 s : dS : 25.4944
Grid : Message : 4622.385073 s : dSpred : 25.4672
Grid : Message : 4622.385074 s : diff : 0.0271414
Grid : Message : 4622.385075 s : *********************************************************

  2f bdy MpcdagMpc 10^6 force   Force ~ 2200
Grid : Message : 4622.385061 s : +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Grid : Message : 4622.385067 s : S1 : 1.52885e+09
Grid : Message : 4622.385071 s : S2 : 1.52885e+09
Grid : Message : 4622.385072 s : dS : 25.4944
Grid : Message : 4622.385073 s : dSpred : 25.4672
Grid : Message : 4622.385074 s : diff : 0.0271414
Grid : Message : 4622.385075 s : *********************************************************
  
  1f Bdy Det
  Optimisation log:  looser rational AND MD tolerances sloppy
MobiusForce.221179 -- same as HMC. dS is mispredicted Forece  ~2.8
Grid : Message : 6582.258991 s : dS : 0.024478
Grid : Message : 6582.258992 s : dSpred : 0.00791876
Grid : Message : 6582.258994 s : diff : 0.0165592

MobiusForce.221193 -- tight rational AND MD tolerances to 1e-8 ~ 2.8 same
Grid : Message : 1964.939209 s : S1 : 7.64404e+08
Grid : Message : 1964.939213 s : S2 : 7.64404e+08
Grid : Message : 1964.939215 s : dS : -0.00775838 <--- too loose even on action
Grid : Message : 1964.939216 s : dSpred : -0.00416793 
Grid : Message : 1964.939217 s : diff : -0.00359045

MobiusForce.221394 -- looser rational, MD tol 1e-8 ~ 2.8 same
Grid : Message : 1198.346720 s : S1 : 764404649.48886
Grid : Message : 1198.346760 s : S2 : 764404649.5133
Grid : Message : 1198.346780 s : dS : 0.024440884590149
Grid : Message : 1198.346800 s : dSpred : 0.0079145154465184
Grid : Message : 1198.346810 s : diff : 0.016526369143631

MobiusForce.221394 -- tight rational, MD tol sloppy Force ~ 2.8
Grid : Message : 2376.921950 s : S1 : 764404436.44069
Grid : Message : 2376.921954 s : S2 : 764404436.43299
Grid : Message : 2376.921956 s : dS : -0.0076971054077148
Grid : Message : 2376.921958 s : dSpred : -0.0041610472282526
Grid : Message : 2376.921959 s : diff : -0.0035360581794623

*/

//
/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_double_ratio.cc

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

using namespace std;
using namespace Grid;

typedef MobiusFermionD FermionAction;
typedef WilsonImplD FimplD;
typedef WilsonImplD FermionImplPolicy;

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

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Refresh "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  Gimpl::generate_momenta(P,sRNG,RNG4);
  Filter.applyFilter(P);

#if  0
  FieldMetaData header;
  std::string file("./ckpoint_lat.2000");
  NerscIO::readConfiguration(U,header,file);
#else
  U = 1.0;
#endif
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
 
  ////////////////////////////////////////////////////////////////
  // Domain decomposed operator
  ////////////////////////////////////////////////////////////////
  Coordinate CommDim(Nd);
  for(int d=0;d<Nd;d++) CommDim[d]= (mpi_layout[d]/shm[d])>1 ? 1 : 0;

  Coordinate NonDirichlet(Nd+1,0);
  Coordinate Dirichlet(Nd+1,0);
  Dirichlet[1] = CommDim[0]*latt_size[0]/mpi_layout[0] * shm[0];
  Dirichlet[2] = CommDim[1]*latt_size[1]/mpi_layout[1] * shm[1];
  Dirichlet[3] = CommDim[2]*latt_size[2]/mpi_layout[2] * shm[2];
  Dirichlet[4] = CommDim[3]*latt_size[3]/mpi_layout[3] * shm[3];

  Coordinate Block4(Nd);
  Block4[0] = Dirichlet[1];
  Block4[1] = Dirichlet[2];
  Block4[2] = Dirichlet[3];
  Block4[3] = Dirichlet[4];

  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);
  FermionAction::ImplParams ParamsDir(boundary);
  Params.dirichlet=NonDirichlet;
  ParamsDir.dirichlet=Dirichlet;
  ParamsDir.partialDirichlet=1;

  ///////////////////// Gauge Field and Gauge Forces ////////////////////////////
  LatticeGaugeField U(UGrid);

  RealD beta=6.0;
  WilsonGaugeActionR PlaqAction(beta);
  IwasakiGaugeActionR RectAction(beta);

  MomentumFilterNone<LatticeGaugeField> FilterNone;
  ForceTest<GimplTypesR>(PlaqAction,U,FilterNone);
  ForceTest<GimplTypesR>(RectAction,U,FilterNone);

  ////////////////////////////////////
  // Action
  ////////////////////////////////////
  RealD mass=0.00078; 
  RealD pvmass=1.0; 
  RealD M5=1.8; 
  RealD b=1.5;
  RealD c=0.5;
  
  // Double versions
  FermionAction DdwfPeriodic(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,Params);
  FermionAction PVPeriodic  (U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pvmass,M5,b,c,Params);
  FermionAction DdwfDirichlet(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,ParamsDir);

  double StoppingCondition = 1.0e-8;
  double MaxCGIterations = 50000;
  ConjugateGradient<LatticeFermion>  CG(StoppingCondition,MaxCGIterations);
  
  //////////////////// Two Flavour Determinant Ratio ///////////////////////////////
  TwoFlavourRatioPseudoFermionAction<FimplD> Nf2(PVPeriodic, DdwfPeriodic,CG,CG);
  //  ForceTest<GimplTypesR>(Nf2,U,FilterNone);

  //////////////////// Two Flavour Determinant force test Even Odd ///////////////////////////////
  TwoFlavourEvenOddRatioPseudoFermionAction<FimplD> Nf2eo(PVPeriodic, DdwfPeriodic,CG,CG);
  //  ForceTest<GimplTypesR>(Nf2eo,U,FilterNone);

  //////////////////// Domain forces ////////////////////
  int Width=4;
  DDHMCFilter<WilsonImplD::Field> DDHMCFilter(Block4,Width);
  
  //////////////////// Two flavour boundary det  ////////////////////
  TwoFlavourRatioPseudoFermionAction<FimplD> BdyNf2(DdwfDirichlet, DdwfPeriodic,CG,CG);
  //  ForceTest<GimplTypesR>(BdyNf2,U,DDHMCFilter);

  //////////////////// Two flavour eo boundary det  ////////////////////
  TwoFlavourEvenOddRatioPseudoFermionAction<FimplD> BdyNf2eo(DdwfDirichlet, DdwfPeriodic,CG,CG);
  //  ForceTest<GimplTypesR>(BdyNf2eo,U,DDHMCFilter);

  //////////////////// One flavour boundary det  ////////////////////
  OneFlavourRationalParams OFRp; // Up/down
  OFRp.lo       = 4.0e-5;
  OFRp.hi       = 90.0;
  OFRp.MaxIter  = 60000;
  OFRp.tolerance= 1.0e-8;
  OFRp.mdtolerance= 1.0e-6;
  OFRp.degree   = 18;
  OFRp.precision= 80;
  OFRp.BoundsCheckFreq=0;
  std::vector<RealD> ActionTolByPole({
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8
    });
  std::vector<RealD> MDTolByPole({
      1.0e-6,3.0e-7,1.0e-7,1.0e-7,  // Orig sloppy
      //      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8,
      1.0e-8,1.0e-8,1.0e-8,1.0e-8
    });
  OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> BdySqrt(DdwfDirichlet,DdwfPeriodic,OFRp);
  ForceTest<GimplTypesR>(BdySqrt,U,DDHMCFilter);

  Grid_finalize();
}
