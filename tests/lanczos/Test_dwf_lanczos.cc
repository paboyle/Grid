    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_lanczos.cc

    Copyright (C) 2015

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

using namespace std;
using namespace Grid;
 ;

template<typename Action>
struct Setup{};

template<>
struct Setup<GparityMobiusFermionF>{
  static GparityMobiusFermionF* getAction(LatticeGaugeFieldF &Umu,
					  GridCartesian* FGrid, GridRedBlackCartesian* FrbGrid, GridCartesian* UGrid, GridRedBlackCartesian* UrbGrid){
    RealD mass=0.00054;
    RealD M5=1.8;
    RealD mob_b=1.5;
    GparityMobiusFermionD ::ImplParams params;
    std::vector<int> twists({1,1,1,0});
    params.twists = twists;
    return new GparityMobiusFermionF(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_b-1.,params);
  }
};

template<>
struct Setup<DomainWallFermionF>{
  static DomainWallFermionF* getAction(LatticeGaugeFieldF &Umu,
struct Setup<DomainWallFermionD>{
  static DomainWallFermionD* getAction(LatticeGaugeField &Umu,
					  GridCartesian* FGrid, GridRedBlackCartesian* FrbGrid, GridCartesian* UGrid, GridRedBlackCartesian* UrbGrid){
    RealD mass=0.00054;
    RealD M5=1.8;
    return new DomainWallFermionF(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  }
};

template<>
struct Setup<MobiusFermionF>{
  static MobiusFermionF* getAction(LatticeGaugeFieldF &Umu,
					  GridCartesian* FGrid, GridRedBlackCartesian* FrbGrid, GridCartesian* UGrid, GridRedBlackCartesian* UrbGrid){
    RealD mass=0.00054;
    RealD M5=1.8;
    RealD mob_b=1.5;
    std::vector<Complex> boundary = {1,1,1,-1};
    MobiusFermionF::ImplParams Params(boundary);

  std::cout << GridLogMessage << "mass "<<mass<<std::endl;
  std::cout << GridLogMessage << "M5 "<<M5<<std::endl;
  std::cout << GridLogMessage << "mob_b "<<mob_b<<std::endl;
    return new MobiusFermionF(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_b-1.,Params);
  }
};



template<typename Action>
void run(){
  typedef typename Action::FermionField FermionField;
  const int Ls=12;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
//  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n",UGrid,UrbGrid,FGrid,FrbGrid);
  
  GridCartesian* UGridF = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian* FGridF = SpaceTimeGrid::makeFiveDimGrid(Ls, UGridF);
  GridRedBlackCartesian* FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGridF);


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGridF);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGridF);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5rb(FrbGridF);  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid); 
//  SU<Nc>::HotConfiguration(RNG4, Umu);
  FieldMetaData header;
  std::string file("./config");

//  int precision32 = 0;
//  int tworow      = 0;
//  NerscIO::writeConfiguration(Umu,file,tworow,precision32);
  NerscIO::readConfiguration(Umu,header,file);

  LatticeGaugeFieldF UmuF(UGridF); 
  precisionChange(UmuF, Umu);

  Action *action = Setup<Action>::getAction(UmuF,FGridF,FrbGridF,UGridF,UrbGridF);
 
  //MdagMLinearOperator<Action,FermionField> HermOp(Ddwf);
//  SchurDiagTwoOperator<Action,FermionField> HermOp(*action);
  SchurDiagOneOperator<Action,FermionField> HermOp(*action);

  const int Nstop = 150;
  const int Nk = 160;
  const int Np = 40;
  const int Nm = Nk+Np;
  const int MaxIt= 10000;
  RealD resid = 1.0e-6;
  std::cout << GridLogMessage << "Nstop "<<Nstop<<std::endl;
  std::cout << GridLogMessage << "Nk "<<Nk<<std::endl;
  std::cout << GridLogMessage << "Np "<<Np<<std::endl;
  std::cout << GridLogMessage << "resid "<<resid<<std::endl;

  std::vector<double> Coeffs { 0.,-1.};
  Polynomial<FermionField> PolyX(Coeffs);
  Chebyshev<FermionField> Cheby(0.0000006,5.5,4001);
  std::cout << GridLogMessage << "Cheby(0.0000006,5.5,4001) "<<std::endl;

  FunctionHermOp<FermionField> OpCheby(Cheby,HermOp);
  PlainHermOp<FermionField> Op     (HermOp);

  ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby,Op,Nstop,Nk,Nm,resid,MaxIt);
 
  std::vector<RealD>          eval(Nm);
  FermionField    src(FrbGridF); 
  gaussian(RNG5rb,src);
  std::vector<FermionField> evec(Nm,FrbGridF);
  for(int i=0;i<1;i++){
    std::cout << GridLogMessage <<i<<" / "<< Nm<< " grid pointer "<<evec[i].Grid()<<std::endl;
  };

  int Nconv;
  IRL.calc(eval,evec,src,Nconv);

  delete action;
}
  
int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::string action = "Mobius";
  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "-action"){
      action = argv[i+1];
    }
  }

  if(action == "GparityMobius"){
    run<GparityMobiusFermionF>();
  }else if(action == "DWF"){
    run<DomainWallFermionF>();
  }else if(action == "Mobius"){
    run<MobiusFermionF>();
  }else{
    std::cout << "Unknown action" << std::endl;
    exit(1);
  }
  
  Grid_finalize();
}
