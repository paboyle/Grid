    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_gparity.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
using namespace Grid::QCD;

//typedef GparityDomainWallFermionD GparityDiracOp;
//typedef DomainWallFermionD StandardDiracOp;
//#define DOP_PARAMS

typedef GparityMobiusFermionD GparityDiracOp;
typedef MobiusFermionD StandardDiracOp;
#define DOP_PARAMS ,1.5, 0.5


typedef typename GparityDiracOp::FermionField GparityFermionField;
typedef typename GparityDiracOp::GaugeField GparityGaugeField;
typedef typename GparityFermionField::vector_type vComplexType;

typedef typename StandardDiracOp::FermionField StandardFermionField;
typedef typename StandardDiracOp::GaugeField StandardGaugeField;

enum{ same_vComplex = std::is_same<vComplexType, typename StandardFermionField::vector_type>::value };
static_assert(same_vComplex == 1, "Dirac Operators must have same underlying SIMD complex type");

int main (int argc, char ** argv)
{
  int nu = 0;

  Grid_init(&argc,&argv);

  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--Gparity-dir"){
      std::stringstream ss; ss << argv[i+1]; ss >> nu;
      std::cout << GridLogMessage << "Set Gparity direction to " << nu << std::endl;
    }
  }

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Testing Gparity Dirac operator                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplexType::Nsimd()<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using UNROLLED Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  const int Ls=4;
  //const int L =4;
  //std::vector<int> latt_2f(Nd,L);

  std::vector<int> latt_2f = GridDefaultLatt();
  std::vector<int> latt_1f(latt_2f); latt_1f[nu] = 2*latt_2f[nu];
  int L = latt_2f[nu];


  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexType::Nsimd());

  std::cout << GridLogMessage << "SIMD layout: ";
  for(int i=0;i<simd_layout.size();i++) std::cout << simd_layout[i] << " ";
  std::cout << std::endl;
  
  std::vector<int> mpi_layout  = GridDefaultMpi(); //node layout

  GridCartesian         * UGrid_1f   = SpaceTimeGrid::makeFourDimGrid(latt_1f, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_1f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_1f);
  GridCartesian         * FGrid_1f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_1f);
  GridRedBlackCartesian * FrbGrid_1f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_1f);


  GridCartesian         * UGrid_2f   = SpaceTimeGrid::makeFourDimGrid(latt_2f, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_2f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_2f);
  GridCartesian         * FGrid_2f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_2f);
  GridRedBlackCartesian * FrbGrid_2f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_2f);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5_2f(FGrid_2f);  RNG5_2f.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4_2f(UGrid_2f);  RNG4_2f.SeedFixedIntegers(seeds4);

  GparityGaugeField Umu_2f(UGrid_2f);
  SU3::HotConfiguration(RNG4_2f,Umu_2f);

  StandardFermionField    src   (FGrid_2f); 
  StandardFermionField    tmpsrc(FGrid_2f); 
  GparityFermionField      src_2f(FGrid_2f); 
  StandardFermionField    src_1f(FGrid_1f); 

  // Replicate fermion source
  random(RNG5_2f,src);
  PokeIndex<0>(src_2f,src,0);
  tmpsrc=src*2.0;
  PokeIndex<0>(src_2f,tmpsrc,1);

  StandardFermionField result_1f(FGrid_1f); result_1f=zero;
  StandardGaugeField Umu_1f(UGrid_1f); 
  Replicate(Umu_2f,Umu_1f);

  //Coordinate grid for reference
  LatticeInteger xcoor_1f(UGrid_1f);
  LatticeCoordinate(xcoor_1f,nu);

  //Copy-conjugate the gauge field
  //First C-shift the lattice by Lx/2
  {
    StandardGaugeField Umu_shift = conjugate( Cshift(Umu_1f,nu,L) );
    Umu_1f = where( xcoor_1f >= Integer(L), Umu_shift, Umu_1f );

    // hack test to check the same
    Replicate(Umu_2f,Umu_shift);
    Umu_shift=Umu_shift-Umu_1f;
    cout << GridLogMessage << "Umu diff " << norm2(Umu_shift)<<std::endl;

    //Make the gauge field antiperiodic in nu-direction
    decltype(PeekIndex<LorentzIndex>(Umu_1f,nu)) Unu(UGrid_1f);
    Unu = PeekIndex<LorentzIndex>(Umu_1f,nu);
    Unu = where(xcoor_1f == Integer(2*L-1), -Unu, Unu);
    PokeIndex<LorentzIndex>(Umu_1f,Unu,nu);
  }

  //Coordinate grid for reference
  LatticeInteger    xcoor_1f5(FGrid_1f);
  LatticeCoordinate(xcoor_1f5,1+nu);
  Replicate(src,src_1f);
  src_1f   = where( xcoor_1f5 >= Integer(L), 2.0*src_1f,src_1f );

  RealD mass=0.0;
  RealD M5=1.8;
  StandardDiracOp Ddwf(Umu_1f,*FGrid_1f,*FrbGrid_1f,*UGrid_1f,*UrbGrid_1f,mass,M5 DOP_PARAMS);

  StandardFermionField    src_o_1f(FrbGrid_1f);
  StandardFermionField result_o_1f(FrbGrid_1f);
  pickCheckerboard(Odd,src_o_1f,src_1f);
  result_o_1f=zero;

  SchurDiagMooeeOperator<StandardDiracOp,StandardFermionField> HermOpEO(Ddwf);
  ConjugateGradient<StandardFermionField> CG(1.0e-8,10000);
  CG(HermOpEO,src_o_1f,result_o_1f);
  
  //  const int nu = 3;
  std::vector<int> twists(Nd,0);
  twists[nu] = 1;
  GparityDiracOp::ImplParams params;
  params.twists = twists;
  GparityDiracOp GPDdwf(Umu_2f,*FGrid_2f,*FrbGrid_2f,*UGrid_2f,*UrbGrid_2f,mass,M5 DOP_PARAMS,params);

  for(int disp=-1;disp<=1;disp+=2)
  for(int mu=0;mu<5;mu++)
  { 
    GparityFermionField Dsrc_2f(FGrid_2f);

    StandardFermionField Dsrc_1f(FGrid_1f);
    StandardFermionField Dsrc_2freplica(FGrid_1f);
    StandardFermionField Dsrc_2freplica0(FGrid_1f);
    StandardFermionField Dsrc_2freplica1(FGrid_1f);

    if ( mu ==0 ) {
      std::cout << GridLogMessage<< " Cross checking entire hopping term"<<std::endl;
      GPDdwf.Dhop(src_2f,Dsrc_2f,DaggerNo);
        Ddwf.Dhop(src_1f,Dsrc_1f,DaggerNo);
    } else { 
      std::cout << GridLogMessage<< " Cross checking mu="<<mu<< " disp="<< disp<<std::endl;
      GPDdwf.DhopDir(src_2f,Dsrc_2f,mu,disp);
        Ddwf.DhopDir(src_1f,Dsrc_1f,mu,disp);
    }

    std::cout << GridLogMessage << "S norms "<< norm2(src_2f) << " " << norm2(src_1f)  <<std::endl;
    std::cout << GridLogMessage << "D norms "<< norm2(Dsrc_2f)<< " " << norm2(Dsrc_1f) <<std::endl;

    StandardFermionField Dsrc_2f0(FGrid_2f); Dsrc_2f0 = PeekIndex<0>(Dsrc_2f,0);
    StandardFermionField Dsrc_2f1(FGrid_2f); Dsrc_2f1 = PeekIndex<0>(Dsrc_2f,1);

    //    Dsrc_2f1 = Dsrc_2f1 - Dsrc_2f0;
    //    std::cout << GridLogMessage << " Cross check two halves " <<norm2(Dsrc_2f1)<<std::endl;

    Replicate(Dsrc_2f0,Dsrc_2freplica0);
    Replicate(Dsrc_2f1,Dsrc_2freplica1);

    Dsrc_2freplica = where( xcoor_1f5 >= Integer(L), Dsrc_2freplica1,Dsrc_2freplica0 );
    Dsrc_2freplica = Dsrc_2freplica - Dsrc_1f ;
    std::cout << GridLogMessage << " Cross check against doubled latt " <<norm2(Dsrc_2freplica)<<std::endl;

    //    std::cout << Dsrc_2f <<std::endl;

  }

  {
    GparityFermionField chi   (FGrid_2f); gaussian(RNG5_2f,chi);
    GparityFermionField phi   (FGrid_2f); gaussian(RNG5_2f,phi);
  
    GparityFermionField chi_e   (FrbGrid_2f);
    GparityFermionField chi_o   (FrbGrid_2f);
    
    GparityFermionField dchi_e  (FrbGrid_2f);
    GparityFermionField dchi_o  (FrbGrid_2f);
    
    GparityFermionField phi_e   (FrbGrid_2f);
    GparityFermionField phi_o   (FrbGrid_2f);
    
    GparityFermionField dphi_e  (FrbGrid_2f);
    GparityFermionField dphi_o  (FrbGrid_2f);

    pickCheckerboard(Even,chi_e,chi);
    pickCheckerboard(Odd ,chi_o,chi);
    pickCheckerboard(Even,phi_e,phi);
    pickCheckerboard(Odd ,phi_o,phi);

    GPDdwf.Meooe(chi_e,dchi_o);
    GPDdwf.Meooe(chi_o,dchi_e);
    GPDdwf.MeooeDag(phi_e,dphi_o);
    GPDdwf.MeooeDag(phi_o,dphi_e);
    
    ComplexD pDce = innerProduct(phi_e,dchi_e);
    ComplexD pDco = innerProduct(phi_o,dchi_o);
    ComplexD cDpe = innerProduct(chi_e,dphi_e);
    ComplexD cDpo = innerProduct(chi_o,dphi_o);
    
    std::cout<<GridLogMessage <<"e "<<pDce<<" "<<cDpe <<std::endl;
    std::cout<<GridLogMessage <<"o "<<pDco<<" "<<cDpo <<std::endl;
    
    std::cout<<GridLogMessage <<"pDce - conj(cDpo) "<< pDce-conj(cDpo) <<std::endl;
    std::cout<<GridLogMessage <<"pDco - conj(cDpe) "<< pDco-conj(cDpe) <<std::endl;

  }

  GparityFermionField result_2f(FGrid_2f); result_2f=zero;
  GparityFermionField    src_o_2f(FrbGrid_2f);
  GparityFermionField result_o_2f(FrbGrid_2f);
  pickCheckerboard(Odd,src_o_2f,src_2f);
  result_o_2f=zero;

  ConjugateGradient<GparityFermionField> CG2f(1.0e-8,10000);
  SchurDiagMooeeOperator<GparityDiracOp,GparityFermionField> HermOpEO2f(GPDdwf);
  CG2f(HermOpEO2f,src_o_2f,result_o_2f);

  std::cout << "2f cb "<<result_o_2f.checkerboard<<std::endl;
  std::cout << "1f cb "<<result_o_1f.checkerboard<<std::endl;

  std::cout << " result norms " <<norm2(result_o_2f)<<" " <<norm2(result_o_1f)<<std::endl;

  StandardFermionField    res0o  (FrbGrid_2f); 
  StandardFermionField    res1o  (FrbGrid_2f); 
  StandardFermionField    res0  (FGrid_2f); 
  StandardFermionField    res1  (FGrid_2f); 

  res0=zero;
  res1=zero;

  res0o = PeekIndex<0>(result_o_2f,0);
  res1o = PeekIndex<0>(result_o_2f,1);

  std::cout << "res cb "<<res0o.checkerboard<<std::endl;
  std::cout << "res cb "<<res1o.checkerboard<<std::endl;

  setCheckerboard(res0,res0o);
  setCheckerboard(res1,res1o);

  StandardFermionField replica (FGrid_1f);
  StandardFermionField replica0(FGrid_1f);
  StandardFermionField replica1(FGrid_1f);
  Replicate(res0,replica0);
  Replicate(res1,replica1);

  replica = where( xcoor_1f5 >= Integer(L), replica1,replica0 );

  replica0 = zero;
  setCheckerboard(replica0,result_o_1f);

  std::cout << "Norm2 solutions is " <<norm2(replica)<<" "<< norm2(replica0)<<std::endl;

  replica = replica - replica0;
  
  std::cout << "Norm2 of difference in solutions is " <<norm2(replica)<<std::endl;
  

  Grid_finalize();
}
