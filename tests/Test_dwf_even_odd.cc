    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_even_odd.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class d>
struct scal {
  d internal;
};

  Gamma::GammaMatrix Gmu [] = {
    Gamma::GammaX,
    Gamma::GammaY,
    Gamma::GammaZ,
    Gamma::GammaT
  };


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;


  const int Ls=8;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);

  LatticeFermion src   (FGrid); random(RNG5,src);
  LatticeFermion phi   (FGrid); random(RNG5,phi);
  LatticeFermion chi   (FGrid); random(RNG5,chi);
  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid);    ref=zero;
  LatticeFermion    tmp(FGrid);    tmp=zero;
  LatticeFermion    err(FGrid);    tmp=zero;
  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);
  std::vector<LatticeColourMatrix> U(4,UGrid);

  // Only one non-zero (y)
  Umu=zero;
  for(int nn=0;nn<Nd;nn++){
    random(RNG4,U[nn]);
    if ( nn>0 ) 
      U[nn]=zero;
    PokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }

  RealD mass=0.1;
  RealD M5  =1.8;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  LatticeFermion src_e (FrbGrid);
  LatticeFermion src_o (FrbGrid);
  LatticeFermion r_e   (FrbGrid);
  LatticeFermion r_o   (FrbGrid);
  LatticeFermion r_eo  (FGrid);
  LatticeFermion r_eeoo(FGrid);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that Meo + Moe + Moo + Mee = Munprec "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  Ddwf.Meooe(src_e,r_o);  std::cout<<GridLogMessage<<"Applied Meo"<<std::endl;
  Ddwf.Meooe(src_o,r_e);  std::cout<<GridLogMessage<<"Applied Moe"<<std::endl;
  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  Ddwf.Mooee(src_e,r_e);  std::cout<<GridLogMessage<<"Applied Mee"<<std::endl;
  Ddwf.Mooee(src_o,r_o);  std::cout<<GridLogMessage<<"Applied Moo"<<std::endl;
  setCheckerboard(r_eeoo,r_e);
  setCheckerboard(r_eeoo,r_o);

  r_eo=r_eo+r_eeoo;
  Ddwf.M(src,ref);  

  //  std::cout<<GridLogMessage << r_eo<<std::endl;
  //  std::cout<<GridLogMessage << ref <<std::endl;

  err= ref - r_eo;
  std::cout<<GridLogMessage << "EO norm diff   "<< norm2(err)<< " "<<norm2(ref)<< " " << norm2(r_eo) <<std::endl;
    
  LatticeComplex cerr(FGrid);
  cerr = localInnerProduct(err,err);
  //  std::cout<<GridLogMessage << cerr<<std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test Ddagger is the dagger of D by requiring                "<<std::endl;
  std::cout<<GridLogMessage<<"=  < phi | Deo | chi > * = < chi | Deo^dag| phi>  "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  
  LatticeFermion chi_e   (FrbGrid);
  LatticeFermion chi_o   (FrbGrid);

  LatticeFermion dchi_e  (FrbGrid);
  LatticeFermion dchi_o  (FrbGrid);

  LatticeFermion phi_e   (FrbGrid);
  LatticeFermion phi_o   (FrbGrid);

  LatticeFermion dphi_e  (FrbGrid);
  LatticeFermion dphi_o  (FrbGrid);


  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);

  Ddwf.Meooe(chi_e,dchi_o);
  Ddwf.Meooe(chi_o,dchi_e);
  Ddwf.MeooeDag(phi_e,dphi_o);
  Ddwf.MeooeDag(phi_o,dphi_e);

  ComplexD pDce = innerProduct(phi_e,dchi_e);
  ComplexD pDco = innerProduct(phi_o,dchi_o);
  ComplexD cDpe = innerProduct(chi_e,dphi_e);
  ComplexD cDpo = innerProduct(chi_o,dphi_o);

  std::cout<<GridLogMessage <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout<<GridLogMessage <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout<<GridLogMessage <<"pDce - conj(cDpo) "<< pDce-conj(cDpo) <<std::endl;
  std::cout<<GridLogMessage <<"pDco - conj(cDpe) "<< pDco-conj(cDpe) <<std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MeeInv Mee = 1                                         "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Ddwf.Mooee(chi_e,src_e);
  Ddwf.MooeeInv(src_e,phi_e);

  Ddwf.Mooee(chi_o,src_o);
  Ddwf.MooeeInv(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MeeInvDag MeeDag = 1                                   "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Ddwf.MooeeDag(chi_e,src_e);
  Ddwf.MooeeInvDag(src_e,phi_e);

  Ddwf.MooeeDag(chi_o,src_o);
  Ddwf.MooeeInvDag(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MpcDagMpc is Hermitian              "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  
  random(RNG5,phi);
  random(RNG5,chi);
  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);
  RealD t1,t2;


  SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> HermOpEO(Ddwf);
  HermOpEO.MpcDagMpc(chi_e,dchi_e,t1,t2);
  HermOpEO.MpcDagMpc(chi_o,dchi_o,t1,t2);

  HermOpEO.MpcDagMpc(phi_e,dphi_e,t1,t2);
  HermOpEO.MpcDagMpc(phi_o,dphi_o,t1,t2);

  pDce = innerProduct(phi_e,dchi_e);
  pDco = innerProduct(phi_o,dchi_o);
  cDpe = innerProduct(chi_e,dphi_e);
  cDpo = innerProduct(chi_o,dphi_o);

  std::cout<<GridLogMessage <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout<<GridLogMessage <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout<<GridLogMessage <<"pDce - conj(cDpo) "<< pDco-conj(cDpo) <<std::endl;
  std::cout<<GridLogMessage <<"pDco - conj(cDpe) "<< pDce-conj(cDpe) <<std::endl;
  
  Grid_finalize();
}
