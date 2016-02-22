    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_even_odd.cc

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

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          pRNG(&Grid);
  //  std::vector<int> seeds({1,2,3,4});
  //  pRNG.SeedFixedIntegers(seeds);
  pRNG.SeedRandomDevice();

  LatticeFermion src   (&Grid); random(pRNG,src);
  LatticeFermion phi   (&Grid); random(pRNG,phi);
  LatticeFermion chi   (&Grid); random(pRNG,chi);
  LatticeFermion result(&Grid); result=zero;
  LatticeFermion    ref(&Grid);    ref=zero;
  LatticeFermion    tmp(&Grid);    tmp=zero;
  LatticeFermion    err(&Grid);    tmp=zero;
  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);
  std::vector<LatticeColourMatrix> U(4,&Grid);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  // Only one non-zero (y)
  Umu=zero;
  for(int nn=0;nn<Nd;nn++){
    random(pRNG,U[nn]);
    std::cout<<GridLogMessage<<"U[nn]"<<norm2(U[nn])<<std::endl;
    PokeIndex<LorentzIndex>(Umu,U[nn],nn);
    std::cout<<GridLogMessage<<"Umu"<<norm2(Umu)<<std::endl;
  }

  RealD mass=0.1;
  RealD mu  = 0.1;
  WilsonTMFermionR Dw(Umu,Grid,RBGrid,mass,mu);

  LatticeFermion src_e   (&RBGrid);
  LatticeFermion src_o   (&RBGrid);
  LatticeFermion r_e   (&RBGrid);
  LatticeFermion r_o   (&RBGrid);
  LatticeFermion r_eo  (&Grid);

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that Deo + Doe = Dunprec "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  Dw.Meooe(src_e,r_o);  std::cout<<GridLogMessage<<"Applied Meo"<<std::endl;
  Dw.Meooe(src_o,r_e);  std::cout<<GridLogMessage<<"Applied Moe"<<std::endl;
  Dw.Dhop (src,ref,DaggerNo);

  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  err= ref - r_eo;
  std::cout<<GridLogMessage << "EO norm diff   "<< norm2(err)<< " "<<norm2(ref)<< " " << norm2(r_eo) <<std::endl;

  LatticeComplex cerr(&Grid);
  cerr = localInnerProduct(err,err);

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test Ddagger is the dagger of D by requiring                "<<std::endl;
  std::cout<<GridLogMessage<<"=  < phi | Deo | chi > * = < chi | Deo^dag| phi>  "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  
  LatticeFermion chi_e   (&RBGrid);
  LatticeFermion chi_o   (&RBGrid);

  LatticeFermion dchi_e  (&RBGrid);
  LatticeFermion dchi_o  (&RBGrid);

  LatticeFermion phi_e   (&RBGrid);
  LatticeFermion phi_o   (&RBGrid);

  LatticeFermion dphi_e  (&RBGrid);
  LatticeFermion dphi_o  (&RBGrid);


  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);

  Dw.Meooe(chi_e,dchi_o);
  Dw.Meooe(chi_o,dchi_e);
  Dw.MeooeDag(phi_e,dphi_o);
  Dw.MeooeDag(phi_o,dphi_e);

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

  Dw.Mooee(chi_e,src_e);
  Dw.MooeeInv(src_e,phi_e);

  Dw.Mooee(chi_o,src_o);
  Dw.MooeeInv(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MeeInvDag MeeDag = 1                                   "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Dw.MooeeDag(chi_e,src_e);
  Dw.MooeeInvDag(src_e,phi_e);

  Dw.MooeeDag(chi_o,src_o);
  Dw.MooeeInvDag(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MpcDagMpc is Hermitian              "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  
  random(pRNG,phi);
  random(pRNG,chi);
  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);
  RealD t1,t2;

  SchurDiagMooeeOperator<WilsonTMFermionR,LatticeFermion> HermOpEO(Dw);
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
