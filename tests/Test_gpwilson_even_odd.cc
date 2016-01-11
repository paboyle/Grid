    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_gpwilson_even_odd.cc

    Copyright (C) 2015

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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

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

  typedef typename GparityWilsonFermionR::FermionField FermionField;

  FermionField src   (&Grid); random(pRNG,src);
  FermionField phi   (&Grid); random(pRNG,phi);
  FermionField chi   (&Grid); random(pRNG,chi);
  FermionField result(&Grid); result=zero;
  FermionField    ref(&Grid);    ref=zero;
  FermionField    tmp(&Grid);    tmp=zero;
  FermionField    err(&Grid);    tmp=zero;
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

  GparityWilsonFermionR::ImplParams params;
  std::vector<int> twists(Nd,0);  twists[1] = 1;
  params.twists = twists;
  GparityWilsonFermionR Dw(Umu,Grid,RBGrid,mass,params);

  FermionField src_e   (&RBGrid);
  FermionField src_o   (&RBGrid);
  FermionField r_e   (&RBGrid);
  FermionField r_o   (&RBGrid);
  FermionField r_eo  (&Grid);

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
  
  FermionField chi_e   (&RBGrid);
  FermionField chi_o   (&RBGrid);

  FermionField dchi_e  (&RBGrid);
  FermionField dchi_o  (&RBGrid);

  FermionField phi_e   (&RBGrid);
  FermionField phi_o   (&RBGrid);

  FermionField dphi_e  (&RBGrid);
  FermionField dphi_o  (&RBGrid);


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

  SchurDiagMooeeOperator<GparityWilsonFermionR,FermionField> HermOpEO(Dw);
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
