    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_wilson.cc

    Copyright (C) 2015

    Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
  std::cout<<GridLogMessage << "Grid floating point word size is REALF"<< sizeof(RealF)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALD"<< sizeof(RealD)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REAL"<< sizeof(Real)<<std::endl;

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);
  //  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9});

  typedef typename WilsonCloverFermionR::FermionField FermionField; 
  typename WilsonCloverFermionR::ImplParams params; 

  FermionField src   (&Grid); random(pRNG,src);
  FermionField result(&Grid); result=zero;
  FermionField    ref(&Grid);    ref=zero;
  FermionField    tmp(&Grid);    tmp=zero;
  FermionField    err(&Grid);    tmp=zero;
  FermionField phi   (&Grid); random(pRNG,phi);
  FermionField chi   (&Grid); random(pRNG,chi);
  LatticeGaugeField Umu(&Grid); SU3::HotConfiguration(pRNG,Umu);
  std::vector<LatticeColourMatrix> U(4,&Grid);


  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  // Only one non-zero (y)
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  /* Debug force unit
    U[mu] = 1.0;
    PokeIndex<LorentzIndex>(Umu,U[mu],mu);
  */
  }

  ref = zero;

  RealD mass=0.1;
  RealD csw = 1.0;

  { // Simple clover implementation
  
    //    ref = ref + mass * src;
  }

  WilsonCloverFermionR Dwc(Umu,Grid,RBGrid,mass,csw,params);
  

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing Dhop against cshift implementation         "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  std::cout<<GridLogMessage << "Calling Dwc"<<std::endl;
  int ncall=1000;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Dwc.Dhop(src,result,0);
  }
  double t1=usecond();
  double t2;
  double flops=(16*(3*(6+8+8)) + 15*3*2)*volume*ncall; // == 66*16 +  == 1146
  
  std::cout<<GridLogMessage << "Called Dwc"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;

  err = ref-result; 
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that Deo + Doe = Dunprec "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  FermionField src_e   (&RBGrid);
  FermionField src_o   (&RBGrid);
  FermionField r_e   (&RBGrid);
  FermionField r_o   (&RBGrid);
  FermionField r_eo  (&Grid);
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  Dwc.Meooe(src_e,r_o);  std::cout<<GridLogMessage<<"Applied Meo"<<std::endl;
  Dwc.Meooe(src_o,r_e);  std::cout<<GridLogMessage<<"Applied Moe"<<std::endl;
  Dwc.Dhop (src,ref,DaggerNo);

  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  err= ref - r_eo;
  std::cout<<GridLogMessage << "EO norm diff   "<< norm2(err)<< " "<<norm2(ref)<< " " << norm2(r_eo) <<std::endl;

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

  Dwc.Meooe(chi_e,dchi_o);
  Dwc.Meooe(chi_o,dchi_e);
  Dwc.MeooeDag(phi_e,dphi_o);
  Dwc.MeooeDag(phi_o,dphi_e);

  ComplexD pDce = innerProduct(phi_e,dchi_e);
  ComplexD pDco = innerProduct(phi_o,dchi_o);
  ComplexD cDpe = innerProduct(chi_e,dphi_e);
  ComplexD cDpo = innerProduct(chi_o,dphi_o);

  std::cout<<GridLogMessage <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout<<GridLogMessage <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout<<GridLogMessage <<"pDce - conj(cDpo) "<< pDce-conj(cDpo) <<std::endl;
  std::cout<<GridLogMessage <<"pDco - conj(cDpe) "<< pDco-conj(cDpe) <<std::endl;
  std::cout<<GridLogMessage <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout<<GridLogMessage <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout<<GridLogMessage <<"pDce - conj(cDpo) "<< pDce-conj(cDpo) <<std::endl;
  std::cout<<GridLogMessage <<"pDco - conj(cDpe) "<< pDco-conj(cDpe) <<std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MeeInv Mee = 1                                         "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Dwc.Mooee(chi_e,src_e);
  Dwc.MooeeInv(src_e,phi_e);

  Dwc.Mooee(chi_o,src_o);
  Dwc.MooeeInv(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MeeInvDag MeeDag = 1                                   "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Dwc.MooeeDag(chi_e,src_e);
  Dwc.MooeeInvDag(src_e,phi_e);

  Dwc.MooeeDag(chi_o,src_o);
  Dwc.MooeeInvDag(src_o,phi_o);
  
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

  SchurDiagMooeeOperator<WilsonCloverFermionR,FermionField> HermOpEO(Dwc);
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
