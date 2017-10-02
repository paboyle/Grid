    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_wilson.cc

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
  GridRedBlackCartesian     RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALF"<< sizeof(RealF)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALD"<< sizeof(RealD)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REAL"<< sizeof(Real)<<std::endl;

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);
  //  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9});

  typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 
  typedef typename ImprovedStaggeredFermionR::ComplexField ComplexField; 
  typename ImprovedStaggeredFermionR::ImplParams params; 

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
  RealD c1=9.0/8.0;
  RealD c2=-1.0/24.0;
  RealD u0=1.0;

  { // Simple improved staggered implementation
    ref = zero;
    RealD c1tad = 0.5*c1/u0;
    RealD c2tad = 0.5*c2/u0/u0/u0;

    Lattice<iScalar<vInteger> > coor(&Grid);

    Lattice<iScalar<vInteger> > x(&Grid); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(&Grid); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(&Grid); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(&Grid); LatticeCoordinate(t,3);

    Lattice<iScalar<vInteger> > lin_z(&Grid); lin_z=x+y;
    Lattice<iScalar<vInteger> > lin_t(&Grid); lin_t=x+y+z;

    for(int mu=0;mu<Nd;mu++){

      // Staggered Phase.
      ComplexField phases(&Grid);	phases=1.0;
  
      if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
      if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
      if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);

      tmp = PeriodicBC::CovShiftForward(U[mu],mu,src);
      ref = ref +c1tad*tmp*phases; // Forward 1 hop

      tmp = PeriodicBC::CovShiftForward(U[mu],mu,tmp); // 2 hop
      tmp = PeriodicBC::CovShiftForward(U[mu],mu,tmp); // 3 hop
      ref = ref +c2tad*tmp*phases; // Forward 3 hop

      tmp = PeriodicBC::CovShiftBackward(U[mu],mu,src);
      ref = ref -c1tad*tmp*phases; // Forward 3 hop

      tmp = PeriodicBC::CovShiftBackward(U[mu],mu,tmp);
      tmp = PeriodicBC::CovShiftBackward(U[mu],mu,tmp);
      ref = ref -c2tad*tmp*phases; // Forward 3 hop
    }
    //    ref = ref + mass * src;
  }

  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,mass,c1,c2,u0,params);
  

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing Dhop against cshift implementation         "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  std::cout<<GridLogMessage << "Calling Ds"<<std::endl;
  int ncall=1000;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Ds.Dhop(src,result,0);
  }
  double t1=usecond();
  double t2;
  double flops=(16*(3*(6+8+8)) + 15*3*2)*volume*ncall; // == 66*16 +  == 1146
  
  std::cout<<GridLogMessage << "Called Ds"<<std::endl;
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

  Ds.Meooe(src_e,r_o);  std::cout<<GridLogMessage<<"Applied Meo"<<std::endl;
  Ds.Meooe(src_o,r_e);  std::cout<<GridLogMessage<<"Applied Moe"<<std::endl;
  Ds.Dhop (src,ref,DaggerNo);

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

  Ds.Meooe(chi_e,dchi_o);
  Ds.Meooe(chi_o,dchi_e);
  Ds.MeooeDag(phi_e,dphi_o);
  Ds.MeooeDag(phi_o,dphi_e);

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

  Ds.Mooee(chi_e,src_e);
  Ds.MooeeInv(src_e,phi_e);

  Ds.Mooee(chi_o,src_o);
  Ds.MooeeInv(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Test MeeInvDag MeeDag = 1                                   "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Ds.MooeeDag(chi_e,src_e);
  Ds.MooeeInvDag(src_e,phi_e);

  Ds.MooeeDag(chi_o,src_o);
  Ds.MooeeInvDag(src_o,phi_o);
  
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

  SchurDiagMooeeOperator<ImprovedStaggeredFermionR,FermionField> HermOpEO(Ds);
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
