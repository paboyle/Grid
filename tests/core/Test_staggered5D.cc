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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  std::cout << GridLogMessage << "Making s innermost grids"<<std::endl;

  const int Ls=16;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  int threads = GridThread::GetThreads();

  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALF"<< sizeof(RealF)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALD"<< sizeof(RealD)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REAL"<< sizeof(Real)<<std::endl;

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG4(UGrid);
  GridParallelRNG          pRNG5(FGrid);
  pRNG4.SeedFixedIntegers(seeds);
  pRNG5.SeedFixedIntegers(seeds);

  typedef typename ImprovedStaggeredFermion5DR::FermionField FermionField; 
  typedef typename ImprovedStaggeredFermion5DR::ComplexField ComplexField; 
  typename ImprovedStaggeredFermion5DR::ImplParams params; 

  FermionField src   (FGrid);

  random(pRNG5,src);

  FermionField result(FGrid); result=Zero();
  FermionField    ref(FGrid);    ref=Zero();
  FermionField    tmp(FGrid);    tmp=Zero();
  FermionField    err(FGrid);    tmp=Zero();
  FermionField phi   (FGrid); random(pRNG5,phi);
  FermionField chi   (FGrid); random(pRNG5,chi);

  LatticeGaugeField Umu(UGrid); SU3::ColdConfiguration(pRNG4,Umu);
  LatticeGaugeField Umua(UGrid); Umua=Umu;

  double volume=Ls;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  ////////////////////////////////////
  // Naive implementation needs to
  // replicate across fifth dimension
  ////////////////////////////////////
  LatticeGaugeField Umu5d(FGrid); 
  auto umu5d = Umu5d.View();
  auto umu   = Umu.View();
  for(int ss=0;ss<Umu.Grid()->oSites();ss++){
    for(int s=0;s<Ls;s++){
      umu5d[Ls*ss+s] = umu[ss];
    }
  }

  std::vector<LatticeColourMatrix> U(4,FGrid);

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu5d,mu);
  }

  RealD mass=0.1;
  RealD c1=9.0/8.0;
  RealD c2=-1.0/24.0;
  RealD u0=1.0;

  { // Simple improved staggered implementation
    ref = Zero();
    RealD c1tad = 0.5*c1/u0;
    RealD c2tad = 0.5*c2/u0/u0/u0;

    Lattice<iScalar<vInteger> > coor(FGrid);

    Lattice<iScalar<vInteger> > x(FGrid); LatticeCoordinate(x,1); // s innermost
    Lattice<iScalar<vInteger> > y(FGrid); LatticeCoordinate(y,2);
    Lattice<iScalar<vInteger> > z(FGrid); LatticeCoordinate(z,3);
    Lattice<iScalar<vInteger> > t(FGrid); LatticeCoordinate(t,4);

    Lattice<iScalar<vInteger> > lin_z(FGrid); lin_z=x+y;
    Lattice<iScalar<vInteger> > lin_t(FGrid); lin_t=x+y+z;

    for(int mu=0;mu<Nd;mu++){

      // Staggered Phase.
      ComplexField phases(FGrid);	phases=1.0;
  
      if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
      if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
      if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);

      tmp = PeriodicBC::CovShiftForward(U[mu],mu+1,src);
      ref = ref +c1tad*tmp*phases; // Forward 1 hop

      tmp = PeriodicBC::CovShiftForward(U[mu],mu+1,tmp); // 2 hop
      tmp = PeriodicBC::CovShiftForward(U[mu],mu+1,tmp); // 3 hop
      ref = ref +c2tad*tmp*phases; // Forward 3 hop

      tmp = PeriodicBC::CovShiftBackward(U[mu],mu+1,src);
      ref = ref -c1tad*tmp*phases; // Forward 3 hop

      tmp = PeriodicBC::CovShiftBackward(U[mu],mu+1,tmp);
      tmp = PeriodicBC::CovShiftBackward(U[mu],mu+1,tmp);
      ref = ref -c2tad*tmp*phases; // Forward 3 hop
    }
  }

  ImprovedStaggeredFermion5DR Ds(Umu,Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,c1,c2,u0,params);

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
  //  std::cout<<GridLogMessage << "result"<< result <<std::endl;
  //  std::cout<<GridLogMessage << "ref   "<< ref   <<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;

  err = ref-result; 
  std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;

  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that Deo + Doe = Dunprec "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  FermionField src_e   (FrbGrid);
  FermionField src_o   (FrbGrid);
  FermionField r_e   (FrbGrid);
  FermionField r_o   (FrbGrid);
  FermionField r_eo  (FGrid);
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

  FermionField chi_e   (FrbGrid);
  FermionField chi_o   (FrbGrid);

  FermionField dchi_e  (FrbGrid);
  FermionField dchi_o  (FrbGrid);

  FermionField phi_e   (FrbGrid);
  FermionField phi_o   (FrbGrid);

  FermionField dphi_e  (FrbGrid);
  FermionField dphi_o  (FrbGrid);

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
  
  random(pRNG5,phi);
  random(pRNG5,chi);
  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);

  SchurDiagMooeeOperator<ImprovedStaggeredFermion5DR,FermionField> HermOpEO(Ds);
  HermOpEO.MpcDagMpc(chi_e,dchi_e);
  HermOpEO.MpcDagMpc(chi_o,dchi_o);

  HermOpEO.MpcDagMpc(phi_e,dphi_e);
  HermOpEO.MpcDagMpc(phi_o,dphi_o);

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
