    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_force_phiMphi.cc

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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

#define parallel_for PARALLEL_FOR_LOOP for

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
  pRNG.SeedRandomDevice();

  LatticeFermion phi        (&Grid); gaussian(pRNG,phi);
  LatticeFermion Mphi       (&Grid); 
  LatticeFermion MphiPrime  (&Grid); 
  LatticeFermion dMphi      (&Grid); 

  LatticeGaugeField U(&Grid);


  SU3::HotConfiguration(pRNG,U);

  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=-4.0; //kills the diagonal term
  WilsonFermionR Dw     (U,     Grid,RBGrid,mass);
  Dw.M(phi,Mphi);

  ComplexD S = innerProduct(phi,Mphi);

  // get the deriv
  LatticeGaugeField UdSdU(&Grid);
  Dw.MDeriv(UdSdU,phi, phi,DaggerNo ); 


  ////////////////////////////////////
  // Modify the gauge field a little in one dir
  ////////////////////////////////////
  RealD dt = 1.0e-3;
  Complex Complex_i(0,1);

  LatticeColourMatrix Umu(&Grid);
  LatticeColourMatrix Umu_save(&Grid);
  LatticeColourMatrix dU (&Grid);
  LatticeColourMatrix mom(&Grid); 
  SU3::GaussianLieAlgebraMatrix(pRNG, mom); // Traceless antihermitian momentum; gaussian in lie alg


  // check mom is as i expect
  LatticeColourMatrix tmpmom(&Grid); 
  tmpmom = mom+adj(mom);
  std::cout << GridLogMessage << "mom anti-herm check "<< norm2(tmpmom)<<std::endl;
  std::cout << GridLogMessage << "mom tr check "<< norm2(trace(mom))<<std::endl;
  
  const int mu=0;
  Umu = PeekIndex<LorentzIndex>(U,mu);
  Umu_save=Umu;
  dU = mom * Umu * dt;
  Umu= Umu+dU;
  PokeIndex<LorentzIndex>(Dw.Umu,Umu,mu);

  Dw.M(phi,MphiPrime);

  ComplexD Sprime = innerProduct(phi,MphiPrime);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;

  Dw.Umu=zero;
  PokeIndex<LorentzIndex>(Dw.Umu,dU,mu);
  Dw.M(phi,dMphi);


  ComplexD deltaS = innerProduct(phi,dMphi);
  std::cout << GridLogMessage << "deltaS      "<<deltaS<<std::endl;

  Dw.Umu=zero;
  PokeIndex<LorentzIndex>(Dw.Umu,Umu_save,mu);
  Dw.Mdir(phi,dMphi,mu,1);
  dMphi = dt*mom*dMphi;

  deltaS = innerProduct(phi,dMphi);
  std::cout << GridLogMessage << "deltaS from inner prod of mom* M[u]     "<<deltaS<<std::endl;

  deltaS = sum(trace(outerProduct(dMphi,phi)));
  std::cout << GridLogMessage << "deltaS from trace outer prod of deltaM      "<<deltaS<<std::endl;

/*
  LatticeComplex lip(&Grid);
  lip  = localInnerProduct(phi,dMphi);
  
  LatticeComplex trop(&Grid);
  trop = trace(outerProduct(dMphi,phi));

  LatticeSpinColourMatrix op(&Grid);
  op = outerProduct(dMphi,phi);

  LatticeSpinColourMatrix hop(&Grid);
  LatticeComplex op_cpt(&Grid);
  for(int s1=0;s1<Ns;s1++){
  for(int s2=0;s2<Ns;s2++){
  for(int c1=0;c1<Nc;c1++){
  for(int c2=0;c2<Nc;c2++){

    op_cpt = peekColour(peekSpin(dMphi,s1),c1) * adj(peekColour(peekSpin(phi,s2),c2));

    parallel_for(auto i=hop.begin();i<hop.end();i++){
      hop[i]()(s1,s2)(c1,c2) = op_cpt[i]()()();
    }
    
  }}}}

  LatticeSpinColourMatrix diffop(&Grid);

  diffop = hop - op;
  std::cout << GridLogMessage << "hand outer prod diff   "<<norm2(diffop)<<std::endl;

  deltaS = sum(trace(hop));
  std::cout << GridLogMessage << "deltaS hop   "<<deltaS<<std::endl;

  std::cout << GridLogMessage<< "  phi[0] : "<<  phi._odata[0]<<std::endl;
  std::cout << GridLogMessage<< "dMphi[0] : "<<dMphi._odata[0]<<std::endl;
  std::cout << GridLogMessage<< "hop[0]   : "<<  hop._odata[0]<<std::endl;
  std::cout << GridLogMessage<< " op[0]   : "<<   op._odata[0]<<std::endl;


  std::cout << GridLogMessage << "lip      "<<lip<<std::endl;
  std::cout << GridLogMessage << "trop     "<<trop<<std::endl;

*/  
  
  //  std::cout << GridLogMessage << " UdSdU " << UdSdU << std::endl;

  LatticeComplex dS(&Grid); dS = zero;
  parallel_for(auto i=mom.begin();i<mom.end();i++){
    dS[i]() = trace(mom[i]() * UdSdU[i](mu) )*dt;
  }
  Complex dSpred = sum(dS);

  std::cout << GridLogMessage << "predict dS    "<< dSpred <<std::endl;


  cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
