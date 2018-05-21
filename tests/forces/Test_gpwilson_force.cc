    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_gpdwf_force.cc

    Copyright (C) 2015

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

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = UGrid;
  GridRedBlackCartesian * FrbGrid = UrbGrid;

  std::vector<int> seeds4({1,2,3,4});
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  typedef typename GparityDomainWallFermionR::FermionField FermionField;
  FermionField phi        (FGrid); gaussian(RNG4,phi);
  FermionField Mphi       (FGrid); 
  FermionField MphiPrime  (FGrid); 

  LatticeGaugeField U(UGrid);

  SU3::HotConfiguration(RNG4,U);
  
  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=0.01; 
  RealD M5=1.8; 

  const int nu = 3;
  std::vector<int> twists(Nd,0);  twists[nu] = 1;
  GparityWilsonFermionR::ImplParams params;  params.twists = twists;
  GparityWilsonFermionR Wil(U,*UGrid,*UrbGrid,mass,params);
  Wil.M   (phi,Mphi);

  ComplexD S    = innerProduct(Mphi,Mphi); // pdag MdagM p

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(UGrid);
  LatticeGaugeField tmp(UGrid);

  Wil.MDeriv(tmp , Mphi,  phi,DaggerNo );  UdSdU=tmp;
  Wil.MDeriv(tmp , phi,  Mphi,DaggerYes ); UdSdU=(UdSdU+tmp);  
  
  FermionField Ftmp      (FGrid);

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 0.01;

  LatticeColourMatrix mommu(UGrid); 
  LatticeColourMatrix forcemu(UGrid); 
  LatticeGaugeField mom(UGrid); 
  LatticeGaugeField Uprime(UGrid); 

  for(int mu=0;mu<Nd;mu++){

    // Traceless antihermitian momentum; gaussian in lie alg
    SU3::GaussianFundamentalLieAlgebraMatrix(RNG4, mommu); 

    PokeIndex<LorentzIndex>(mom,mommu,mu);

    // fourth order exponential approx
    parallel_for(auto i=mom.begin();i<mom.end();i++){
      Uprime[i](mu) =
	  U[i](mu)
	+ mom[i](mu)*U[i](mu)*dt 
	+ mom[i](mu) *mom[i](mu) *U[i](mu)*(dt*dt/2.0)
	+ mom[i](mu) *mom[i](mu) *mom[i](mu) *U[i](mu)*(dt*dt*dt/6.0)
	+ mom[i](mu) *mom[i](mu) *mom[i](mu) *mom[i](mu) *U[i](mu)*(dt*dt*dt*dt/24.0)
	+ mom[i](mu) *mom[i](mu) *mom[i](mu) *mom[i](mu) *mom[i](mu) *U[i](mu)*(dt*dt*dt*dt*dt/120.0)
	+ mom[i](mu) *mom[i](mu) *mom[i](mu) *mom[i](mu) *mom[i](mu) *mom[i](mu) *U[i](mu)*(dt*dt*dt*dt*dt*dt/720.0)
	;
    }
  }
  
  Wil.ImportGauge(Uprime);
  Wil.M          (phi,MphiPrime);

  ComplexD Sprime    = innerProduct(MphiPrime   ,MphiPrime);

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////

  for(int mu=0;mu<Nd;mu++){
    mommu   = PeekIndex<LorentzIndex>(UdSdU,mu);
    mommu=Ta(mommu)*2.0;
    PokeIndex<LorentzIndex>(UdSdU,mommu,mu);
  }

  LatticeComplex dS(UGrid); dS = zero;
  for(int mu=0;mu<Nd;mu++){
    forcemu = PeekIndex<LorentzIndex>(UdSdU,mu);
    mommu   = PeekIndex<LorentzIndex>(mom,mu);

    // Update PF action density
    dS = dS+trace(mommu*forcemu)*dt;
  }

  ComplexD dSpred    = sum(dS);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "predict dS    "<< dSpred <<std::endl;

  assert( fabs(real(Sprime-S-dSpred)) < 2.0 ) ;

  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
