    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_force.cc

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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  LatticeFermion phi        (FGrid); gaussian(RNG5,phi);
  LatticeFermion Mphi       (FGrid); 
  LatticeFermion MphiPrime  (FGrid); 

  LatticeGaugeField U(UGrid);

  SU3::HotConfiguration(RNG4,U);
  
  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=0.01; 
  RealD M5=1.8; 
  RealD b=0.5;
  RealD c=0.5;

  std::vector < std::complex<double>  > omegas;
  omegas.push_back( std::complex<double>(1.45806438985048,-0) );
  omegas.push_back( std::complex<double>(1.18231318389348,-0) );
  omegas.push_back( std::complex<double>(0.830951166685955,-0) );
  omegas.push_back( std::complex<double>(0.542352409156791,-0) );
  omegas.push_back( std::complex<double>(0.341985020453729,-0) );
  omegas.push_back( std::complex<double>(0.21137902619029,-0) );
  omegas.push_back( std::complex<double>(0.126074299502912,-0) );
  omegas.push_back( std::complex<double>(0.0990136651962626,-0) );
  omegas.push_back( std::complex<double>(0.0686324988446592,0.0550658530827402) );
  omegas.push_back( std::complex<double>(0.0686324988446592,-0.0550658530827402) );

  ZMobiusFermionR Ddwf(U, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, omegas,b,c);

  Ddwf.M   (phi,Mphi);

  ComplexD S    = innerProduct(Mphi,Mphi); // pdag MdagM p

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(UGrid);
  LatticeGaugeField tmp(UGrid);
  

  Ddwf.MDeriv(tmp , Mphi,  phi,DaggerNo );  UdSdU=tmp;
  Ddwf.MDeriv(tmp , phi,  Mphi,DaggerYes ); UdSdU=(UdSdU+tmp);  
  
  LatticeFermion Ftmp      (FGrid);

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 0.0001;

  LatticeColourMatrix mommu(UGrid); 
  LatticeColourMatrix forcemu(UGrid); 
  LatticeGaugeField mom(UGrid); 
  LatticeGaugeField Uprime(UGrid); 

  for(int mu=0;mu<Nd;mu++){

    SU3::GaussianFundamentalLieAlgebraMatrix(RNG4, mommu); // Traceless antihermitian momentum; gaussian in lie alg

    PokeIndex<LorentzIndex>(mom,mommu,mu);

    // fourth order exponential approx
    auto mom_v    = mom.View();
    auto U_v      = U.View();
    auto Uprime_v = Uprime.View();

    thread_loop( (auto i=mom_v.begin();i<mom_v.end();i++),{
      Uprime_v[i](mu) =	  U_v[i](mu)
	+ mom_v[i](mu)*U_v[i](mu)*dt 
	+ mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt/2.0)
	+ mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt*dt/6.0)
	+ mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt*dt*dt/24.0)
	+ mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt*dt*dt*dt/120.0)
	+ mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt*dt*dt*dt*dt/720.0)
	;
    });
  }
  
  Ddwf.ImportGauge(Uprime);
  Ddwf.M          (phi,MphiPrime);

  ComplexD Sprime    = innerProduct(MphiPrime   ,MphiPrime);

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////


  LatticeComplex dS(UGrid); dS = Zero();
  for(int mu=0;mu<Nd;mu++){
    mommu   = PeekIndex<LorentzIndex>(UdSdU,mu);
    mommu=Ta(mommu)*2.0;
    PokeIndex<LorentzIndex>(UdSdU,mommu,mu);
  }

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

  assert( fabs(real(Sprime-S-dSpred)) < 3.0 ) ;

  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
