    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_gpforce.cc

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

  const int Ls=8;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  typedef typename GparityDomainWallFermionR::FermionField FermionField;

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          RNG5(FGrid);  RNG5.SeedRandomDevice();
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedRandomDevice();
  
  FermionField phi        (FGrid); gaussian(RNG5,phi);
  FermionField Mphi       (FGrid); 
  FermionField MphiPrime  (FGrid); 

  LatticeGaugeField U(UGrid);

  SU3::HotConfiguration(RNG4,U);
  //  SU3::ColdConfiguration(pRNG,U);
  
  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=0.2; //kills the diagonal term
  RealD M5=1.8;
  //  const int nu = 3;
  //  std::vector<int> twists(Nd,0); // twists[nu] = 1;
  //  GparityDomainWallFermionR::ImplParams params;  params.twists = twists;
  //  GparityDomainWallFermionR Ddwf(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,params);

  //  DomainWallFermionR Dw     (U,     Grid,RBGrid,mass,M5);

  const int nu = 3;
  std::vector<int> twists(Nd,0);
  twists[nu] = 1;
  GparityDomainWallFermionR::ImplParams params;
  params.twists = twists;

  GparityDomainWallFermionR Dw(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,params);

  Dw.M   (phi,Mphi);

  ComplexD S    = innerProduct(Mphi,Mphi); // pdag MdagM p

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(UGrid);
  LatticeGaugeField tmp(UGrid);

  Dw.MDeriv(tmp , Mphi,  phi,DaggerNo );  UdSdU=tmp;
  Dw.MDeriv(tmp , phi,  Mphi,DaggerYes ); UdSdU=(UdSdU+tmp);
  
  FermionField Ftmp      (FGrid);

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 0.0001;
  RealD Hmom = 0.0;
  RealD Hmomprime = 0.0;
  RealD Hmompp    = 0.0;
  LatticeColourMatrix mommu(UGrid); 
  LatticeColourMatrix forcemu(UGrid); 
  LatticeGaugeField mom(UGrid); 
  LatticeGaugeField Uprime(UGrid); 

  for(int mu=0;mu<Nd;mu++){

    SU3::GaussianLieAlgebraMatrix(RNG4, mommu); // Traceless antihermitian momentum; gaussian in lie alg

    Hmom -= real(sum(trace(mommu*mommu)));

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

  std::cout << GridLogMessage <<"Initial mom hamiltonian is "<< Hmom <<std::endl;
  Dw.ImportGauge(Uprime);
  Dw.M          (phi,MphiPrime);

  ComplexD Sprime    = innerProduct(MphiPrime   ,MphiPrime);

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////


  for(int mu=0;mu<Nd;mu++){
    std::cout << "" <<std::endl;
    mommu   = PeekIndex<LorentzIndex>(mom,mu);
    std::cout << GridLogMessage<< " Mommu  " << norm2(mommu)<<std::endl;
    mommu   = mommu+adj(mommu);
    std::cout << GridLogMessage<< " Mommu + Mommudag " << norm2(mommu)<<std::endl;
    mommu   = PeekIndex<LorentzIndex>(UdSdU,mu);
    std::cout << GridLogMessage<< " dsdumu  " << norm2(mommu)<<std::endl;
    mommu   = mommu+adj(mommu);
    std::cout << GridLogMessage<< " dsdumu + dag  " << norm2(mommu)<<std::endl;
  }

  LatticeComplex dS(UGrid); dS = zero;
  LatticeComplex dSmom(UGrid); dSmom = zero;
  LatticeComplex dSmom2(UGrid); dSmom2 = zero;
  for(int mu=0;mu<Nd;mu++){
    mommu   = PeekIndex<LorentzIndex>(UdSdU,mu);
    mommu=Ta(mommu)*2.0;
    PokeIndex<LorentzIndex>(UdSdU,mommu,mu);
  }

  for(int mu=0;mu<Nd;mu++){
    mommu   = PeekIndex<LorentzIndex>(mom,mu);
    std::cout << GridLogMessage<< " Mommu  " << norm2(mommu)<<std::endl;
    mommu   = mommu+adj(mommu);
    std::cout << GridLogMessage<< " Mommu + Mommudag " << norm2(mommu)<<std::endl;
    mommu   = PeekIndex<LorentzIndex>(UdSdU,mu);
    std::cout << GridLogMessage<< " dsdumu  " << norm2(mommu)<<std::endl;
    mommu   = mommu+adj(mommu);
    std::cout << GridLogMessage<< " dsdumu + dag  " << norm2(mommu)<<std::endl;
  }

  for(int mu=0;mu<Nd;mu++){
    forcemu = PeekIndex<LorentzIndex>(UdSdU,mu);
    mommu   = PeekIndex<LorentzIndex>(mom,mu);

    // Update PF action density
    dS = dS+trace(mommu*forcemu)*dt;

    dSmom  = dSmom  - trace(mommu*forcemu) * dt;
    dSmom2 = dSmom2 - trace(forcemu*forcemu) *(0.25* dt*dt);

    // Update mom action density
    mommu = mommu + forcemu*(dt*0.5);

    Hmomprime -= real(sum(trace(mommu*mommu)));

  }

  Complex dSpred    = sum(dS);
  Complex dSm       = sum(dSmom);
  Complex dSm2      = sum(dSmom2);


  std::cout << GridLogMessage <<"Initial mom hamiltonian is "<< Hmom <<std::endl;
  std::cout << GridLogMessage <<"Final   mom hamiltonian is "<< Hmomprime <<std::endl;
  std::cout << GridLogMessage <<"Delta   mom hamiltonian is "<< Hmomprime-Hmom <<std::endl;

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "predict dS    "<< dSpred <<std::endl;
  std::cout << GridLogMessage <<"dSm "<< dSm<<std::endl;
  std::cout << GridLogMessage <<"dSm2"<< dSm2<<std::endl;

  std::cout << GridLogMessage << "Total dS    "<< Hmomprime - Hmom + Sprime - S <<std::endl;


  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
