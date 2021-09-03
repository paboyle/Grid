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

  SU<Nc>::HotConfiguration(RNG4,U);
  
  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=0.01; 
  RealD M5=1.8; 
  RealD b=0.5;
  RealD c=0.5;

  WilsonImplParams p;
  p.boundary_phases[0] = 1.0;
  p.boundary_phases[1] = 1.0;
  p.boundary_phases[2] = 1.0;
  p.boundary_phases[3] =- 1.0;
  
  MobiusFermionR Ddwf(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,p);
  Ddwf.M   (phi,Mphi);

  ComplexD S    = innerProduct(Mphi,Mphi); // pdag MdagM p

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(UGrid);
  LatticeGaugeField tmp(UGrid);
  

  Ddwf.MDeriv(tmp , Mphi,  phi,DaggerNo );  UdSdU=tmp;
  Ddwf.MDeriv(tmp , phi,  Mphi,DaggerYes ); UdSdU=(UdSdU+tmp);  
  
  // *****************************************************************************************
  // *** There is a funny negative sign in all derivatives. This is - UdSdU.               ***
  // ***                                                                                   ***
  // *** Deriv in both Wilson gauge action and the TwoFlavour.h seems to miss a minus sign ***
  // *** UdSdU is negated relative to what I think - call what is returned mUdSdU,         ***
  // *** and insert minus sign                                                             ***
  // *****************************************************************************************

  UdSdU = - UdSdU ; // Follow sign convention of actions in Grid. Seems crazy.

  LatticeFermion Ftmp      (FGrid);

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 0.001;
  RealD Hmom = 0.0;
  RealD Hmomprime = 0.0;

  LatticeColourMatrix mommu(UGrid); 
  LatticeColourMatrix mUdSdUmu(UGrid); 
  LatticeGaugeField mom(UGrid); 
  LatticeGaugeField Uprime(UGrid); 

  for(int mu=0;mu<Nd;mu++){

    SU<Nc>::GaussianFundamentalLieAlgebraMatrix(RNG4, mommu); // Traceless antihermitian momentum; gaussian in lie alg
    PokeIndex<LorentzIndex>(mom,mommu,mu);

    // Momentum Hamiltonian is - trace(p^2)/HMC_MOM_DENOMINATOR
    //
    // Integrator.h:   RealD H = - FieldImplementation::FieldSquareNorm(P)/HMC_MOMENTUM_DENOMINATOR; // - trace (P*P)/denom                                                                                       //     GaugeImplTypes.h:        Hloc += trace(Pmu * Pmu);
    //                          Sign comes from a sneaky multiply by "i" in GaussianFundemantalLie algebra
    //                          P is i P^a_\mu T^a, not Pa Ta
    // 
    // Integrator.h: H =  Hmom + sum S(action)
    Hmom -= real(sum(trace(mommu*mommu)))/ HMC_MOMENTUM_DENOMINATOR;
    
    // fourth order exponential approx
    autoView( U_v , U, CpuRead);
    autoView( mom_v, mom, CpuRead);
    autoView(Uprime_v, Uprime, CpuWrite);
    thread_foreach(i,mom_v,{
      Uprime_v[i](mu) =  U_v[i](mu)
	+ mom_v[i](mu)*U_v[i](mu)*dt 
	+ mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt/2.0)
	+ mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt*dt/6.0)
	+ mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt*dt*dt/24.0)
	+ mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt*dt*dt*dt/120.0)
	+ mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *mom_v[i](mu) *U_v[i](mu)*(dt*dt*dt*dt*dt*dt/720.0)
	;
    });
  }
  std::cout << GridLogMessage <<"Initial mom hamiltonian is "<< Hmom <<std::endl;
  
  Ddwf.ImportGauge(Uprime);
  Ddwf.M          (phi,MphiPrime);

  ComplexD Sprime    = innerProduct(MphiPrime   ,MphiPrime);

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////

  LatticeComplex dS(UGrid); dS = Zero();
  LatticeComplex dSmom(UGrid); dSmom = Zero();
  LatticeComplex dSmom2(UGrid); dSmom2 = Zero();

  for(int mu=0;mu<Nd;mu++){
    mommu   = PeekIndex<LorentzIndex>(UdSdU,mu);
    mommu=Ta(mommu);
    PokeIndex<LorentzIndex>(UdSdU,mommu,mu);
  }

  for(int mu=0;mu<Nd;mu++){

    mUdSdUmu= PeekIndex<LorentzIndex>(UdSdU,mu);
    mommu   = PeekIndex<LorentzIndex>(mom,mu);

  //
  // Derive HMC eom:
  //
  // Sdot =  - 2 trace( p p^dot ) / D - trace( p [ mUdSdU - h.c. ] ) = 0 
  //
  //
  // Sdot = 0 = - 2 trace( p p^dot ) / D - 2 trace( p Ta( mUdSdU ) = 0
  //
  // EOM: 
  //
  // pdot = - D Ta( mUdSdU ) -- source of sign is the "funny sign" above
  //
  // dSqcd_dt  = - 2.0*trace(mommu* Ta(mUdSdU) )*dt -- i.e. mUdSdU with adjoint term -> force has a 2x implicit
  //
  // dH_mom/dt = - 2 trace (p pdot)/Denom  
  //
  // dH_tot / dt = 0 <= pdot = - Denom * mUdSdU 
  //
  // dH_mom/dt = 2 trace (p mUdSdU ) 
  //
  // True Momentum delta H has a dt^2 piece
  //
  // dSmom = [ trace mom*mom - trace ( (mom-Denom*f*dt)(mom-Denom*f*dt) ) ] / Denom
  //       = 2*trace(mom*f) dt  - Denom*dt*dt * trace(f*f).
  //       = dSmom + dSmom2
  //
    
    dS     = dS - 2.0*trace(mommu*mUdSdUmu)*dt;   // U and Udagger derivs hence 2x.

    dSmom  = dSmom  + 2.0*trace(mommu*mUdSdUmu) * dt;  // this 2.0 coms from derivative of p^2 
    
    dSmom2 = dSmom2 - trace(mUdSdUmu*mUdSdUmu) * dt*dt* HMC_MOMENTUM_DENOMINATOR; // Remnant

    mommu = mommu - mUdSdUmu * dt* HMC_MOMENTUM_DENOMINATOR;; 

    Hmomprime -= real(sum(trace(mommu*mommu))) / HMC_MOMENTUM_DENOMINATOR;

  }

  ComplexD dSpred    = sum(dS);
  ComplexD dSm       = sum(dSmom);
  ComplexD dSm2      = sum(dSmom2);

  std::cout << GridLogMessage <<"dSm "<< dSm<<std::endl;
  std::cout << GridLogMessage <<"dSm2 "<< dSm2<<std::endl;

  std::cout << GridLogMessage <<"Initial mom hamiltonian is "<< Hmom <<std::endl;
  std::cout << GridLogMessage <<"Final   mom hamiltonian is "<< Hmomprime <<std::endl;

  std::cout << GridLogMessage <<"Delta   mom hamiltonian is "<< Hmomprime-Hmom <<std::endl;
  std::cout << GridLogMessage <<"predict Delta mom hamiltonian is "<< dSm+dSm2 <<std::endl;
  
  std::cout << GridLogMessage << "Initial S      "<<S<<std::endl;
  std::cout << GridLogMessage << "Final   S      "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "Delta   S      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "predict delta S"<< dSpred <<std::endl;
  std::cout << GridLogMessage << "defect "<< Sprime-S-dSpred <<std::endl;

  std::cout << GridLogMessage << "Total dS    "<< Hmomprime - Hmom + Sprime - S <<std::endl;

  std::cout << GridLogMessage << "dS - dt^2 term "<< Hmomprime - Hmom + Sprime - S - dSm2 <<std::endl;

  assert( fabs(real(Sprime-S-dSpred)) < 1.0 ) ;


  
  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
