    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/forces/Test_gpdwf_Xconj_force.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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

typedef GparityDomainWallFermionD::FermionField FermionField2f;
typedef XconjugateDomainWallFermionD::FermionField FermionField1f;

const Gamma & Xmatrix(){
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  return X;
}

void XconjBoost(FermionField2f &to, const FermionField1f &from){
  FermionField1f tmp = -(Xmatrix()*conjugate(from));
  PokeIndex<GparityFlavourIndex>(to, from, 0);
  PokeIndex<GparityFlavourIndex>(to, tmp, 1);
  to.Checkerboard() = from.Checkerboard();
}

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
  
  std::vector<int> seeds4({1,2,3,5});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  LatticeGaugeField U(UGrid);

  SU<Nc>::HotConfiguration(RNG4,U);
  
  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=0.01; 
  RealD M5=1.8; 
  
  XconjugateDomainWallFermionD::ImplParams xparams;
  std::vector<int> twists({1,1,1,1}); //GPBC in 3 dirs, antiperiodic in time
  xparams.twists = twists;
  xparams.boundary_phase = 1.0;

  XconjugateDomainWallFermionD Ddwf(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,xparams);

  GparityDomainWallFermionD::ImplParams params;  params.twists = twists;
  GparityDomainWallFermionD Ddwf_gp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,params);
  typedef GparityDomainWallFermionD::FermionField FermionField2f;

  //Action of X-conjugate operator
  LatticeFermion phi        (FGrid); gaussian(RNG5,phi);
  LatticeFermion Mphi       (FGrid); 
  LatticeFermion MphiPrime  (FGrid); 
  Ddwf.M   (phi,Mphi);
  ComplexD S    = innerProduct(Mphi,Mphi); // pdag MdagM p

  //Check against action of G-parity operator
  FermionField2f phi_gp        (FGrid);
  FermionField2f Mphi_gp       (FGrid); 
  FermionField2f MphiPrime_gp  (FGrid);
  XconjBoost(phi_gp, phi);
  phi_gp = phi_gp * sqrt(0.5); //rho from rho'
 
  Ddwf_gp.M   (phi_gp,Mphi_gp);
  ComplexD S_gp    = innerProduct(Mphi_gp,Mphi_gp); // pdag MdagM p

  std::cout << GridLogMessage << "Initial action Xconjugate: " << S << " Gparity: " << S_gp << " diff: " << S_gp - S << std::endl;


  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(UGrid);
  LatticeGaugeField tmp(UGrid);
  Ddwf.MDeriv(tmp , Mphi,  phi,DaggerNo );  UdSdU=tmp;
  Ddwf.MDeriv(tmp , phi,  Mphi,DaggerYes ); UdSdU=(UdSdU+tmp);  
  
  //and the same for the G-parity operator
  LatticeGaugeField UdSdU_gp(UGrid);
  LatticeGaugeField tmp_gp(UGrid);
  Ddwf_gp.MDeriv(tmp_gp , Mphi_gp,  phi_gp,DaggerNo );  UdSdU_gp=tmp_gp;
  Ddwf_gp.MDeriv(tmp_gp , phi_gp,  Mphi_gp,DaggerYes ); UdSdU_gp=(UdSdU_gp+tmp_gp);  
  
  LatticeGaugeField force_diff = UdSdU_gp -  UdSdU;
  std::cout << GridLogMessage << "Force Xconjugate: " << norm2(UdSdU) << " G-parity: " << norm2(UdSdU_gp) << " diff: " << norm2(force_diff) << std::endl;
  


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

    SU<Nc>::GaussianFundamentalLieAlgebraMatrix(RNG4, mommu); // Traceless antihermitian momentum; gaussian in lie alg

    PokeIndex<LorentzIndex>(mom,mommu,mu);

    // fourth order exponential approx

    autoView( mom_v, mom, CpuRead);
    autoView( U_v , U, CpuRead);
    autoView(Uprime_v, Uprime, CpuWrite);

    thread_foreach( i,mom_v,{
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

  Ddwf_gp.ImportGauge(Uprime);
  Ddwf_gp.M          (phi_gp,MphiPrime_gp);
  ComplexD Sprime_gp    = innerProduct(MphiPrime_gp   ,MphiPrime_gp);

  std::cout << GridLogMessage << "Final action Xconjugate: " << Sprime << " Gparity: " << Sprime_gp << " diff: " << Sprime_gp - Sprime << std::endl;


  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////

  LatticeComplex dS(UGrid); dS = Zero();
  LatticeComplex dS_gp(UGrid); dS_gp = Zero();

  for(int mu=0;mu<Nd;mu++){
    mommu   = PeekIndex<LorentzIndex>(UdSdU,mu);
    mommu=Ta(mommu)*2.0;
    PokeIndex<LorentzIndex>(UdSdU,mommu,mu);

    mommu   = PeekIndex<LorentzIndex>(UdSdU_gp,mu);
    mommu=Ta(mommu)*2.0;
    PokeIndex<LorentzIndex>(UdSdU_gp,mommu,mu);
  }

  for(int mu=0;mu<Nd;mu++){
    mommu   = PeekIndex<LorentzIndex>(mom,mu);

    // Update PF action density
    forcemu = PeekIndex<LorentzIndex>(UdSdU,mu);
    dS = dS+trace(mommu*forcemu)*dt;

    forcemu = PeekIndex<LorentzIndex>(UdSdU_gp,mu);
    dS_gp = dS_gp+trace(mommu*forcemu)*dt;
  }

  ComplexD dSpred    = sum(dS);
  ComplexD dSpred_gp    = sum(dS_gp);
  std::cout << std::setprecision(14);

  std::cout << GridLogMessage << "              Xconj               Gparity                 Diff" << std::endl;      
  std::cout << GridLogMessage << "S            "<< S<< " " << S_gp << " " << S_gp-S << std::endl;
  std::cout << GridLogMessage << "Sprime       "<< Sprime << " " << Sprime_gp << " " << Sprime_gp - Sprime <<  std::endl;
  std::cout << GridLogMessage << "dS            "<< Sprime-S << " " << Sprime_gp - S_gp << " " << (Sprime-S)-(Sprime_gp-S_gp) << std::endl;
  std::cout << GridLogMessage << "predict dS    "<< dSpred << " " << dSpred_gp << " " << dSpred_gp - dSpred << std::endl;

  assert( fabs(real(Sprime-S-dSpred)) < 1.0 ) ;

  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
