    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_rect_force.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>

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

#define parallel_for PARALLEL_FOR_LOOP for

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

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(std::vector<int>({15,91,21,3}));

  LatticeGaugeField U(&Grid);
  LatticeGaugeField P(&Grid);
  LatticeColourMatrix P_mu(&Grid);
  // Matrix in the algebra
  for (int mu = 0; mu < Nd; mu++) {
    SU<Nc>::GaussianFundamentalLieAlgebraMatrix(pRNG, P_mu);
    PokeIndex<LorentzIndex>(P, P_mu, mu);
  }

  SU3::HotConfiguration(pRNG,U);
  

  ConjugateGradient<LatticeGaugeField> CG(1.0e-8, 10000);
  LaplacianParams LapPar(0.001, 1.0, 1000, 1e-8, 10, 64);
  RealD Kappa = 0.99;
  LaplacianAdjointField<PeriodicGimplR> Laplacian(&Grid, CG, LapPar, Kappa);
  GeneralisedMomenta<PeriodicGimplR> LaplacianMomenta(&Grid, Laplacian);
  LaplacianMomenta.M.ImportGauge(U);
  LaplacianMomenta.MomentaDistribution(pRNG);// fills the Momenta with the correct distr
  

  std::cout << std::setprecision(15);
  std::cout << GridLogMessage << "MomentaAction" << std::endl;
  ComplexD S    = LaplacianMomenta.MomentaAction();

  // get the deriv with respect to "U"
  LatticeGaugeField UdSdU(&Grid);
  LatticeGaugeField AuxDer(&Grid);
  std::cout << GridLogMessage<< "DerivativeU" << std::endl;
  LaplacianMomenta.DerivativeU(LaplacianMomenta.Mom, UdSdU);
  LaplacianMomenta.AuxiliaryFieldsDerivative(AuxDer);
  UdSdU += AuxDer;

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 0.0001;

  LatticeColourMatrix mommu(&Grid); 
  LatticeColourMatrix forcemu(&Grid); 
  LatticeGaugeField mom(&Grid); 
  LatticeGaugeField Uprime(&Grid); 

  std::cout << GridLogMessage << "Update the U " << std::endl;
  for(int mu=0;mu<Nd;mu++){
  // Traceless antihermitian momentum; gaussian in lie algebra
    SU3::GaussianFundamentalLieAlgebraMatrix(pRNG, mommu); 
    auto Umu = PeekIndex<LorentzIndex>(U, mu);
    PokeIndex<LorentzIndex>(mom,mommu,mu);
    Umu = expMat(mommu, dt, 12) * Umu;
    PokeIndex<LorentzIndex>(Uprime, ProjectOnGroup(Umu), mu);
  
  }

  std::cout << GridLogMessage << "New action " << std::endl;
  LaplacianMomenta.M.ImportGauge(Uprime);
  ComplexD Sprime    = LaplacianMomenta.MomentaAction();

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////

  LatticeComplex dS(&Grid); dS = zero;

  for(int mu=0;mu<Nd;mu++){
    auto UdSdUmu = PeekIndex<LorentzIndex>(UdSdU,mu);
         mommu   = PeekIndex<LorentzIndex>(mom,mu);

    // Update gauge action density
    // U = exp(p dt) U
    // dU/dt = p U
    // so dSdt = trace( dUdt dSdU) = trace( p UdSdUmu ) 

    dS = dS + trace(mommu*UdSdUmu)*dt*2.0;
  }

  ComplexD dSpred    = sum(dS);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "pred dS "<< dSpred <<std::endl;


  // P derivative
  // Increment p 
  dt = 0.0001;
  LaplacianMomenta.M.ImportGauge(U);
  LatticeGaugeField UdSdP(&Grid);
  LaplacianMomenta.DerivativeP(UdSdP);


  LaplacianMomenta.Mom += dt*P;
  
  Sprime    = LaplacianMomenta.MomentaAction();

  // Prediciton

  dS = zero;
   for(int mu=0;mu<Nd;mu++){
    auto dSdPmu = PeekIndex<LorentzIndex>(UdSdP,mu);
    auto Pmu = PeekIndex<LorentzIndex>(P,mu);
    // Update gauge action density
    // 
    // dMom/dt = P
    // so dSdt = trace( dPdt dSdP) = trace( P dSdP ) 

    dS = dS + trace(Pmu*dSdPmu)*dt*2.0;
  } 
  dSpred    = sum(dS);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "pred dS "<< dSpred <<std::endl; 

  assert( fabs(real(Sprime-S-dSpred)) < 1.0 ) ;

  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
