    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_gp_rect_force.cc

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

 

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

  LatticeGaugeField U(&Grid);

  SU<Nc>::HotConfiguration(pRNG,U);
  
  double beta = 1.0;
  double c1   = 0.331;

  const int nu = 1;
  std::vector<int> twists(Nd,0);
  twists[nu] = 1;
  ConjugateGimplD::setDirections(twists);
  ConjugatePlaqPlusRectangleActionR Action(beta,c1);
  //ConjugateWilsonGaugeActionR Action(beta);
  //WilsonGaugeActionR Action(beta);

  ComplexD S    = Action.S(U);

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(&Grid);

  Action.deriv(U,UdSdU);

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 0.01;

  LatticeColourMatrix mommu(&Grid); 
  LatticeColourMatrix forcemu(&Grid); 
  LatticeGaugeField mom(&Grid); 
  LatticeGaugeField Uprime(&Grid); 

  for(int mu=0;mu<Nd;mu++){

    SU<Nc>::GaussianFundamentalLieAlgebraMatrix(pRNG, mommu); // Traceless antihermitian momentum; gaussian in lie alg

    PokeIndex<LorentzIndex>(mom,mommu,mu);

    // fourth order exponential approx
    autoView( mom_v, mom, CpuRead);
    autoView(Uprime_v, Uprime, CpuWrite);
    autoView( U_v , U, CpuRead);
    thread_foreach(i,mom_v,{ // exp(pmu dt) * Umu
      Uprime_v[i](mu) = U_v[i](mu) + mom_v[i](mu)*U_v[i](mu)*dt ;
    });
  }

  ComplexD Sprime    = Action.S(Uprime);

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////

  LatticeComplex dS(&Grid); dS = Zero();

  for(int mu=0;mu<Nd;mu++){

    auto UdSdUmu = PeekIndex<LorentzIndex>(UdSdU,mu);
         mommu   = PeekIndex<LorentzIndex>(mom,mu);

    // Update gauge action density
    // U = exp(p dt) U
    // dU/dt = p U
    // so dSdt = trace( dUdt dSdU) = trace( p UdSdUmu ) 

    dS = dS - trace(mommu*UdSdUmu)*dt*2.0;

  }
  ComplexD dSpred    = sum(dS);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "pred dS "<< dSpred <<std::endl;
  assert( fabs(real(Sprime-S-dSpred)) < 1.0e-1 ) ;
  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
