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

  LatticeGaugeField U(&Grid);

  SU3::HotConfiguration(pRNG,U);
  
  double beta = 1.0;
  double c1   = 0.331;

  //GparityPlaqPlusRectangleActionR Action(beta,c1);
  ConjugateWilsonGaugeActionR Action(beta);
  //WilsonGaugeActionR Action(beta);

  ComplexD S    = Action.S(U);

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(&Grid);

  Action.deriv(U,UdSdU);

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 0.0001;

  LatticeColourMatrix mommu(&Grid); 
  LatticeColourMatrix forcemu(&Grid); 
  LatticeGaugeField mom(&Grid); 
  LatticeGaugeField Uprime(&Grid); 

  for(int mu=0;mu<Nd;mu++){

    SU3::GaussianLieAlgebraMatrix(pRNG, mommu); // Traceless antihermitian momentum; gaussian in lie alg

    PokeIndex<LorentzIndex>(mom,mommu,mu);

    // fourth order exponential approx
    parallel_for(auto i=mom.begin();i<mom.end();i++){ // exp(pmu dt) * Umu
      Uprime[i](mu) = U[i](mu) + mom[i](mu)*U[i](mu)*dt ;
    }
  }

  ComplexD Sprime    = Action.S(Uprime);

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

    dS = dS - trace(mommu*UdSdUmu)*dt*2.0;

  }
  Complex dSpred    = sum(dS);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "pred dS "<< dSpred <<std::endl;

  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
