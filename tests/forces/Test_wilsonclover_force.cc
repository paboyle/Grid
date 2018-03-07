/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/Test_wilson_force.cc

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
using namespace Grid::QCD;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  std::vector<int> latt_size = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  std::vector<int> mpi_layout = GridDefaultMpi();

  GridCartesian Grid(latt_size, simd_layout, mpi_layout);
  GridRedBlackCartesian RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  std::vector<int> seeds({1, 2, 30, 50});

  GridParallelRNG pRNG(&Grid);

  std::vector<int> vrand(4);
  std::srand(std::time(0));
  std::generate(vrand.begin(), vrand.end(), std::rand);
  std::cout << GridLogMessage << vrand << std::endl;
  pRNG.SeedFixedIntegers(vrand);
  //pRNG.SeedFixedIntegers(seeds);

  LatticeFermion phi(&Grid);
  gaussian(pRNG, phi);
  LatticeFermion Mphi(&Grid);
  LatticeFermion MphiPrime(&Grid);

  LatticeGaugeField U(&Grid);

  std::vector<int> site = {0, 0, 0, 0};
  SU3::HotConfiguration(pRNG, U);
  //SU3::ColdConfiguration(pRNG, U);// Clover term zero

  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass = 0.1;
  Real csw = 1.0;
  WilsonCloverFermionR Dw(U, Grid, RBGrid, mass, csw, csw);
  Dw.ImportGauge(U);
  Dw.M(phi, Mphi);
  ComplexD S = innerProduct(Mphi, Mphi); // Action : pdag MdagM p

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(&Grid);
  LatticeGaugeField tmp(&Grid);

  ////////////////////////////////////////////
  Dw.MDeriv(tmp, Mphi, phi, DaggerNo);
  UdSdU = tmp;
  Dw.MDeriv(tmp, phi, Mphi, DaggerYes);
  UdSdU += tmp;
  /////////////////////////////////////////////

  ////////////////////////////////////
  // Modify the gauge field a little
  ////////////////////////////////////
  RealD dt = 0.00005;
  RealD Hmom = 0.0;
  RealD Hmomprime = 0.0;
  RealD Hmompp = 0.0;
  LatticeColourMatrix mommu(&Grid);
  LatticeColourMatrix forcemu(&Grid);
  LatticeGaugeField mom(&Grid);
  LatticeGaugeField Uprime(&Grid);

  for (int mu = 0; mu < Nd; mu++)
  {
    // Traceless antihermitian momentum; gaussian in lie alg
    SU3::GaussianFundamentalLieAlgebraMatrix(pRNG, mommu);
    Hmom -= real(sum(trace(mommu * mommu)));
    PokeIndex<LorentzIndex>(mom, mommu, mu);

    parallel_for(int ss = 0; ss < mom._grid->oSites(); ss++)
    {
      Uprime[ss]._internal[mu] = ProjectOnGroup(Exponentiate(mom[ss]._internal[mu], dt, 12) * U[ss]._internal[mu]);
    }
  }

  std::cout << GridLogMessage << "Initial mom hamiltonian is " << Hmom << std::endl;

  // New action
  Dw.ImportGauge(Uprime);
  Dw.M(phi, MphiPrime);
  ComplexD Sprime = innerProduct(MphiPrime, MphiPrime);

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////

  LatticeComplex dS(&Grid);
  dS = zero;
  LatticeComplex dSmom(&Grid);
  dSmom = zero;
  LatticeComplex dSmom2(&Grid);
  dSmom2 = zero;

  for (int mu = 0; mu < Nd; mu++)
  {
    mommu = PeekIndex<LorentzIndex>(UdSdU, mu); // P_mu =
    mommu = Ta(mommu) * 2.0;                    // Mom = (P_mu - P_mu^dag) - trace(P_mu - P_mu^dag)
    PokeIndex<LorentzIndex>(UdSdU, mommu, mu);  // UdSdU_mu = Mom
  }

  std::cout << GridLogMessage << "Antihermiticity tests" << std::endl;
  for (int mu = 0; mu < Nd; mu++)
  {
    mommu = PeekIndex<LorentzIndex>(mom, mu);
    std::cout << GridLogMessage << " Mommu  " << norm2(mommu) << std::endl;
    mommu = mommu + adj(mommu);
    std::cout << GridLogMessage << " Mommu + Mommudag " << norm2(mommu) << std::endl;
    mommu = PeekIndex<LorentzIndex>(UdSdU, mu);
    std::cout << GridLogMessage << " dsdumu  " << norm2(mommu) << std::endl;
    mommu = mommu + adj(mommu);
    std::cout << GridLogMessage << " dsdumu + dag  " << norm2(mommu) << std::endl;
    std::cout << "" << std::endl;
  }
  /////////////////////////////////////////////////////

  for (int mu = 0; mu < Nd; mu++)
  {
    forcemu = PeekIndex<LorentzIndex>(UdSdU, mu);
    mommu = PeekIndex<LorentzIndex>(mom, mu);

    // Update PF action density
    dS = dS + trace(mommu * forcemu) * dt;

    dSmom = dSmom - trace(mommu * forcemu) * dt;
    dSmom2 = dSmom2 - trace(forcemu * forcemu) * (0.25 * dt * dt);

    // Update mom action density
    mommu = mommu + forcemu * (dt * 0.5);

    Hmomprime -= real(sum(trace(mommu * mommu)));
  }

  ComplexD dSpred = sum(dS);
  ComplexD dSm = sum(dSmom);
  ComplexD dSm2 = sum(dSmom2);

  std::cout << GridLogMessage << "Initial mom hamiltonian is " << Hmom << std::endl;
  std::cout << GridLogMessage << "Final   mom hamiltonian is " << Hmomprime << std::endl;
  std::cout << GridLogMessage << "Delta   mom hamiltonian is " << Hmomprime - Hmom << std::endl;

  std::cout << GridLogMessage << " S      " << S << std::endl;
  std::cout << GridLogMessage << " Sprime " << Sprime << std::endl;
  std::cout << GridLogMessage << "dS (S' - S)          :" << Sprime - S << std::endl;
  std::cout << GridLogMessage << "predict dS (force)   :" << dSpred << std::endl;
  std::cout << GridLogMessage << "dSm " << dSm << std::endl;
  std::cout << GridLogMessage << "dSm2" << dSm2 << std::endl;

  std::cout << GridLogMessage << "Total dS    " << Hmomprime - Hmom + Sprime - S << std::endl;

  assert(fabs(real(Sprime - S - dSpred)) < 1.0);

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
}
