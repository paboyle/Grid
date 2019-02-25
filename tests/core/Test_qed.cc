/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: tests/core/Test_qed.cc

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: James Harrison <J.Harrison@soton.ac.uk>

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

#include <Grid/Grid.h>

using namespace Grid;
using namespace QCD;

typedef PeriodicGaugeImpl<QedGImplR>  QedPeriodicGImplR;
typedef PhotonR::GaugeField           EmField;
typedef PhotonR::GaugeLinkField       EmComp;

const int NCONFIGS = 20;
const int NWILSON  = 10;

int main(int argc, char *argv[])
{
  // initialization
  Grid_init(&argc, &argv);
  std::cout << GridLogMessage << "Grid initialized" << std::endl;
  
  // QED stuff
  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4, vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian    grid(latt_size,simd_layout,mpi_layout);
  GridParallelRNG  pRNG(&grid);
  PhotonR          photon(&grid, PhotonR::Gauge::coulomb, PhotonR::ZmScheme::qedL);
  EmField          a(&grid);
  EmField          expA(&grid);

  Complex imag_unit(0, 1);

  Real wlA;
  std::vector<Real> logWlAvg(NWILSON, 0.0), logWlTime(NWILSON, 0.0), logWlSpace(NWILSON, 0.0);

  pRNG.SeedFixedIntegers({1, 2, 3, 4});

  std::cout << GridLogMessage << "Wilson loop calculation beginning" << std::endl;
  for(int ic = 0; ic < NCONFIGS; ic++){
      std::cout << GridLogMessage << "Configuration " << ic <<std::endl;
      photon.StochasticField(a, pRNG);

      // Exponentiate photon field
      expA = exp(imag_unit*a);

      // Calculate zero-modes
      std::vector<EmField::vector_object::scalar_object> zm;

      std::cout << GridLogMessage << "Total zero-mode norm 2 " 
                << std::sqrt(norm2(sum(a))) << std::endl;

      std::cout << GridLogMessage << "Spatial zero-mode norm 2" << std::endl;
      sliceSum(a, zm, grid.Nd() - 1);
      for (unsigned int t = 0; t < latt_size.back(); ++t)
      {
        std::cout << GridLogMessage << "t = " << t << " " << std::sqrt(norm2(zm[t])) << std::endl;
      }

      // Calculate divergence
      EmComp diva(&grid), amu(&grid);

      diva = zero;
      for (unsigned int mu = 0; mu < grid.Nd(); ++mu)
      {
        amu   = peekLorentz(a, mu);
        diva += amu - Cshift(amu, mu, -1);
        if (mu == grid.Nd() - 2)
        {
          std::cout << GridLogMessage << "Spatial divergence norm 2 " << std::sqrt(norm2(diva)) << std::endl;
        }
      }
      std::cout << GridLogMessage << "Total divergence norm 2 " << std::sqrt(norm2(diva)) << std::endl;

      // Calculate Wilson loops
      for(int iw=1; iw<=NWILSON; iw++){
          wlA = WilsonLoops<QedPeriodicGImplR>::avgWilsonLoop(expA, iw, iw) * 3;
          logWlAvg[iw-1] -= 2*log(wlA);
          wlA = WilsonLoops<QedPeriodicGImplR>::avgTimelikeWilsonLoop(expA, iw, iw) * 3;
          logWlTime[iw-1] -= 2*log(wlA);
          wlA = WilsonLoops<QedPeriodicGImplR>::avgSpatialWilsonLoop(expA, iw, iw) * 3;
          logWlSpace[iw-1] -= 2*log(wlA);
      }
  }
  std::cout << GridLogMessage << "Wilson loop calculation completed" << std::endl;
  
  // Calculate Wilson loops
  // From A. Portelli's PhD thesis:
  // size  -2*log(W)
  // 1     0.500000000(1)
  // 2     1.369311535(1) 
  // 3     2.305193057(1) 
  // 4     3.261483854(1) 
  // 5     4.228829967(1) 
  // 6     5.203604529(1) 
  // 7     6.183728249(1) 
  // 8     7.167859805(1) 
  // 9     8.155091868(1) 
  // 10    9.144788116(1)

  for(int iw=1; iw<=10; iw++){
      std::cout << GridLogMessage << iw << 'x' << iw << " Wilson loop" << std::endl;
      std::cout << GridLogMessage << "-2*log(W) average: " << logWlAvg[iw-1]/NCONFIGS << std::endl;
      std::cout << GridLogMessage << "-2*log(W) timelike: " << logWlTime[iw-1]/NCONFIGS << std::endl;
      std::cout << GridLogMessage << "-2*log(W) spatial: " << logWlSpace[iw-1]/NCONFIGS << std::endl;
  }

  // epilogue
  std::cout << GridLogMessage << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
