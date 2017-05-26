/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/hmc/Test_WilsonFlow.cc

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  GridLogLayout();

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size, simd_layout, mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size, simd_layout, mpi_layout);

  std::vector<int> seeds({1, 2, 3, 4, 5});
  GridSerialRNG sRNG;
  GridParallelRNG pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(&Grid), Uflow(&Grid);
  SU<Nc>::HotConfiguration(pRNG, Umu);
  CheckpointerParameters CPPar("ckpoint_lat", "ckpoint_rng");
  BinaryHmcCheckpointer<PeriodicGimplR> CPBin(CPPar);

  CPBin.CheckpointRestore(3000, Umu, sRNG, pRNG);

  std::cout << std::setprecision(15);
  std::cout << GridLogMessage << "Plaquette: "
    << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << std::endl;

  WilsonFlow<PeriodicGimplR> WF(200, 0.01, 50);

  WF.smear(Uflow, Umu);

  RealD WFlow_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Uflow);
  RealD WFlow_TC   = WilsonLoops<PeriodicGimplR>::TopologicalCharge(Uflow);
  std::cout << GridLogMessage << "Plaquette         : "<< WFlow_plaq << std::endl;
  std::cout << GridLogMessage << "TopologicalCharge : "<< WFlow_TC   << std::endl;

  std::cout<< GridLogMessage << " Admissibility check:\n";
  const double sp_adm = 0.067;                // admissible threshold
  const double pl_adm = 1.0-sp_adm/Nc;
  std::cout << GridLogMessage << "   (pl_adm =" << pl_adm << ")\n";

  // Need min and reduce min for this function
  //double sp_max = NC_*(1.0-stpl.plaq_min(U,pl_adm));
  double sp_ave = Nc*(1.0-WFlow_plaq);

  //std::cout<< GridLogMessage << "   sp_max = "        << sp_max <<"\n";
  std::cout<< GridLogMessage << "   sp_ave = "        << sp_ave <<"\n";
  std::cout<< GridLogMessage << "   (sp_admissible = "<< sp_adm <<")\n";
  //std::cout<< GridLogMessage << "   sp_admissible - sp_max = "<<sp_adm-sp_max <<"\n";
  std::cout<< GridLogMessage << "   sp_admissible - sp_ave = "<<sp_adm-sp_ave <<"\n";

  Grid_finalize();
}  // main
