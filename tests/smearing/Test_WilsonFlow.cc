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

namespace Grid{
  struct WFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(WFParameters,
            int, steps,
            double, step_size,
            int, meas_interval,
            double, maxTau); // for the adaptive algorithm
       

    template <class ReaderClass >
    WFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "WilsonFlow", *this);
    }

  };

  struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConfParameters,
           std::string, conf_prefix,
            std::string, rng_prefix,
				    int, StartConfiguration,
				    int, EndConfiguration,
            int, Skip);
  
    template <class ReaderClass >
    ConfParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Configurations", *this);
    }

  };
}

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  GridLogLayout();

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size, simd_layout, mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  std::vector<int> seeds({1, 2, 3, 4, 5});
  GridSerialRNG sRNG;
  GridParallelRNG pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(&Grid), Uflow(&Grid);
  SU<Nc>::HotConfiguration(pRNG, Umu);
  
  typedef Grid::JSONReader       Serialiser;
  Serialiser Reader("input.json");
  WFParameters WFPar(Reader);
  ConfParameters CPar(Reader);
  CheckpointerParameters CPPar(CPar.conf_prefix, CPar.rng_prefix);
  BinaryHmcCheckpointer<PeriodicGimplR> CPBin(CPPar);

  for (int conf = CPar.StartConfiguration; conf <= CPar.EndConfiguration; conf+= CPar.Skip){

  CPBin.CheckpointRestore(conf, Umu, sRNG, pRNG);

  std::cout << std::setprecision(15);
  std::cout << GridLogMessage << "Initial plaquette: "
    << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << std::endl;

  WilsonFlow<PeriodicGimplR> WF(WFPar.steps, WFPar.step_size, WFPar.meas_interval);

  WF.smear_adaptive(Uflow, Umu, WFPar.maxTau);

  RealD WFlow_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Uflow);
  RealD WFlow_TC   = WilsonLoops<PeriodicGimplR>::TopologicalCharge(Uflow);
  RealD WFlow_T0   = WF.energyDensityPlaquette(Uflow);
  std::cout << GridLogMessage << "Plaquette          "<< conf << "   " << WFlow_plaq << std::endl;
  std::cout << GridLogMessage << "T0                 "<< conf << "   " << WFlow_T0 << std::endl;
  std::cout << GridLogMessage << "TopologicalCharge  "<< conf << "   " << WFlow_TC   << std::endl;

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
  }
  Grid_finalize();
}  // main


/*
Input file example


JSON

{
    "WilsonFlow":{
	"steps": 200,
	"step_size": 0.01,
	"meas_interval": 50,
  "maxTau": 2.0
    },
    "Configurations":{
	"conf_prefix": "ckpoint_lat",
	"rng_prefix": "ckpoint_rng",
	"StartConfiguration": 3000,
	"EndConfiguration": 3000,
	"Skip": 5
    }
}


*/
