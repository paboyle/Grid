/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/hmc/Test_WilsonFlow_adaptive.cc

Copyright (C) 2017

Author: Christopher Kelly <ckelly@bnl.gov>

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

using namespace Grid;

//Linearly interpolate between two nearest times
RealD interpolate(const RealD t_int, const std::vector<std::pair<RealD,RealD> > &data){
  RealD tdiff1=1e32; int t1_idx=-1;
  RealD tdiff2=1e32; int t2_idx=-1;

  for(int i=0;i<data.size();i++){
    RealD diff = fabs(data[i].first-t_int);
    //std::cout << "targ " << t_int << " cur " << data[i].first << " diff " << diff << " best diff1 " << tdiff1 << " diff2 " << tdiff2 << std::endl;

    if(diff < tdiff1){ 
      if(tdiff1 < tdiff2){ //swap out tdiff2
	tdiff2 = tdiff1; t2_idx = t1_idx;
      }
      tdiff1 = diff; t1_idx = i; 
    }
    else if(diff < tdiff2){ tdiff2 = diff; t2_idx = i; }
  }
  assert(t1_idx != -1 && t2_idx != -1);
  
  RealD t2 = data[t2_idx].first,  v2 = data[t2_idx].second;
  RealD t1 = data[t1_idx].first,  v1 = data[t1_idx].second;
  
  //v = a + bt
  //v2-v1 = b(t2-t1)
  RealD b = (v2-v1)/(t2-t1);
  RealD a = v1 - b*t1;
  RealD vout = a + b*t_int;

  //std::cout << "Interpolate to " << t_int << " two closest points " << t1  << " " << t2 
  //<< " with values " << v1 << " "<< v2 << " : got " << vout <<  std::endl;
  return vout;
}


int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size, simd_layout, mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  std::vector<int> seeds({1, 2, 3, 4, 5});
  GridSerialRNG sRNG;
  GridParallelRNG pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField U(&Grid);
  SU<Nc>::HotConfiguration(pRNG, U);

  int Nstep = 300;
  RealD epsilon = 0.01;
  RealD maxTau = Nstep*epsilon;
  RealD tolerance = 1e-4;

  for(int i=1;i<argc;i++){
    std::string sarg(argv[i]);
    if(sarg == "--tolerance"){
      std::stringstream ss; ss << argv[i+1]; ss >> tolerance;
    }
  }
  std::cout << "Adaptive smear tolerance " << tolerance << std::endl;

  //Setup iterative Wilson flow
  WilsonFlow<PeriodicGimplD> wflow(epsilon,Nstep);
  wflow.resetActions();    

  std::vector<std::pair<RealD, RealD> > meas_orig;

  wflow.addMeasurement(1, [&wflow,&meas_orig](int step, RealD t, const LatticeGaugeField &U){ 
      std::cout << GridLogMessage << "[WilsonFlow] Computing Cloverleaf energy density for step " << step << std::endl;
      meas_orig.push_back( {t, wflow.energyDensityCloverleaf(t,U)} );
    });

  //Setup adaptive Wilson flow
  WilsonFlowAdaptive<PeriodicGimplD> wflow_ad(epsilon,maxTau,tolerance);
  wflow_ad.resetActions();    

  std::vector<std::pair<RealD, RealD> > meas_adaptive;

  wflow_ad.addMeasurement(1, [&wflow_ad,&meas_adaptive](int step, RealD t, const LatticeGaugeField &U){ 
      std::cout << GridLogMessage << "[WilsonFlow] Computing Cloverleaf energy density for step " << step << std::endl;
      meas_adaptive.push_back( {t, wflow_ad.energyDensityCloverleaf(t,U)} );
    });

  //Run
  LatticeGaugeFieldD Vtmp(U.Grid());
  wflow.smear(Vtmp, U); //basic smear

  Vtmp = Zero();
  wflow_ad.smear(Vtmp, U);

  //Output values for plotting
  {
    std::ofstream out("wflow_t2E_orig.dat");
    out.precision(16);
    for(auto const &e: meas_orig){
      out << e.first << " " << e.second << std::endl;
    }
  }
  {
    std::ofstream out("wflow_t2E_adaptive.dat");
    out.precision(16);
    for(auto const &e: meas_adaptive){
      out << e.first << " " << e.second << std::endl;
    }
  }
  
  //Compare at times available with adaptive smearing
  for(int i=0;i<meas_adaptive.size();i++){
    RealD t = meas_adaptive[i].first;
    RealD v_adaptive = meas_adaptive[i].second;
    RealD v_orig = interpolate(t,  meas_orig); //should be very precise due to fine timestep
    std::cout << t << " orig: " << v_orig << " adaptive: " << v_adaptive << " reldiff: " << (v_adaptive-v_orig)/v_orig << std::endl;
  }
  
  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
}
