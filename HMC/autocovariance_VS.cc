/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/hmc/Test_WilsonFlow.cc

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Shuhei Yamamoto <syamamoto@bnl.gov>

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
#include <string>
#include "ACC.hpp"


int main(int argc, char **argv) {
  using namespace Grid;
  
  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);

  typedef typename PeriodicGimplR::ComplexField ComplexField;
  
  typedef Grid::XmlReader Serialiser;
  Serialiser              Reader("input_ACF.xml", false, "root");
  WFParameters            WFPar(Reader);
  ConfParameters          CPar(Reader);
  ACFParameters           APar(Reader);

  std::string fname;
  std::string file_path  = WFPar.path;
  int W = APar.MScut, tau = WFPar.tau; // W: data[t>W] is not used in Madras-Sokal approximation of the error
  ComplexD coeff0, coeff1;
  int total_configs = CPar.EndConfiguration - CPar.StartConfiguration + 1;
  int T             = total_configs/APar.MDtime_div_fac; // #processed configs per bin
  int arr_size      = T*APar.MDtime_div_fac;  
  if ( APar.MDtime_div_fac == 1 ) assert(APar.R < 0 || W >= 100 );

  ComplexField A0(&Grid), A1(&Grid), one(&Grid); one=ComplexField::scalar_type(1.0,0.0);
  std::vector<ComplexField> F(total_configs,&Grid), G(total_configs,&Grid), G2(total_configs,&Grid);

  
#if 0 // DEBUG: Test to see Cshift shifts the lattice cyclically
  Coordinate latt_size0({4,4,4,4});
  GridCartesian Grid0(latt_size0, simd_layout, mpi_layout);
  ComplexField one0(&Grid0); one0=ComplexField::scalar_type(1.0,0.0);
  LatticeInteger coor(&Grid0);
  ComplexField tmp(&Grid0), filter(&Grid0), zero(&Grid0); filter = one0; zero=Zero();
  int d[4] = {1,0,0,0};
  for (int d=0; d<Nd; d++) {
    LatticeCoordinate(coor,d);
    filter= where(mod(coor,2)==Integer(0),filter,zero);
  }
  tmp = filter*one0;
  for(int mu=0; mu<Nd; mu++)
    tmp = Cshift(tmp, mu, d[mu]);
  std::cout<< GridLogMessage << "Sfhited Sparsed Field " << tmp << std::endl;
#endif
  
  std::cout << std::setprecision(15);

  ////////////////////////////////
  // Retrieve the fields
  ///////////////////////////////
  for (int i=0; i<total_configs; i++) {
    int conf = CPar.StartConfiguration + i;
    fname = file_path + WFPar.data_name + "_" + std::to_string(tau) + "_" + CPar.conf_prefix + ".";
    readFile(A0, fname + std::to_string(conf));
    F[i] = A0;
#if 1 // DEBUG
    RealD out = real(sum(F[i]));
    std::cout << GridLogMessage << "LVS " + WFPar.data_name + " (conf, tau, val):   " << " " << conf << " " << tau << " "
	      <<  out/Real(Grid.gSites()) << std::endl;
#endif
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute VS autocovariance function in two ways (<- corresponds to 2 diff estimation methods)     
  ////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i_bin=0, it=0; it<arr_size; i_bin++){
    for (int t=0; t<T; t++, it++){
	
      G[it] = Zero(); G2[it] = Zero();

      int n_src_configs = APar.isFullTimeAvg ? T-t : 1;
      for (int i=0; i<n_src_configs; i++){

	A0 = F[i_bin*T + i];
	A1 = F[i_bin*T + i + t];

	coeff0 = TensorRemove(sum(A0))/RealD(Grid.gSites());
	coeff1 = TensorRemove(sum(A1))/RealD(Grid.gSites());
	G2[it] = G2[it] + A0*A1 - (coeff0*coeff1)*one;
	  
	A0 = A0 - coeff0*one;
	A1 = A1 - coeff1*one;
	G[it] = G[it] + A0*A1;

      }

      G[it] = (1.0/((RealD) n_src_configs )) * G[it]; G2[it] = (1.0/((RealD) n_src_configs )) * G2[it];

#if 0 // DEBUG: Monitor variation of COV
      Coordinate scoor; for (int mu=0; mu < Nd; mu++) scoor[mu] = 0;
      RealD Gt0 = real(peekSite(G[it],scoor));
      coeff0 = real(sum(G[it]))/RealD(Grid.gSites());
      A0 = G[it] - coeff0*one;
      A0 = A0*A0;
      std::cout << GridLogMessage << "LVS COV at origin for " + WFPar.data_name + " " << t << " " << Gt0 << " " << real(sum(A0))/RealD(Grid.gSites()) << std::endl;
#endif
      
    }// END(t): loop within a bin for MD time
  }// END(cong_s): loop over bins
	
  /////////////////////////////////
  // Data for Error Estimate
  ////////////////////////////////
  for(const int& bs : APar.space_block_sizes){
    ///////////////////////////////////////////////////////////////////////////////////////
    // Apply blocking on autocovariance function fields (block average or sparse sampling)
    ///////////////////////////////////////////////////////////////////////////////////////
    Coordinate clatt_size(Nd);
    for(int i=0;i<Nd;i++) clatt_size[i] = Grid.FullDimensions()[i]/bs;
    GridCartesian Coarse(clatt_size,simd_layout,mpi_layout);
    
    // blocking
    std::vector<ComplexField> G_B(arr_size,&Coarse), G2_B(arr_size,&Coarse);
    for (int i=0; i<arr_size; i++){
      blockSum(G_B[i],G[i]); blockSum(G2_B[i],G2[i]);
      G_B[i] = (1/RealD(std::pow(bs,Nd)))*G_B[i]; G2_B[i] = (1/RealD(std::pow(bs,Nd)))*G2_B[i];
    }
    
    // sparse sampling
    std::vector<ComplexField> G_s(arr_size,&Coarse), G2_s(arr_size,&Coarse);
    for (int i=0; i<arr_size; i++){
      
      LatticeInteger coor(&Grid);
      ComplexField tmp(&Grid), filter(&Grid), zero(&Grid); filter = one; zero = Zero();
      for (int d=0; d<Nd; d++) {
	LatticeCoordinate(coor,d);
	filter = where(mod(coor,bs)==Integer(0),filter,zero);
      }
      tmp = filter*G[i]; blockSum(G_s[i],tmp);
      tmp = filter*G2[i]; blockSum(G2_s[i],tmp);
    }

    /***********   Generate data  ****************************************************************
      No Binning => Use Madras-Sokal approximation or Master-field technique for error estimation
      Otherwise  => error estimate via sample variance by binning over MD time
    **********************************************************************************************/
    int n_bin = total_configs/T;
    if ( total_configs == T ){

      // Madras-Sokal Approximation
      if ( APar.R < 0 ){
	//   Valid: when t << T
	assert( T > W );
	MS_approx(&Coarse, G_B, T, W, bs, tau, "Blocked " + WFPar.data_name + " ACC", "LVS");
	MS_approx(&Coarse, G2_B, T, W, bs, tau, "Blocked " + WFPar.data_name + " ACC2", "LVS");
	
	MS_approx(&Coarse, G_s, T, W, bs, tau, "Sparsed " + WFPar.data_name + " ACC", "LVS");
	MS_approx(&Coarse, G2_s, T, W, bs, tau, "Sparsed " + WFPar.data_name + " ACC2", "LVS");
      }

      // Master-Field Approximation
      MF_approx(&Coarse, G_B, T, APar.R, bs, tau, "Blocked " + WFPar.data_name + " ACC", "LVS");
      //MF_approx(&Coarse, G2_B, T, APar.R, bs, tau, "Blocked " + WFPar.data_name + " ACC2", "LVS");
    }
    else {
      // Binning
      binning(&Coarse, G_B, T, n_bin, bs, tau, "Blocked " + WFPar.data_name + " ACC", "LVS");
      binning(&Coarse, G2_B, T, n_bin, bs, tau, "Blocked " + WFPar.data_name + " ACC2", "LVS");
      
      binning(&Coarse, G_s, T, n_bin, bs, tau, "Sparsed " + WFPar.data_name + " ACC", "LVS");
      binning(&Coarse, G2_s, T, n_bin, bs, tau, "Sparsed " + WFPar.data_name + " ACC2", "LVS");

      
      binning2(&Coarse, G_B, T, n_bin, bs, tau, "Blocked " + WFPar.data_name + " ACC", "LVS");
      binning2(&Coarse, G2_B, T, n_bin, bs, tau, "Blocked " + WFPar.data_name + " ACC2", "LVS");
      
      binning2(&Coarse, G_s, T, n_bin, bs, tau, "Sparsed " + WFPar.data_name + " ACC", "LVS");
      binning2(&Coarse, G2_s, T, n_bin, bs, tau, "Sparsed " + WFPar.data_name + " ACC2", "LVS");
    }
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
