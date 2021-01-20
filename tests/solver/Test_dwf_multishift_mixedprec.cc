    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_multishift_mixedprec.cc

    Copyright (C) 2015

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

using namespace Grid;

template<typename SpeciesD, typename SpeciesF, typename GaugeStatisticsType>
void run_test(int argc, char ** argv, const typename SpeciesD::ImplParams &params){
  const int Ls = 16;
  GridCartesian* UGrid_d = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid_d = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_d);
  GridCartesian* FGrid_d = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid_d);
  GridRedBlackCartesian* FrbGrid_d = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid_d);

  GridCartesian* UGrid_f = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian* FGrid_f = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid_f);
  GridRedBlackCartesian* FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid_f);

  typedef typename SpeciesD::FermionField FermionFieldD;
  typedef typename SpeciesF::FermionField FermionFieldF;
  
  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid_d);
  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid_d);
  RNG4.SeedFixedIntegers(seeds4);

  FermionFieldD src_d(FGrid_d);
  random(RNG5, src_d);

  LatticeGaugeFieldD Umu_d(UGrid_d);

  //CPS-created G-parity ensembles have a factor of 2 error in the plaquette that causes the read to fail unless we workaround it
  bool gparity_plaquette_fix = false;
  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--gparity_plaquette_fix"){
      gparity_plaquette_fix=true;
      break;
    }
  }

  bool cfg_loaded=false;
  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--load_config"){
      assert(i != argc-1);
      std::string file = argv[i+1];
      NerscIO io;
      FieldMetaData metadata;

      if(gparity_plaquette_fix) NerscIO::exitOnReadPlaquetteMismatch() = false;

      io.readConfiguration<GaugeStatisticsType>(Umu_d, metadata, file);

      if(gparity_plaquette_fix){
	metadata.plaquette *= 2.; //correct header value

	//Get the true plaquette
	FieldMetaData tmp;
	GaugeStatisticsType gs; gs(Umu_d, tmp);
	
	std::cout << "After correction: plaqs " << tmp.plaquette << " " << metadata.plaquette << std::endl;
	assert(fabs(tmp.plaquette -metadata.plaquette ) < 1.0e-5 );
      }

      cfg_loaded=true;
      break;
    }
  }

  if(!cfg_loaded)
    SU<Nc>::HotConfiguration(RNG4, Umu_d);

  LatticeGaugeFieldF Umu_f(UGrid_f);
  precisionChange(Umu_f, Umu_d);

  std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt() << "   Ls: " << Ls << std::endl;

  RealD mass = 0.01;
  RealD M5 = 1.8;
  SpeciesD Ddwf_d(Umu_d, *FGrid_d, *FrbGrid_d, *UGrid_d, *UrbGrid_d, mass, M5, params);
  SpeciesF Ddwf_f(Umu_f, *FGrid_f, *FrbGrid_f, *UGrid_f, *UrbGrid_f, mass, M5, params);

  FermionFieldD src_o_d(FrbGrid_d);
  pickCheckerboard(Odd, src_o_d, src_d);

  SchurDiagMooeeOperator<SpeciesD, FermionFieldD> HermOpEO_d(Ddwf_d);
  SchurDiagMooeeOperator<SpeciesF, FermionFieldF> HermOpEO_f(Ddwf_f);

  AlgRemez remez(1e-4, 64, 50);
  int order = 15;
  remez.generateApprox(order, 1, 2); //sqrt

  MultiShiftFunction shifts(remez, 1e-10, false);

  int relup_freq = 50;
  double t1=usecond();
  ConjugateGradientMultiShiftMixedPrec<FermionFieldD,FermionFieldF> mcg(10000, shifts, FrbGrid_f, HermOpEO_f, relup_freq);

  std::vector<FermionFieldD> results_o_d(order, FrbGrid_d);
  mcg(HermOpEO_d, src_o_d, results_o_d);
  double t2=usecond();

  std::cout<<GridLogMessage << "Test: Total usec    =   "<< (t2-t1)<<std::endl;
}





int main (int argc, char ** argv)
{
  Grid_init(&argc, &argv);

  bool gparity = false;
  int gpdir;

  for(int i=1;i<argc;i++){
    std::string arg(argv[i]);
    if(arg == "--Gparity"){
      assert(i!=argc-1);
      gpdir = std::stoi(argv[i+1]);
      assert(gpdir >= 0 && gpdir <= 2); //spatial!
      gparity = true;
    }
  }
  if(gparity){
    std::cout << "Running test with G-parity BCs in " << gpdir << " direction" << std::endl;
    GparityWilsonImplParams params;
    params.twists[gpdir] = 1;
    
    std::vector<int> conj_dirs(Nd,0);
    conj_dirs[gpdir] = 1;
    ConjugateGimplD::setDirections(conj_dirs);

    run_test<GparityDomainWallFermionD, GparityDomainWallFermionF, ConjugateGaugeStatistics>(argc,argv,params);
  }else{
    std::cout << "Running test with periodic BCs" << std::endl;
    WilsonImplParams params;
    run_test<DomainWallFermionD, DomainWallFermionF, PeriodicGaugeStatistics>(argc,argv,params);
  }

  Grid_finalize();
}
