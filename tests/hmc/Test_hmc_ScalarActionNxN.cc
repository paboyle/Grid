/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2016

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>
namespace Grid {
class ScalarActionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ScalarActionParameters,
    double, mass_squared,
    double, lambda,
    double, g);

  ScalarActionParameters() = default;

    template <class ReaderClass >
  ScalarActionParameters(Reader<ReaderClass>& Reader){
    read(Reader, "ScalarAction", *this);
  }

};
}

using namespace Grid;
using namespace Grid::QCD;

template <class Impl>
class MagMeas : public HmcObservable<typename Impl::Field> {
public:
  typedef typename Impl::Field Field;
  typedef typename Impl::Simd::scalar_type Trace;
  
  void TrajectoryComplete(int traj,
                          Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {
    
    int def_prec = std::cout.precision();
    
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10 + 1);
    std::cout << GridLogMessage
              << "m= " << TensorRemove(trace(sum(U))) << std::endl;
    std::cout << GridLogMessage
              << "m^2= " << TensorRemove(trace(sum(U)*sum(U))) << std::endl;
    std::cout << GridLogMessage
    << "phi^2= " << TensorRemove(sum(trace(U*U))) << std::endl;
    std::cout.precision(def_prec);
    
  }
private:
  
};

template <class Impl>
class MagMod: public ObservableModule<MagMeas<Impl>, NoParameters>{
  typedef ObservableModule<MagMeas<Impl>, NoParameters> ObsBase;
  using ObsBase::ObsBase; // for constructors
  
  // acquire resource
  virtual void initialize(){
    this->ObservablePtr.reset(new MagMeas<Impl>());
  }
public:
  MagMod(): ObsBase(NoParameters()){}
};

int main(int argc, char **argv) {
  typedef Grid::JSONReader       Serialiser;
  
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  // Typedefs to simplify notation
  constexpr int Ncolours    = 2;
  constexpr int Ndimensions = 3;
  typedef ScalarNxNAdjGenericHMCRunner<Ncolours> HMCWrapper;  // Uses the default minimum norm, real scalar fields
  typedef ScalarAdjActionR<Ncolours, Ndimensions> ScalarAction;
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;
  TheHMC.ReadCommandLine(argc, argv);

  if (TheHMC.ParameterFile.empty()){
    std::cout << "Input file not specified."
              << "Use --ParameterFile option in the command line.\nAborting" 
              << std::endl;
    exit(1);
  }
  Serialiser Reader(TheHMC.ParameterFile);

  // Grid from the command line
  GridModule ScalarGrid;
  if (GridDefaultLatt().size() != Ndimensions){
    std::cout << "Incorrect dimension of the grid\n. Expected dim="<< Ndimensions << std::endl;
    exit(1);
  }
  if (GridDefaultMpi().size() != Ndimensions){
    std::cout << "Incorrect dimension of the mpi grid\n. Expected dim="<< Ndimensions << std::endl;
    exit(1);
  }
  ScalarGrid.set_full(new GridCartesian(GridDefaultLatt(),GridDefaultSimd(Ndimensions, vComplex::Nsimd()),GridDefaultMpi()));
  ScalarGrid.set_rb(new GridRedBlackCartesian(ScalarGrid.get_full()));
  TheHMC.Resources.AddGrid("scalar", ScalarGrid);
  std::cout << "Lattice size : " << GridDefaultLatt() << std::endl;
  
  ScalarActionParameters SPar(Reader);

  // Checkpointer definition
  CheckpointerParameters CPparams(Reader);
  //TheHMC.Resources.LoadBinaryCheckpointer(CPparams);
  TheHMC.Resources.LoadScidacCheckpointer(CPparams, SPar);

  RNGModuleParameters RNGpar(Reader);
  TheHMC.Resources.SetRNGSeeds(RNGpar);
  
  // Construct observables
  typedef MagMod<HMCWrapper::ImplPolicy> MagObs;
  TheHMC.Resources.AddObservable<MagObs>();
  
  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation

  // Scalar action in adjoint representation
  ScalarAction Saction(SPar.mass_squared, SPar.lambda, SPar.g);

  // Collect actions
  ActionLevel<ScalarAction::Field, ScalarNxNMatrixFields<Ncolours>> Level1(1);
  Level1.push_back(&Saction);
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////
  TheHMC.Parameters.initialize(Reader);

  TheHMC.Run();

  Grid_finalize();
}  // main

/* Examples for input files

JSON

{
    "Checkpointer": {
    "config_prefix": "ckpoint_scalar_lat",
    "rng_prefix": "ckpoint_scalar_rng",
    "saveInterval": 1,
    "format": "IEEE64BIG"
    },
    "RandomNumberGenerator": {
    "serial_seeds": "1 2 3 4 6",
    "parallel_seeds": "6 7 8 9 11"
    },
    "ScalarAction":{
      "mass_squared": 0.5,
      "lambda": 0.1
    },
    "HMC":{
    "StartTrajectory": 0,
    "Trajectories": 100,
    "MetropolisTest": true,
    "NoMetropolisUntil": 10,
    "StartingType": "HotStart",
    "MD":{
        "name": "MinimumNorm2",
	      "MDsteps": 15,
	      "trajL": 2.0
	    }
    }
}


XML example not provided yet

*/
