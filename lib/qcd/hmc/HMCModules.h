/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/GenericHmcRunner.h

Copyright (C) 2015
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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_HMC_MODULES
#define GRID_HMC_MODULES


#include "HMC_GridModules.h"

namespace Grid {
namespace QCD {

////////////////////////////////////////////////////////////////////
class RNGModuleParameters: Serializable {

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(RNGModuleParameters,
  std::string, serial_seeds,
  std::string, parallel_seeds,);
  std::vector<int> SerialSeed;
  std::vector<int> ParallelSeed;

  RNGModuleParameters(const std::vector<int> S = std::vector<int>(),
                      const std::vector<int> P = std::vector<int>())
      : SerialSeed(S), ParallelSeed(P) {}


  template <class ReaderClass >
  RNGModuleParameters(Reader<ReaderClass>& Reader){
    read(Reader, "RandomNumberGenerator", *this); 
    SerialSeed = strToVec<int>(serial_seeds);
    ParallelSeed = strToVec<int>(parallel_seeds);
  }
  
};

// Random number generators module
class RNGModule{
   GridSerialRNG sRNG_;
   std::unique_ptr<GridParallelRNG> pRNG_;
   RNGModuleParameters Params_;

public:

  RNGModule(){};

  void set_pRNG(GridParallelRNG* pRNG){
    pRNG_.reset(pRNG);
  }

  void set_RNGSeeds(RNGModuleParameters& Params) {
    Params_ = Params;
  }

  GridSerialRNG& get_sRNG() { return sRNG_; }
  GridParallelRNG& get_pRNG() { return *pRNG_.get(); }

  void seed() {
    if (Params_.SerialSeed.size() == 0 && Params_.ParallelSeed.size() == 0) {
      std::cout << "Seeds not initialized" << std::endl;
      exit(1);
    }
    sRNG_.SeedFixedIntegers(Params_.SerialSeed);
    pRNG_->SeedFixedIntegers(Params_.ParallelSeed);
  }
};


/*
///////////////////////////////////////////////////////////////////
/// Smearing module
template <class ImplementationPolicy>
class SmearingModule{
   virtual void get_smearing();
};

template <class ImplementationPolicy>
class StoutSmearingModule: public SmearingModule<ImplementationPolicy>{
   SmearedConfiguration<ImplementationPolicy> SmearingPolicy;
};

*/



}  // namespace QCD
}  // namespace Grid

#endif  // GRID_HMC_MODULES
