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
struct RNGModuleParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(RNGModuleParameters,
  std::string, serial_seeds,
  std::string, parallel_seeds,);

  std::vector<int> getSerialSeeds(){return strToVec<int>(serial_seeds);}
  std::vector<int> getParallelSeeds(){return strToVec<int>(parallel_seeds);}

  RNGModuleParameters(): serial_seeds("1"), parallel_seeds("1"){}

  template <class ReaderClass >
  RNGModuleParameters(Reader<ReaderClass>& Reader){
    read(Reader, "RandomNumberGenerator", *this); 
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
    auto SerialSeeds   = Params_.getSerialSeeds();
    auto ParallelSeeds = Params_.getParallelSeeds();
    if (SerialSeeds.size() == 0 && ParallelSeeds.size() == 0) {
      std::cout << GridLogError << "Seeds not initialized" << std::endl;
      exit(1);
    }
    sRNG_.SeedFixedIntegers(SerialSeeds);
    pRNG_->SeedFixedIntegers(ParallelSeeds);
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
