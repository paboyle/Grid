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

namespace Grid {
namespace QCD {

///////////////////////////////////////////////////
// Modules
class GridModule {
 public:
  GridCartesian* get_full() { return grid_.get(); }
  GridRedBlackCartesian* get_rb() { return rbgrid_.get(); }

  void set_full(GridCartesian* grid) { grid_.reset(grid); }
  void set_rb(GridRedBlackCartesian* rbgrid) { rbgrid_.reset(rbgrid); }

 protected:
  std::unique_ptr<GridCartesian> grid_;
  std::unique_ptr<GridRedBlackCartesian> rbgrid_;
};

// helpers
class GridFourDimModule : public GridModule {
 public:
  GridFourDimModule() {
    set_full(SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(4, vComplex::Nsimd()),
        GridDefaultMpi()));
    set_rb(SpaceTimeGrid::makeFourDimRedBlackGrid(grid_.get()));
  }

};

class RNGModule{
   // Random number generators
   GridSerialRNG sRNG_;
   std::unique_ptr<GridParallelRNG> pRNG_;
   std::vector<int> SerialSeed_;
   std::vector<int> ParallelSeed_;

public:
   void set_pRNG(GridParallelRNG* pRNG){
      pRNG_.reset(pRNG);
   }

   void set_RNGSeeds(const std::vector<int> S, const std::vector<int> P) {
    SerialSeed_   = S;
    ParallelSeed_ = P;
  }

  GridSerialRNG& get_sRNG(){return sRNG_;}
  GridParallelRNG& get_pRNG(){return *pRNG_.get();}
  void seed(){
    sRNG_.SeedFixedIntegers(SerialSeed_);
    pRNG_->SeedFixedIntegers(ParallelSeed_);
  }
};

/// Smearing module
template <class ImplementationPolicy>
class SmearingModule{
   virtual void get_smearing();
};

template <class ImplementationPolicy>
class StoutSmearingModule: public SmearingModule<ImplementationPolicy>{
   SmearedConfiguration<ImplementationPolicy> SmearingPolicy;
};

// Checkpoint module, owns the Checkpointer
template <class ImplementationPolicy>
class CheckPointModule {
  std::unique_ptr<BaseHmcCheckpointer<ImplementationPolicy> > cp_;

 public:
  void set_Checkpointer(BaseHmcCheckpointer<ImplementationPolicy>* cp) {
    cp_.reset(cp);
  };

  BaseHmcCheckpointer<ImplementationPolicy>* get_CheckPointer() {
    return cp_.get();
  }

  void initialize(CheckpointerParameters& P) { cp_.initialize(P); }
};

}  // namespace QCD
}  // namespace Grid

#endif  // GRID_HMC_MODULES