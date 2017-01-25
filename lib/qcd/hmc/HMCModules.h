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

// call these: resources

// Some modules for the basic setup


///////////////////////////////////////////////////
// Modules
class GridModuleParameters: Serializable{
  
  GRID_SERIALIZABLE_CLASS_MEMBERS(GridModuleParameters,
  std::string, lattice,
  std::string,  mpi);

public: 
  // these namings are ugly
  // also ugly the distinction between the serializable members
  // and this
  std::vector<int> lattice_v;
  std::vector<int> mpi_v;

  GridModuleParameters(const std::vector<int> l_ = std::vector<int>(),
    const std::vector<int> mpi_ = std::vector<int>()):lattice_v(l_), mpi_v(mpi_){}

  template <class ReaderClass>
  GridModuleParameters(Reader<ReaderClass>& Reader) {
    read(Reader, "LatticeGrid", *this);
    lattice_v = strToVec<int>(lattice);
    mpi_v = strToVec<int>(mpi);
    if (mpi_v.size() != lattice_v.size()) {
      std::cout << "Error in GridModuleParameters: lattice and mpi dimensions "
                   "do not match"
                << std::endl;
      exit(1);
    }
  }
};

class GridModule {
 public:
  GridCartesian* get_full() { 
    std::cout << GridLogDebug << "Getting cartesian in module"<< std::endl;
    return grid_.get(); }
  GridRedBlackCartesian* get_rb() { 
    std::cout << GridLogDebug << "Getting rb-cartesian in module"<< std::endl;
    return rbgrid_.get(); }

  void set_full(GridCartesian* grid) { grid_.reset(grid); }
  void set_rb(GridRedBlackCartesian* rbgrid) { rbgrid_.reset(rbgrid); }

 protected:
  std::unique_ptr<GridCartesian> grid_;
  std::unique_ptr<GridRedBlackCartesian> rbgrid_;
  
};

// helpers
class GridFourDimModule : public GridModule {
 public:
  // add a function to create the module from a Reader
  GridFourDimModule() {
    set_full(SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(4, vComplex::Nsimd()),
        GridDefaultMpi()));
    set_rb(SpaceTimeGrid::makeFourDimRedBlackGrid(grid_.get()));
  }

  template <class vector_type = vComplex>
  GridFourDimModule(GridModuleParameters Params) {
    if (Params.lattice_v.size() == 4) {
      set_full(SpaceTimeGrid::makeFourDimGrid(
          Params.lattice_v, GridDefaultSimd(4, vector_type::Nsimd()),
          Params.mpi_v));
      set_rb(SpaceTimeGrid::makeFourDimRedBlackGrid(grid_.get()));
    } else {
      std::cout
          << "Error in GridFourDimModule: lattice dimension different from 4"
          << std::endl;
      exit(1);
    }
  }
};


////////////////////////////////////////////////////////////////////
class RNGModuleParameters: Serializable {

  GRID_SERIALIZABLE_CLASS_MEMBERS(RNGModuleParameters,
  std::string, serial_seeds,
  std::string, parallel_seeds,);
public:
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