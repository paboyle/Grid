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
#ifndef HMC_RESOURCE_MANAGER_H
#define HMC_RESOURCE_MANAGER_H

#include <unordered_map>

// One function per Checkpointer, use a macro to simplify
  #define RegisterLoadCheckPointerFunction(NAME)                                       \
  void Load##NAME##Checkpointer(CheckpointerParameters& Params_) {   \
    if (!have_CheckPointer) {                                        \
      std::cout << GridLogDebug << "Loading Checkpointer " << #NAME  \
                << std::endl;                                        \
      CP.set_Checkpointer(                                           \
          new NAME##HmcCheckpointer<ImplementationPolicy>(Params_)); \
      have_CheckPointer = true;                                      \
    } else {                                                         \
      std::cout << GridLogError << "Checkpointer already loaded "    \
                << std::endl;                                        \
      exit(1);                                                       \
    }                                                                \
  }



namespace Grid {
namespace QCD {

// HMC Resource manager
  template <class ImplementationPolicy>
class HMCResourceManager{
  // Storage for grid pairs (std + red-black)
  std::unordered_map<std::string, GridModule> Grids;
  RNGModule RNGs;

  //SmearingModule<ImplementationPolicy> Smearing;
  CheckPointModule<ImplementationPolicy> CP;

  bool have_RNG;
  bool have_CheckPointer;

 public:
  HMCResourceManager() : have_RNG(false), have_CheckPointer(false) {}
  void AddGrid(std::string s, GridModule& M) {
    // Check for name clashes
    auto search = Grids.find(s);
    if (search != Grids.end()) {
      std::cout << GridLogError << "Grid with name \"" << search->first
                << "\" already present. Terminating\n";
      exit(1);
    }
    Grids[s] = std::move(M);
  }

  // Add a named grid set
  void AddFourDimGrid(std::string s) {
    GridFourDimModule Mod;
    AddGrid(s, Mod);
  }

  GridCartesian* GetCartesian(std::string s = "") {
    if (s.empty()) s = Grids.begin()->first;
    std::cout << GridLogDebug << "Getting cartesian grid from: " << s
              << std::endl;
    return Grids[s].get_full();
  }

  GridRedBlackCartesian* GetRBCartesian(std::string s = "") {
    if (s.empty()) s = Grids.begin()->first;
    std::cout << GridLogDebug << "Getting rb-cartesian grid from: " << s
              << std::endl;
    return Grids[s].get_rb();
  }

  void AddRNGs(std::string s = "") {
    // Couple the RNGs to the GridModule tagged by s
    // the default is the first grid registered
    assert(Grids.size() > 0 && !have_RNG);
    if (s.empty()) s = Grids.begin()->first;
    std::cout << GridLogDebug << "Adding RNG to grid: " << s << std::endl;
    RNGs.set_pRNG(new GridParallelRNG(GetCartesian(s)));
    have_RNG = true;
  }

  void AddRNGSeeds(const std::vector<int> S, const std::vector<int> P) {
    RNGs.set_RNGSeeds(S, P);
  }

  GridSerialRNG& GetSerialRNG() { return RNGs.get_sRNG(); }
  GridParallelRNG& GetParallelRNG() {
    assert(have_RNG);
    return RNGs.get_pRNG();
  }
  
  void SeedFixedIntegers() {
    assert(have_RNG);
    RNGs.seed();
  }



  //////////////////////////////////////////////////////
  // Checkpointers
  //////////////////////////////////////////////////////

  BaseHmcCheckpointer<ImplementationPolicy>* get_CheckPointer(){
    if (have_CheckPointer)
    return CP.get_CheckPointer();
    else{
      std::cout << GridLogError << "Error: no checkpointer defined" << std::endl;
      exit(1);
    }
  }

  RegisterLoadCheckPointerFunction (Binary);
  RegisterLoadCheckPointerFunction (Nersc);
  RegisterLoadCheckPointerFunction (ILDG)

};
}
}

#endif  // HMC_RESOURCE_MANAGER_H