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
#define RegisterLoadCheckPointerFunction(NAME)                           \
  void Load##NAME##Checkpointer(const CheckpointerParameters& Params_) { \
    if (!have_CheckPointer) {                                            \
      std::cout << GridLogDebug << "Loading Checkpointer " << #NAME      \
                << std::endl;                                            \
      CP = std::unique_ptr<CheckpointerBaseModule>(                      \
        new NAME##CPModule<ImplementationPolicy>(Params_));              \
      have_CheckPointer = true;                                          \
    } else {                                                             \
      std::cout << GridLogError << "Checkpointer already loaded "        \
                << std::endl;                                            \
      exit(1);                                                           \
    }                                                                    \
  }

namespace Grid {
namespace QCD {

// HMC Resource manager
template <class ImplementationPolicy>
class HMCResourceManager {
  typedef HMCModuleBase< QCD::BaseHmcCheckpointer<ImplementationPolicy> > CheckpointerBaseModule;
  typedef HMCModuleBase< QCD::HmcObservable<typename ImplementationPolicy::Field> > ObservableBaseModule;
  typedef ActionModuleBase< QCD::Action<typename ImplementationPolicy::Field>, GridModule > ActionBaseModule;

  // Named storage for grid pairs (std + red-black)
  std::unordered_map<std::string, GridModule> Grids;
  RNGModule RNGs;

  // SmearingModule<ImplementationPolicy> Smearing;
  std::unique_ptr<CheckpointerBaseModule> CP;

  // A vector of HmcObservable modules
  std::vector<std::unique_ptr<ObservableBaseModule> > ObservablesList;


  // A vector of HmcObservable modules
  std::multimap<int, std::unique_ptr<ActionBaseModule> > ActionsList;
  std::vector<int> multipliers;

  bool have_RNG;
  bool have_CheckPointer;

  // NOTE: operator << is not overloaded for std::vector<string> 
  // so thsi function is necessary
  void output_vector_string(const std::vector<std::string> &vs){
    for (auto &i: vs)
      std::cout << i << " ";
    std::cout << std::endl;
  }


 public:
  HMCResourceManager() : have_RNG(false), have_CheckPointer(false) {}

  template <class ReaderClass, class vector_type = vComplex >
  void initialize(ReaderClass &Read){
    // assumes we are starting from the main node

    // Geometry
    GridModuleParameters GridPar(Read);
    GridFourDimModule<vector_type> GridMod( GridPar) ;
    AddGrid("gauge", GridMod);

    // Checkpointer
    auto &CPfactory = HMC_CPModuleFactory<cp_string, ImplementationPolicy, ReaderClass >::getInstance();
    Read.push("Checkpointer");
    std::string cp_type;
    read(Read,"name", cp_type);
    std::cout << "Registered types " << std::endl;
    output_vector_string(CPfactory.getBuilderList());


    CP = CPfactory.create(cp_type, Read);
    CP->print_parameters();
    Read.pop();    
    have_CheckPointer = true;  

    RNGModuleParameters RNGpar(Read);
    SetRNGSeeds(RNGpar);

    // Observables
    auto &ObsFactory = HMC_ObservablesModuleFactory<observable_string, typename ImplementationPolicy::Field, ReaderClass>::getInstance(); 
    Read.push(observable_string);// here must check if existing...
    do {
      std::string obs_type;
      read(Read,"name", obs_type);
      std::cout << "Registered types " << std::endl;
      output_vector_string(ObsFactory.getBuilderList() );

      ObservablesList.emplace_back(ObsFactory.create(obs_type, Read));
      ObservablesList[ObservablesList.size() - 1]->print_parameters();
    } while (Read.nextElement(observable_string));
    Read.pop();

    // Loop on levels
    if(!Read.push("Actions")){
      std::cout << "Actions not found" << std::endl; 
      exit(1);
    }

    if(!Read.push("Level")){// push must check if the node exist
         std::cout << "Level not found" << std::endl; 
      exit(1);
    }
    do
    {
      fill_ActionsLevel(Read); 
    }
    while(Read.push("Level"));

    Read.pop();
  }


 
  template <class RepresentationPolicy>
  void GetActionSet(ActionSet<typename ImplementationPolicy::Field, RepresentationPolicy>& Aset){
    Aset.resize(multipliers.size());
 
    for(auto it = ActionsList.begin(); it != ActionsList.end(); it++){
      (*it).second->acquireResource(Grids["gauge"]);
      Aset[(*it).first-1].push_back((*it).second->getPtr());
    }
  }



  //////////////////////////////////////////////////////////////
  // Grids
  //////////////////////////////////////////////////////////////

  void AddGrid(const std::string s, GridModule& M) {
    // Check for name clashes
    auto search = Grids.find(s);
    if (search != Grids.end()) {
      std::cout << GridLogError << "Grid with name \"" << search->first
                << "\" already present. Terminating\n";
      exit(1);
    }
    Grids[s] = std::move(M);
    std::cout << GridLogMessage << "::::::::::::::::::::::::::::::::::::::::" <<std::endl;
    std::cout << GridLogMessage << "HMCResourceManager:" << std::endl;
    std::cout << GridLogMessage << "Created grid set with name '" << s << "' and decomposition for the full cartesian " << std::endl;
    Grids[s].show_full_decomposition();
    std::cout << GridLogMessage << "::::::::::::::::::::::::::::::::::::::::" <<std::endl;
  }

  // Add a named grid set, 4d shortcut
  void AddFourDimGrid(const std::string s) {
    GridFourDimModule<vComplex> Mod;
    AddGrid(s, Mod);
  }

  // Add a named grid set, 4d shortcut + tweak simd lanes
  void AddFourDimGrid(const std::string s, const std::vector<int> simd_decomposition) {
    GridFourDimModule<vComplex> Mod(simd_decomposition);
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

  //////////////////////////////////////////////////////
  // Random number generators
  //////////////////////////////////////////////////////

  void AddRNGs(std::string s = "") {
    // Couple the RNGs to the GridModule tagged by s
    // the default is the first grid registered
    assert(Grids.size() > 0 && !have_RNG);
    if (s.empty()) s = Grids.begin()->first;
    std::cout << GridLogDebug << "Adding RNG to grid: " << s << std::endl;
    RNGs.set_pRNG(new GridParallelRNG(GetCartesian(s)));
    have_RNG = true;
  }

  void SetRNGSeeds(RNGModuleParameters& Params) { RNGs.set_RNGSeeds(Params); }

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

  BaseHmcCheckpointer<ImplementationPolicy>* GetCheckPointer() {
    if (have_CheckPointer)
      return CP->getPtr();
    else {
      std::cout << GridLogError << "Error: no checkpointer defined"
                << std::endl;
      exit(1);
    }
  }

  RegisterLoadCheckPointerFunction(Binary);
  RegisterLoadCheckPointerFunction(Nersc);
  #ifdef HAVE_LIME
  RegisterLoadCheckPointerFunction(ILDG);
  #endif

  ////////////////////////////////////////////////////////
  // Observables
  ////////////////////////////////////////////////////////

  template<class T, class... Types>
  void AddObservable(Types&&... Args){
    ObservablesList.push_back(std::unique_ptr<T>(new T(std::forward<Types>(Args)...)));
    ObservablesList.back()->print_parameters();
  }

  std::vector<HmcObservable<typename ImplementationPolicy::Field>* > GetObservables(){
    std::vector<HmcObservable<typename ImplementationPolicy::Field>* > out;
    for (auto &i : ObservablesList){
      out.push_back(i->getPtr());
    }

    // Add the checkpointer to the observables
    out.push_back(GetCheckPointer());
    return out;
  }



private:
   // this private
  template <class ReaderClass >
  void fill_ActionsLevel(ReaderClass &Read){
    // Actions set
    int m;
    Read.readDefault("multiplier",m);
    multipliers.push_back(m);
    std::cout << "Level : " << multipliers.size()  << " with multiplier : " << m << std::endl; 
    // here gauge
    Read.push("Action");
    do{
      auto &ActionFactory = HMC_ActionModuleFactory<gauge_string, typename ImplementationPolicy::Field, ReaderClass>::getInstance(); 
      std::string action_type;
      Read.readDefault("name", action_type); 
      output_vector_string(ActionFactory.getBuilderList() );
      ActionsList.emplace(m, ActionFactory.create(action_type, Read));
    } while (Read.nextElement("Action"));
    ActionsList.find(m)->second->print_parameters();    
    Read.pop();

  }



};
}
}

#endif  // HMC_RESOURCE_MANAGER_H
