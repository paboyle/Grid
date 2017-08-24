/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/ObservableModules.h

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
#ifndef HMC_OBSERVABLE_MODULES_H
#define HMC_OBSERVABLE_MODULES_H

namespace Grid {

/////////////////////////////
// Observables
/////////////////////////////
template <class ObservableType, class OPar>
class ObservableModule
    : public Parametrized<OPar>,
      public HMCModuleBase< QCD::HmcObservable<typename ObservableType::Field> > {
 public:
  typedef HMCModuleBase< QCD::HmcObservable< typename ObservableType::Field> > Base;
  typedef typename Base::Product Product;
  typedef OPar Parameters;

  std::unique_ptr<ObservableType> ObservablePtr;

  ObservableModule(OPar Par) : Parametrized<OPar>(Par) {}

  virtual void print_parameters(){
    Parametrized<OPar>::print_parameters();
  }

  template <class ReaderClass>
  ObservableModule(Reader<ReaderClass>& Reader) : Parametrized<OPar>(Reader){};

  Product* getPtr() {
    if (!ObservablePtr) initialize();

    return ObservablePtr.get();
  }

 private:
  virtual void initialize() = 0;
};



////////////////
// Modules
////////////////

namespace QCD{

//// Observables module
class PlaquetteObsParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(PlaquetteObsParameters, 
    std::string, output_prefix);
};

template < class Impl >
class PlaquetteMod: public ObservableModule<PlaquetteLogger<Impl>, NoParameters>{
  typedef ObservableModule<PlaquetteLogger<Impl>, NoParameters> ObsBase;
  using ObsBase::ObsBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ObservablePtr.reset(new PlaquetteLogger<Impl>());
  }
  public:
  PlaquetteMod(): ObsBase(NoParameters()){}
};


template < class Impl >
class TopologicalChargeMod: public ObservableModule<TopologicalCharge<Impl>, TopologyObsParameters>{
  typedef ObservableModule<TopologicalCharge<Impl>, TopologyObsParameters> ObsBase;
  using ObsBase::ObsBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ObservablePtr.reset(new TopologicalCharge<Impl>(this->Par_));
  }
  public:
  TopologicalChargeMod(TopologyObsParameters Par): ObsBase(Par){}
  TopologicalChargeMod(): ObsBase(){}
};


}// QCD temporarily here


////////////////////////////////////////
// Factories specialisations
////////////////////////////////////////
// explicit ref to LatticeGaugeField must be changed or put in the factory
//typedef HMCModuleBase< QCD::HmcObservable<QCD::LatticeGaugeField> > HMC_ObsModBase;

template <char const *str, class Field, class ReaderClass >
class HMC_ObservablesModuleFactory
    : public Factory < HMCModuleBase< QCD::HmcObservable<Field> >, Reader<ReaderClass> > {
 public:
  typedef Reader<ReaderClass> TheReader; 
  // use SINGLETON FUNCTOR MACRO HERE
  HMC_ObservablesModuleFactory(const HMC_ObservablesModuleFactory& e) = delete;
  void operator=(const HMC_ObservablesModuleFactory& e) = delete;
  static HMC_ObservablesModuleFactory& getInstance(void) {
    static HMC_ObservablesModuleFactory e;
    return e;
  }

 private:
  HMC_ObservablesModuleFactory(void) = default;
    std::string obj_type() const {
    return std::string(str);
  }
};

extern char observable_string[];

}


#endif //HMC_OBSERVABLE_MODULES_H