/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

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
#ifndef HMC_MODULES_H
#define HMC_MODULES_H

/*
Define loadable, serializable modules
for the HMC execution
*/

namespace Grid {

// Empty class for no parameters
class NoParameters{};


/*
Base class for modules with parameters
*/
template < class P >
class Parametrized{
public:
  typedef P Parameters;

  Parametrized(Parameters Par):Par_(Par){};

  template <class ReaderClass>
  Parametrized(Reader<ReaderClass> & R, std::string section_name = "parameters"){
    read(R, section_name, Par_);
  }

  void set_parameters(Parameters Par){
        Par_ = Par;
  }

  void print_parameters(){
    std::cout << Par_ << std::endl;
  }

protected:
  Parameters Par_;
private:
  std::string section_name;
};


template <>
class Parametrized<NoParameters>{
        public:
  typedef NoParameters Parameters;

  Parametrized(Parameters Par){};

  template <class ReaderClass>
  Parametrized(Reader<ReaderClass> & Reader){};

  void set_parameters(Parameters Par){}

  void print_parameters(){}

};



/*
Lowest level abstract module class
*/
template <class Prod>
class HMCModuleBase {
 public:
  typedef Prod Product;

  virtual Prod* getPtr() = 0;

  // add a getReference? 
  
  virtual void print_parameters(){};  // default to nothing
};

/////////////////////////////
// Observables
/////////////////////////////
// explicit gauge field here....
template <class ObservableType, class OPar>
class ObservableModule
    : public Parametrized<OPar>,
      public HMCModuleBase< QCD::HmcObservable<typename ObservableType::GaugeField> > {
 public:
  typedef HMCModuleBase< QCD::HmcObservable< typename ObservableType::GaugeField> > Base;
  typedef typename Base::Product Product;

  std::unique_ptr<ObservableType> ObservablePtr;

  ObservableModule(OPar Par) : Parametrized<OPar>(Par) {}

  virtual void print_parameters(){
    std::cout << this->Par_ << std::endl;
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
class PlaquetteMod: public ObservableModule<PlaquetteLogger<Impl>, PlaquetteObsParameters>{
  typedef ObservableModule<PlaquetteLogger<Impl>, PlaquetteObsParameters> ObsBase;
  using ObsBase::ObsBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ObservablePtr.reset(new PlaquetteLogger<Impl>(this->Par_.output_prefix));
  }
};


}// QCD temporarily here







////////////////////////////////////////
// Factories specialisations
////////////////////////////////////////


// explicit ref to LatticeGaugeField must be changed of put in the factory
typedef HMCModuleBase< QCD::HmcObservable<QCD::LatticeGaugeField> > HMC_ObsModBase;

template <char const *str, class ReaderClass >
class HMC_ObservablesModuleFactory
    : public Factory < HMC_ObsModBase , Reader<ReaderClass> > {
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



template <class T, class TheFactory>
class Registrar {
 public:
  Registrar(std::string className) {
    // register the class factory function
    TheFactory::getInstance().registerBuilder(className, 
        [&](typename TheFactory::TheReader Reader)
        { 
          return std::unique_ptr<T>(new T(Reader));
        }
        );
  }
};



extern char observable_string[];
static Registrar<QCD::PlaquetteMod<QCD::PeriodicGimplR>, HMC_ObservablesModuleFactory<observable_string, XmlReader> > __OBSPLmodXMLInit("Plaquette"); 

}


#endif //HMC_MODULES_H