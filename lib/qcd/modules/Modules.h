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

/*
Base class for modules with parameters
*/
template < class P >
class Parametrized{
public:
  typedef P Parameters;

  Parametrized(Parameters Par):Par_(Par){};

  template <class ReaderClass>
  Parametrized(Reader<ReaderClass> & Reader){
    read(Reader, section_name(), Par_);
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
  // identifies the section name
  // override in derived classes if needed 
  virtual std::string section_name(){
    return std::string("parameters"); //default
  }
};


/*
Lowest level abstract module class
*/
template < class Prod >
class HMCModuleBase{
public:
  typedef Prod Product;
virtual Prod* getPtr() = 0;  
virtual void print_parameters(){}; //default to nothing
};




//////////////////////////////////////////////
//		Actions
//////////////////////////////////////////////

template <class ActionType, class APar>
class ActionModule
    : public Parametrized<APar>,
      public HMCModuleBase<QCD::Action<typename ActionType::GaugeField> > {
 public:
  typedef HMCModuleBase< QCD::Action<typename ActionType::GaugeField> > Base;
  typedef typename Base::Product Product;



  std::unique_ptr<ActionType> ActionPtr;

  ActionModule(APar Par) : Parametrized<APar>(Par) {}

  template <class ReaderClass>
  ActionModule(Reader<ReaderClass>& Reader) : Parametrized<APar>(Reader){};

  virtual void print_parameters(){
  	std::cout << this->Par_ << std::endl;
  }

  Product* getPtr() {
    if (!ActionPtr) initialize();

    return ActionPtr.get();
  }

 private:
  virtual void initialize() = 0;

};








namespace QCD{

class WilsonGaugeActionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonGaugeActionParameters, 
    RealD, beta);
};




template<class Impl>
class WilsonGModule: public ActionModule<WilsonGaugeAction<Impl>, WilsonGaugeActionParameters> {
  typedef ActionModule<WilsonGaugeAction<Impl>, WilsonGaugeActionParameters> ActionBase;
  using ActionBase::ActionBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ActionPtr.reset(new WilsonGaugeAction<Impl>(this->Par_.beta));
  }

};

typedef WilsonGModule<PeriodicGimplR> WilsonGMod;


}// QCD temporarily here







////////////////////////////////////////
// Factories specialisations
////////////////////////////////////////



// use the same classed defined by Antonin, does not make sense to rewrite
// Factory is perfectly fine
// Registar must be changed because I do not want to use the ModuleFactory

// explicit ref to LatticeGaugeField must be changed
typedef HMCModuleBase< QCD::Action< QCD::LatticeGaugeField > > HMC_LGTActionModBase;

template <char const *str, class ReaderClass >
class HMC_LGTActionModuleFactory
    : public Factory < HMC_LGTActionModBase ,	Reader<ReaderClass> > {
 public:
 	typedef Reader<ReaderClass> TheReader; 
 	// use SINGLETON FUNCTOR MACRO HERE
  HMC_LGTActionModuleFactory(const HMC_LGTActionModuleFactory& e) = delete;
  void operator=(const HMC_LGTActionModuleFactory& e) = delete;
  static HMC_LGTActionModuleFactory& getInstance(void) {
    static HMC_LGTActionModuleFactory e;
    return e;
  }

 private:
  HMC_LGTActionModuleFactory(void) = default;
    std::string obj_type() const {
  	return std::string(str);
  }
};










template <class T, class TheFactory>
class Registrar {
 public:
  Registrar(std::string className) {
    // register the class factory function
    TheFactory::getInstance().registerBuilder(className, [&](typename TheFactory::TheReader Reader)
        { return std::unique_ptr<T>(new T(Reader));});
  }
};



extern char gauge_string[];
static Registrar<QCD::WilsonGMod, HMC_LGTActionModuleFactory<gauge_string, XmlReader> > __WGmodXMLInit("Wilson"); 
// add here the registration for other implementations and readers


}


#endif //HMC_MODULES_H