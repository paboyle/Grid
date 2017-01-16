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
  using ActionBase::ActionBase;

  // acquire resource
  virtual void initialize(){
    ActionBase::ActionPtr.reset(new WilsonGaugeAction<Impl>(ActionBase::Par_.beta));
  }

};

typedef WilsonGModule<PeriodicGimplR> WilsonGMod;


}// QCD temporarily here






// use the same classed defined by Antonin, does not make sense to rewrite
// Factory is perfectly fine
// Registar must be changed because I do not want to use the ModuleFactory
/*
define
*/


typedef HMCModuleBase< QCD::Action< QCD::LatticeGaugeField > > HMCModBase;

template <class ReaderClass >
class HMCActionModuleFactory
    : public Factory < HMCModBase ,	Reader<ReaderClass> > {
 public:
 	typedef Reader<ReaderClass> TheReader;
 	// use SINGLETON FUNCTOR MACRO HERE
  HMCActionModuleFactory(const HMCActionModuleFactory& e) = delete;
  void operator=(const HMCActionModuleFactory& e) = delete;
  static HMCActionModuleFactory& getInstance(void) {
    static HMCActionModuleFactory e;
    return e;
  }

 private:
  HMCActionModuleFactory(void) = default;
};

/*
then rewrite the registar

when this is done we have all the modules that contain the pointer to the objects
(actions, integrators, checkpointers, solvers)

factory will create only the modules and prepare the parameters
when needed a pointer is released


*/



template <class T, class TheFactory>
class Registrar {
 public:
  Registrar(std::string className) {
    // register the class factory function
    TheFactory::getInstance().registerBuilder(className, [&](typename TheFactory::TheReader Reader)
        { return std::unique_ptr<T>(new T(Reader));});
  }
};

Registrar<QCD::WilsonGMod, HMCActionModuleFactory<XmlReader> > __WGmodInit("WilsonGaugeAction"); 



}


#endif //HMC_MODULES_H