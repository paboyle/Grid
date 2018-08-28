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

#ifndef CP_MODULES_H
#define CP_MODULES_H


// FIXME  Reorganize QCD namespace 


namespace Grid {

////////////////////////////////////////////////////////////////////////
// Checkpoint module, owns the Checkpointer
////////////////////////////////////////////////////////////////////////

template <class ImplementationPolicy>
class CheckPointerModule: public Parametrized<QCD::CheckpointerParameters>, public HMCModuleBase< QCD::BaseHmcCheckpointer<ImplementationPolicy> >  {
 public:
 	std::unique_ptr<QCD::BaseHmcCheckpointer<ImplementationPolicy> > CheckPointPtr;
 	typedef QCD::CheckpointerParameters APar;
  typedef HMCModuleBase< QCD::BaseHmcCheckpointer<ImplementationPolicy> > Base;
  typedef typename Base::Product Product;

  CheckPointerModule(APar Par): Parametrized<APar>(Par) {}
  template <class ReaderClass>
  CheckPointerModule(Reader<ReaderClass>& Reader) : Parametrized<APar>(Reader){};

  virtual void print_parameters(){
  	std::cout << this->Par_ << std::endl;
  }

  Product* getPtr() {
    if (!CheckPointPtr) initialize();

    return CheckPointPtr.get();
  }

 private:
  virtual void initialize() = 0;

};



template <char const *str, class ImplementationPolicy, class ReaderClass >
class HMC_CPModuleFactory
    : public Factory < HMCModuleBase< QCD::BaseHmcCheckpointer<ImplementationPolicy> > ,	Reader<ReaderClass> > {
 public:
 	typedef Reader<ReaderClass> TheReader; 
 	// use SINGLETON FUNCTOR MACRO HERE
  HMC_CPModuleFactory(const HMC_CPModuleFactory& e) = delete;
  void operator=(const HMC_CPModuleFactory& e) = delete;
  static HMC_CPModuleFactory& getInstance(void) {
    static HMC_CPModuleFactory e;
    return e;
  }

 private:
  HMC_CPModuleFactory(void) = default;
  std::string obj_type() const {
  	return std::string(str);
  }
};



/////////////////////////////////////////////////////////////////////
// Concrete classes
/////////////////////////////////////////////////////////////////////
namespace QCD{

template<class ImplementationPolicy>
class BinaryCPModule: public CheckPointerModule< ImplementationPolicy> {
  typedef CheckPointerModule< ImplementationPolicy> CPBase;
  using CPBase::CPBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->CheckPointPtr.reset(new BinaryHmcCheckpointer<ImplementationPolicy>(this->Par_));
  }

};


template<class ImplementationPolicy>
class NerscCPModule: public CheckPointerModule< ImplementationPolicy> {
  typedef CheckPointerModule< ImplementationPolicy> CPBase;
  using CPBase::CPBase; // for constructors inheritance

  // acquire resource
  virtual void initialize(){
     this->CheckPointPtr.reset(new NerscHmcCheckpointer<ImplementationPolicy>(this->Par_));
  }

};


#ifdef HAVE_LIME
  
template<class ImplementationPolicy>
class ILDGCPModule: public CheckPointerModule< ImplementationPolicy> {
  typedef CheckPointerModule< ImplementationPolicy> CPBase;
  using CPBase::CPBase; // for constructors

  // acquire resource
  virtual void initialize(){
     this->CheckPointPtr.reset(new ILDGHmcCheckpointer<ImplementationPolicy>(this->Par_));
  }

};

template<class ImplementationPolicy, class Metadata>
class ScidacCPModule: public CheckPointerModule< ImplementationPolicy> {
  typedef CheckPointerModule< ImplementationPolicy> CPBase;
  Metadata M;

  // acquire resource
  virtual void initialize(){
     this->CheckPointPtr.reset(new ScidacHmcCheckpointer<ImplementationPolicy, Metadata>(this->Par_, M));
  }
public:
  ScidacCPModule(typename CPBase::APar Par, Metadata M_):M(M_), CPBase(Par) {}
  template <class ReaderClass>
  ScidacCPModule(Reader<ReaderClass>& Reader) : Parametrized<typename CPBase::APar>(Reader), M(Reader){};
};
#endif


}// QCD temporarily here


extern char cp_string[];

/*
// use macros?
static Registrar<QCD::BinaryCPModule<QCD::PeriodicGimplR>, HMC_CPModuleFactory<cp_string, QCD::PeriodicGimplR, XmlReader> > __CPBinarymodXMLInit("Binary");
static Registrar<QCD::NerscCPModule<QCD::PeriodicGimplR> , HMC_CPModuleFactory<cp_string, QCD::PeriodicGimplR, XmlReader> > __CPNerscmodXMLInit("Nersc");

#ifdef HAVE_LIME
static Registrar<QCD::ILDGCPModule<QCD::PeriodicGimplR>  , HMC_CPModuleFactory<cp_string, QCD::PeriodicGimplR, XmlReader> > __CPILDGmodXMLInit("ILDG");
#endif
*/

}// Grid
#endif //CP_MODULES_H
