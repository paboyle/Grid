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
#ifndef HMC_RUNNER_MODULE
#define HMC_RUNNER_MODULE

NAMESPACE_BEGIN(Grid);

// the reader class is necessary here for the automatic initialization of the resources
// if we had a virtual reader would have been unecessary
template <class HMCType, class ReaderClass >
class HMCModule
  : public Parametrized< HMCparameters >,
    public HMCModuleBase< HMCRunnerBase<ReaderClass> > {
public:
  typedef HMCModuleBase< HMCRunnerBase<ReaderClass> > Base;
  typedef typename Base::Product Product;

  std::unique_ptr<HMCType> HMCPtr;

  HMCModule(HMCparameters Par) : Parametrized<HMCparameters>(Par) {}

  template <class ReaderCl>
  HMCModule(Reader<ReaderCl>& R) : Parametrized<HMCparameters>(R, "HMC"){};

  Product* getPtr() {
    if (!HMCPtr) initialize();
 
    return HMCPtr.get();
  }

private:
  virtual void initialize() = 0;
};

// Factory
template <char const *str, class ReaderClass >
class HMCRunnerModuleFactory
  : public Factory < HMCModuleBase< HMCRunnerBase<ReaderClass> > ,	Reader<ReaderClass> > {
public:
  typedef Reader<ReaderClass> TheReader; 
  // use SINGLETON FUNCTOR MACRO HERE
  HMCRunnerModuleFactory(const HMCRunnerModuleFactory& e) = delete;
  void operator=(const HMCRunnerModuleFactory& e) = delete;
  static HMCRunnerModuleFactory& getInstance(void) {
    static HMCRunnerModuleFactory e;
    return e;
  }

private:
  HMCRunnerModuleFactory(void) = default;
  std::string obj_type() const {
    return std::string(str);
  }
};





///////////////
// macro for these

template < class ImplementationPolicy, class RepresentationPolicy, class ReaderClass >
class HMCLeapFrog: public HMCModule< GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, LeapFrog>, ReaderClass >{
  typedef HMCModule< GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, LeapFrog>, ReaderClass  > HMCBaseMod;
  using HMCBaseMod::HMCBaseMod;

  // aquire resource
  virtual void initialize(){
    this->HMCPtr.reset(new GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, LeapFrog>(this->Par_) );
  }
};

template < class ImplementationPolicy, class RepresentationPolicy, class ReaderClass >
class HMCMinimumNorm2: public HMCModule< GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, MinimumNorm2>, ReaderClass  >{
  typedef HMCModule< GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, MinimumNorm2>, ReaderClass  > HMCBaseMod;
  using HMCBaseMod::HMCBaseMod;

  // aquire resource
  virtual void initialize(){
    this->HMCPtr.reset(new GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, MinimumNorm2>(this->Par_));
  }
};


template < class ImplementationPolicy, class RepresentationPolicy, class ReaderClass >
class HMCForceGradient: public HMCModule< GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, ForceGradient>, ReaderClass  >{
  typedef HMCModule< GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, ForceGradient>, ReaderClass   > HMCBaseMod;
  using HMCBaseMod::HMCBaseMod;

  // aquire resource
  virtual void initialize(){
    this->HMCPtr.reset(new GenericHMCRunnerTemplate<ImplementationPolicy, RepresentationPolicy, ForceGradient>(this->Par_) );
  }
};

extern char hmc_string[];

//////////////////////////////////////////////////////////////

NAMESPACE_END(Grid);

#endif
