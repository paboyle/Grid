/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/FermionOperatorModules.h

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
#ifndef FERMIONOPERATOR_MODULES_H
#define FERMIONOPERATOR_MODULES_H

namespace Grid {

////////////////////////////////////
//  Fermion operators
/////////////////////////////////////
template < class Product>
class FermionOperatorModuleBase : public HMCModuleBase<Product>{
public:
  virtual void AddGridPair(QCD::GridModule&) = 0;
};

template <template <typename> class FOType, class FermionImpl, class FOPar>
class FermionOperatorModule
    : public Parametrized<FOPar>,
      public FermionOperatorModuleBase<QCD::FermionOperator<FermionImpl> > {

protected:
  std::unique_ptr< FOType<FermionImpl> > FOPtr;
  std::vector< std::reference_wrapper<QCD::GridModule> >    GridRefs;
 public:
  typedef HMCModuleBase< QCD::FermionOperator<FermionImpl> > Base;
  typedef typename Base::Product Product;

  FermionOperatorModule(FOPar Par) : Parametrized<FOPar>(Par) {}

  template <class ReaderClass>
  FermionOperatorModule(Reader<ReaderClass>& Reader) : Parametrized<FOPar>(Reader){};

  void AddGridPair(QCD::GridModule &Mod){
    if (GridRefs.size()>2){
      std::cout << GridLogError << "Adding too many Grids to the FermionOperatorModule" << std::endl;
      exit(1);
    }
    GridRefs.push_back(Mod);
  }

  virtual void print_parameters(){
    std::cout << this->Par_ << std::endl;
  }

  Product* getPtr() {
    if (!FOPtr) initialize();

    return FOPtr.get();
  }

 private:
  virtual void initialize() = 0;
};



// Factory
template <char const *str, class FermionImpl, class ReaderClass >
class HMC_FermionOperatorModuleFactory
    : public Factory < FermionOperatorModuleBase<QCD::FermionOperator<FermionImpl> > ,  Reader<ReaderClass> > {
 public:
  // use SINGLETON FUNCTOR MACRO HERE
  typedef Reader<ReaderClass> TheReader; 

  HMC_FermionOperatorModuleFactory(const HMC_FermionOperatorModuleFactory& e) = delete;
  void operator=(const HMC_FermionOperatorModuleFactory& e) = delete;
  static HMC_FermionOperatorModuleFactory& getInstance(void) {
    static HMC_FermionOperatorModuleFactory e;
    return e;
  }

 private:
  HMC_FermionOperatorModuleFactory(void) = default;
    std::string obj_type() const {
        return std::string(str);
  }
};




extern char fermionop_string[];
namespace QCD{

// Modules
class WilsonFermionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonFermionParameters, 
    RealD, mass);
};


template <class FermionImpl >
class WilsonFermionModule: public FermionOperatorModule<WilsonFermion, FermionImpl, WilsonFermionParameters> {
  typedef FermionOperatorModule<WilsonFermion, FermionImpl, WilsonFermionParameters> FermBase;
  using FermBase::FermBase; // for constructors

  // acquire resource
  virtual void initialize(){
    typename FermionImpl::GaugeField U(this->GridRefs[0].get().get_full());
    this->FOPtr.reset(new WilsonFermion<FermionImpl>(U, *(this->GridRefs[0].get().get_full()), *(this->GridRefs[0].get().get_rb()), this->Par_.mass));
  }
};


// Now a specific registration with a fermion field
static Registrar< WilsonFermionModule<WilsonImplR>,   
                  HMC_FermionOperatorModuleFactory<fermionop_string, WilsonImplR, XmlReader> > __WilsonFOPmodXMLInit("Wilson"); 


} // QCD




} // Grid


#endif //SOLVER_MODULES_H