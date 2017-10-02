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
  std::vector< QCD::GridModule* >    GridRefs;
 public:
  typedef HMCModuleBase< QCD::FermionOperator<FermionImpl> > Base;
  typedef typename Base::Product Product;

  FermionOperatorModule(FOPar Par) : Parametrized<FOPar>(Par) {}

  template <class ReaderClass>
  FermionOperatorModule(Reader<ReaderClass>& Reader) : Parametrized<FOPar>(Reader){};

  void AddGridPair(QCD::GridModule &Mod){
    if (GridRefs.size()>1){
      std::cout << GridLogError << "Adding too many Grids to the FermionOperatorModule" << std::endl;
      exit(1);
    }
    GridRefs.push_back(&Mod);

    if (Ls()){
      GridRefs.push_back(new QCD::GridModule());
      GridRefs[1]->set_full(QCD::SpaceTimeGrid::makeFiveDimGrid(Ls(),GridRefs[0]->get_full()));
      GridRefs[1]->set_rb(QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls(),GridRefs[0]->get_full()));
    }
  }

  virtual unsigned int Ls(){
    return 0;
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
    auto GridMod = this->GridRefs[0];
    typename FermionImpl::GaugeField U(GridMod->get_full());
    this->FOPtr.reset(new WilsonFermion<FermionImpl>(U, *(GridMod->get_full()), *(GridMod->get_rb()), this->Par_.mass));
  }
};



class MobiusFermionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MobiusFermionParameters,
    RealD, mass,
    RealD, M5,
    RealD, b,
    RealD, c,
    unsigned int, Ls);
};

template <class FermionImpl >
class MobiusFermionModule: public FermionOperatorModule<MobiusFermion, FermionImpl, MobiusFermionParameters> {
  typedef FermionOperatorModule<MobiusFermion, FermionImpl, MobiusFermionParameters> FermBase;
  using FermBase::FermBase; // for constructors

  virtual unsigned int Ls(){
    return this->Par_.Ls;
  }

  // acquire resource
  virtual void initialize(){
    auto GridMod = this->GridRefs[0];
    auto GridMod5d = this->GridRefs[1];
    typename FermionImpl::GaugeField U(GridMod->get_full());
    this->FOPtr.reset(new MobiusFermion<FermionImpl>( U, *(GridMod->get_full()), *(GridMod->get_rb()),
                                                      *(GridMod5d->get_full()), *(GridMod5d->get_rb()),
                                                      this->Par_.mass, this->Par_.M5, this->Par_.b, this->Par_.c));
  }
};


class DomainWallFermionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(DomainWallFermionParameters,
    RealD, mass,
    RealD, M5,
    unsigned int, Ls);
};

template <class FermionImpl >
class DomainWallFermionModule: public FermionOperatorModule<DomainWallFermion, FermionImpl, DomainWallFermionParameters> {
  typedef FermionOperatorModule<DomainWallFermion, FermionImpl, DomainWallFermionParameters> FermBase;
  using FermBase::FermBase; // for constructors

  virtual unsigned int Ls(){
    return this->Par_.Ls;
  }

  // acquire resource
  virtual void initialize(){
    auto GridMod = this->GridRefs[0];
    auto GridMod5d = this->GridRefs[1];
    typename FermionImpl::GaugeField U(GridMod->get_full());
    this->FOPtr.reset(new DomainWallFermion<FermionImpl>( U, *(GridMod->get_full()), *(GridMod->get_rb()),
                                                      *(GridMod5d->get_full()), *(GridMod5d->get_rb()),
                                                      this->Par_.mass, this->Par_.M5));
  }
};


class DomainWallEOFAFermionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(DomainWallEOFAFermionParameters,
    RealD, mq1,
    RealD, mq2,
    RealD, mq3,
    RealD, shift,
    int, pm,
    RealD, M5,
    unsigned int, Ls);
};

template <class FermionImpl >
class DomainWallEOFAFermionModule: public FermionOperatorModule<DomainWallEOFAFermion, FermionImpl, DomainWallEOFAFermionParameters> {
  typedef FermionOperatorModule<DomainWallEOFAFermion, FermionImpl, DomainWallEOFAFermionParameters> FermBase;
  using FermBase::FermBase; // for constructors

  virtual unsigned int Ls(){
    return this->Par_.Ls;
  }

  // acquire resource
  virtual void initialize(){
    auto GridMod = this->GridRefs[0];
    auto GridMod5d = this->GridRefs[1];
    typename FermionImpl::GaugeField U(GridMod->get_full());
    this->FOPtr.reset(new DomainWallEOFAFermion<FermionImpl>( U, *(GridMod->get_full()), *(GridMod->get_rb()),
                                                      *(GridMod5d->get_full()), *(GridMod5d->get_rb()),
                                                      this->Par_.mq1, this->Par_.mq2, this->Par_.mq3,
                                                      this->Par_.shift, this->Par_.pm, this->Par_.M5));
  }
};


} // QCD
} // Grid


#endif //FERMIONOPERATOR_MODULES_H
