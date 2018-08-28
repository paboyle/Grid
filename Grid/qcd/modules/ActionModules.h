/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/ActionModules.h

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
#ifndef ACTION_MODULES_H
#define ACTION_MODULES_H

/*
Define loadable, serializable modules
for the HMC execution
*/

namespace Grid {

//////////////////////////////////////////////
//              Actions
//////////////////////////////////////////////

template <class Product, class R>
class ActionModuleBase: public HMCModuleBase<Product>{
public:
  typedef R Resource;
  virtual void acquireResource(R& ){};

};


template <class ActionType, class APar>
class ActionModule
    : public Parametrized<APar>,
      public ActionModuleBase< QCD::Action<typename ActionType::GaugeField> , QCD::GridModule > {
 public:
  typedef ActionModuleBase< QCD::Action<typename ActionType::GaugeField>, QCD::GridModule > Base;
  typedef typename Base::Product Product;
  typedef APar Parameters;

  std::unique_ptr<ActionType> ActionPtr;

  ActionModule(APar Par) : Parametrized<APar>(Par) {}

  template <class ReaderClass>
  ActionModule(Reader<ReaderClass>& Reader) : Parametrized<APar>(Reader){};


  virtual void print_parameters(){
    Parametrized<APar>::print_parameters();
  }

  Product* getPtr() {
    if (!ActionPtr) initialize();

    return ActionPtr.get();
  }

 private:
  virtual void initialize() = 0;

};

//////////////////////////
// Modules
//////////////////////////

namespace QCD{

class PlaqPlusRectangleGaugeActionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(PlaqPlusRectangleGaugeActionParameters, 
    RealD, c_plaq,
    RealD, c_rect);

};

class RBCGaugeActionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(RBCGaugeActionParameters, 
    RealD, beta,
    RealD, c1);

};

class BetaGaugeActionParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(BetaGaugeActionParameters, 
    RealD, beta);
};




template <class Impl >
class WilsonGModule: public ActionModule<WilsonGaugeAction<Impl>, BetaGaugeActionParameters> {
  typedef ActionModule<WilsonGaugeAction<Impl>, BetaGaugeActionParameters> ActionBase;
  using ActionBase::ActionBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ActionPtr.reset(new WilsonGaugeAction<Impl>(this->Par_.beta));
  }

};

template <class Impl >
class PlaqPlusRectangleGModule: public ActionModule<PlaqPlusRectangleAction<Impl>, PlaqPlusRectangleGaugeActionParameters> {
  typedef ActionModule<PlaqPlusRectangleAction<Impl>, PlaqPlusRectangleGaugeActionParameters> ActionBase;
  using ActionBase::ActionBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ActionPtr.reset(new PlaqPlusRectangleAction<Impl>(this->Par_.c_plaq, this->Par_.c_rect));
  }

};

template <class Impl >
class RBCGModule: public ActionModule<RBCGaugeAction<Impl>, RBCGaugeActionParameters> {
  typedef ActionModule<RBCGaugeAction<Impl>, RBCGaugeActionParameters> ActionBase;
  using ActionBase::ActionBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ActionPtr.reset(new RBCGaugeAction<Impl>(this->Par_.beta, this->Par_.c1));
  }

};




template <class Impl >
class SymanzikGModule: public ActionModule<SymanzikGaugeAction<Impl>, BetaGaugeActionParameters> {
  typedef ActionModule<SymanzikGaugeAction<Impl>, BetaGaugeActionParameters> ActionBase;
  using ActionBase::ActionBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ActionPtr.reset(new SymanzikGaugeAction<Impl>(this->Par_.beta));
  }

};

template <class Impl >
class IwasakiGModule: public ActionModule<IwasakiGaugeAction<Impl>, BetaGaugeActionParameters> {
  typedef ActionModule<IwasakiGaugeAction<Impl>, BetaGaugeActionParameters> ActionBase;
  using ActionBase::ActionBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ActionPtr.reset(new IwasakiGaugeAction<Impl>(this->Par_.beta));
  }

};


template <class Impl >
class DBW2GModule: public ActionModule<DBW2GaugeAction<Impl>, BetaGaugeActionParameters> {
  typedef ActionModule<DBW2GaugeAction<Impl>, BetaGaugeActionParameters> ActionBase;
  using ActionBase::ActionBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ActionPtr.reset(new DBW2GaugeAction<Impl>(this->Par_.beta));
  }

};

/////////////////////////////////////////
// Fermion Actions
/////////////////////////////////////////


template <class Impl, template <typename> class FermionA, class Params = NoParameters >
class PseudoFermionModuleBase: public ActionModule<FermionA<Impl>, Params> {
protected:
  typedef ActionModule<FermionA<Impl>, Params> ActionBase;
  using ActionBase::ActionBase; // for constructors

  typedef std::unique_ptr<FermionOperatorModuleBase<FermionOperator<Impl>> > operator_type;
  typedef std::unique_ptr<HMCModuleBase<OperatorFunction<typename Impl::FermionField> > > solver_type;

  template <class ReaderClass>
  void getFermionOperator(Reader<ReaderClass>& Reader, operator_type &fo, std::string section_name){
    auto &FOFactory = HMC_FermionOperatorModuleFactory<fermionop_string, Impl, ReaderClass>::getInstance();
    Reader.push(section_name);
    std::string op_name;
    read(Reader,"name", op_name);
    fo = FOFactory.create(op_name, Reader);
    Reader.pop();  
  }

  template <class ReaderClass>
  void getSolverOperator(Reader<ReaderClass>& Reader, solver_type &so, std::string section_name){
    auto& SolverFactory = HMC_SolverModuleFactory<solver_string, typename Impl::FermionField, ReaderClass>::getInstance();
    Reader.push(section_name);
    std::string solv_name;
    read(Reader,"name", solv_name);
    so = SolverFactory.create(solv_name, Reader);
    Reader.pop();    
  }
};


template <class Impl >
class TwoFlavourFModule: public PseudoFermionModuleBase<Impl, TwoFlavourPseudoFermionAction>{
  typedef PseudoFermionModuleBase<Impl, TwoFlavourPseudoFermionAction> Base;
  using Base::Base;

  typename Base::operator_type fop_mod;
  typename Base::solver_type   solver_mod;

 public:
  virtual void acquireResource(typename Base::Resource& GridMod){
    fop_mod->AddGridPair(GridMod);
  }

   // constructor
   template <class ReaderClass>
   TwoFlavourFModule(Reader<ReaderClass>& R): Base(R) {
    this->getSolverOperator(R, solver_mod, "Solver");
    this->getFermionOperator(R, fop_mod, "Operator");
   } 

  // acquire resource
  virtual void initialize() {
    // here temporarily assuming that the force and action solver are the same
    this->ActionPtr.reset(new TwoFlavourPseudoFermionAction<Impl>(*(this->fop_mod->getPtr()), *(this->solver_mod->getPtr()), *(this->solver_mod->getPtr())));
  }

};

// very similar, I could have templated this but it is overkilling
template <class Impl >
class TwoFlavourEOFModule: public PseudoFermionModuleBase<Impl, TwoFlavourEvenOddPseudoFermionAction>{
  typedef PseudoFermionModuleBase<Impl, TwoFlavourEvenOddPseudoFermionAction> Base;
  using Base::Base;

  typename Base::operator_type fop_mod;
  typename Base::solver_type   solver_mod;

 public:
  virtual void acquireResource(typename Base::Resource& GridMod){
    fop_mod->AddGridPair(GridMod);
  }

   // constructor
   template <class ReaderClass>
   TwoFlavourEOFModule(Reader<ReaderClass>& R): PseudoFermionModuleBase<Impl, TwoFlavourEvenOddPseudoFermionAction>(R) {
    this->getSolverOperator(R, solver_mod, "Solver");
    this->getFermionOperator(R, fop_mod, "Operator");
   } 

  // acquire resource
  virtual void initialize() {
    // here temporarily assuming that the force and action solver are the same
    this->ActionPtr.reset(new TwoFlavourEvenOddPseudoFermionAction<Impl>(*(this->fop_mod->getPtr()), *(this->solver_mod->getPtr()), *(this->solver_mod->getPtr())));
  }

};


template <class Impl >
class TwoFlavourRatioFModule: public PseudoFermionModuleBase<Impl, TwoFlavourRatioPseudoFermionAction>{
  typedef PseudoFermionModuleBase<Impl, TwoFlavourRatioPseudoFermionAction> Base;
  using Base::Base;

  typename Base::operator_type fop_numerator_mod;
  typename Base::operator_type fop_denominator_mod;
  typename Base::solver_type   solver_mod;

 public:
  virtual void acquireResource(typename Base::Resource& GridMod){
    fop_numerator_mod->AddGridPair(GridMod);
    fop_denominator_mod->AddGridPair(GridMod);
  }

   // constructor
   template <class ReaderClass>
   TwoFlavourRatioFModule(Reader<ReaderClass>& R): PseudoFermionModuleBase<Impl, TwoFlavourRatioPseudoFermionAction>(R) {
    this->getSolverOperator(R, solver_mod, "Solver");
    this->getFermionOperator(R, fop_numerator_mod, "Numerator");
    this->getFermionOperator(R, fop_denominator_mod, "Denominator");
   } 

  // acquire resource
  virtual void initialize() {
    // here temporarily assuming that the force and action solver are the same
    this->ActionPtr.reset(new TwoFlavourRatioPseudoFermionAction<Impl>(*(this->fop_numerator_mod->getPtr()), 
      *(this->fop_denominator_mod->getPtr()), *(this->solver_mod->getPtr()), *(this->solver_mod->getPtr())));
  }

};

template <class Impl >
class TwoFlavourRatioEOFModule: public PseudoFermionModuleBase<Impl, TwoFlavourEvenOddRatioPseudoFermionAction>{
  typedef PseudoFermionModuleBase<Impl, TwoFlavourEvenOddRatioPseudoFermionAction> Base;
  using Base::Base;

  typename Base::operator_type fop_numerator_mod;
  typename Base::operator_type fop_denominator_mod;
  typename Base::solver_type   solver_mod;

 public:
  virtual void acquireResource(typename Base::Resource& GridMod){
    fop_numerator_mod->AddGridPair(GridMod);
    fop_denominator_mod->AddGridPair(GridMod);
  }

   // constructor
   template <class ReaderClass>
   TwoFlavourRatioEOFModule(Reader<ReaderClass>& R): Base(R) {
    this->getSolverOperator(R, solver_mod, "Solver");
    this->getFermionOperator(R, fop_numerator_mod, "Numerator");
    this->getFermionOperator(R, fop_denominator_mod, "Denominator");
   } 

  // acquire resource
  virtual void initialize() {
    // here temporarily assuming that the force and action solver are the same
    this->ActionPtr.reset(new TwoFlavourEvenOddRatioPseudoFermionAction<Impl>(*(this->fop_numerator_mod->getPtr()), 
      *(this->fop_denominator_mod->getPtr()), *(this->solver_mod->getPtr()), *(this->solver_mod->getPtr())));
  }

};


template <class Impl >
class OneFlavourFModule: public PseudoFermionModuleBase<Impl, OneFlavourRationalPseudoFermionAction, OneFlavourRationalParams>{
  typedef PseudoFermionModuleBase<Impl, OneFlavourRationalPseudoFermionAction, OneFlavourRationalParams> Base;
  using Base::Base;

  typename Base::operator_type fop_mod;

 public:
  virtual void acquireResource(typename Base::Resource& GridMod){
    fop_mod->AddGridPair(GridMod);
  }

   // constructor
   template <class ReaderClass>
   OneFlavourFModule(Reader<ReaderClass>& R): Base(R) {
    this->getFermionOperator(R, fop_mod, "Operator");
   } 

  // acquire resource
  virtual void initialize() {
    this->ActionPtr.reset(new OneFlavourRationalPseudoFermionAction<Impl>(*(this->fop_mod->getPtr()), this->Par_ ));
  }

};

template <class Impl >
class OneFlavourEOFModule: 
  public PseudoFermionModuleBase<Impl, OneFlavourEvenOddRationalPseudoFermionAction, OneFlavourRationalParams>
  {
  typedef PseudoFermionModuleBase<Impl, OneFlavourEvenOddRationalPseudoFermionAction, OneFlavourRationalParams> Base;
  using Base::Base;

  typename Base::operator_type fop_mod;

 public:
  virtual void acquireResource(typename Base::Resource& GridMod){
    fop_mod->AddGridPair(GridMod);
  }

   // constructor
   template <class ReaderClass>
   OneFlavourEOFModule(Reader<ReaderClass>& R): Base(R) {
    this->getFermionOperator(R, fop_mod, "Operator");
   } 

  // acquire resource
  virtual void initialize() {
    this->ActionPtr.reset(new OneFlavourEvenOddRationalPseudoFermionAction<Impl>(*(this->fop_mod->getPtr()), this->Par_ ));
  }

};


template <class Impl >
class OneFlavourRatioFModule: 
  public PseudoFermionModuleBase<Impl, OneFlavourRatioRationalPseudoFermionAction, OneFlavourRationalParams>
  {

  typedef PseudoFermionModuleBase<Impl, OneFlavourRatioRationalPseudoFermionAction, OneFlavourRationalParams> Base;
  using Base::Base;

  typename Base::operator_type fop_numerator_mod;
  typename Base::operator_type fop_denominator_mod;

 public:
  virtual void acquireResource(typename Base::Resource& GridMod){
    fop_numerator_mod->AddGridPair(GridMod);
    fop_denominator_mod->AddGridPair(GridMod);
  }

   // constructor
   template <class ReaderClass>
   OneFlavourRatioFModule(Reader<ReaderClass>& R): Base(R) {
    this->getFermionOperator(R, fop_numerator_mod, "Numerator");
    this->getFermionOperator(R, fop_denominator_mod, "Denominator");
   } 

  // acquire resource
  virtual void initialize() {
    this->ActionPtr.reset(new OneFlavourRatioRationalPseudoFermionAction<Impl>( *(this->fop_numerator_mod->getPtr()), 
                                                                                *(this->fop_denominator_mod->getPtr()), 
                                                                                this->Par_ ));
  }

};


template <class Impl >
class OneFlavourRatioEOFModule: 
  public PseudoFermionModuleBase<Impl, OneFlavourEvenOddRatioRationalPseudoFermionAction, OneFlavourRationalParams>
  {

  typedef PseudoFermionModuleBase<Impl, OneFlavourEvenOddRatioRationalPseudoFermionAction, OneFlavourRationalParams> Base;
  using Base::Base;

  typename Base::operator_type fop_numerator_mod;
  typename Base::operator_type fop_denominator_mod;

 public:
  virtual void acquireResource(typename Base::Resource& GridMod){
    fop_numerator_mod->AddGridPair(GridMod);
    fop_denominator_mod->AddGridPair(GridMod);
  }

   // constructor
   template <class ReaderClass>
   OneFlavourRatioEOFModule(Reader<ReaderClass>& R): Base(R) {
    this->getFermionOperator(R, fop_numerator_mod, "Numerator");
    this->getFermionOperator(R, fop_denominator_mod, "Denominator");
   } 

  // acquire resource
  virtual void initialize() {
    this->ActionPtr.reset(new OneFlavourEvenOddRatioRationalPseudoFermionAction<Impl>(*(this->fop_numerator_mod->getPtr()), 
                                                                                      *(this->fop_denominator_mod->getPtr()), 
                                                                                      this->Par_ ));
  }

};

}// QCD temporarily here







////////////////////////////////////////
// Factories specialisations
////////////////////////////////////////



// use the same classed defined by Antonin, does not make sense to rewrite
// Factory is perfectly fine
// Registar must be changed because I do not want to use the ModuleFactory

// explicit ref to LatticeGaugeField must be changed or put in the factory
//typedef ActionModuleBase< QCD::Action< QCD::LatticeGaugeField >, QCD::GridModule > HMC_LGTActionModBase;
//typedef ActionModuleBase< QCD::Action< QCD::LatticeReal >, QCD::GridModule > HMC_ScalarActionModBase;

template <char const *str, class Field, class ReaderClass >
class HMC_ActionModuleFactory
    : public Factory < ActionModuleBase< QCD::Action< Field >, QCD::GridModule > , Reader<ReaderClass> > {
 public:
  typedef Reader<ReaderClass> TheReader; 
  // use SINGLETON FUNCTOR MACRO HERE
  HMC_ActionModuleFactory(const HMC_ActionModuleFactory& e) = delete;
  void operator=(const HMC_ActionModuleFactory& e) = delete;
  static HMC_ActionModuleFactory& getInstance(void) {
    static HMC_ActionModuleFactory e;
    return e;
  }

 private:
  HMC_ActionModuleFactory(void) = default;
    std::string obj_type() const {
        return std::string(str);
  }
};


extern char gauge_string[];
} // Grid


#endif //HMC_MODULES_H