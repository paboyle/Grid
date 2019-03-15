/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/SolverModules.h

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
#ifndef SOLVER_MODULES_H
#define SOLVER_MODULES_H

namespace Grid {

//////////////////////////////////////////////
//       Operator Functions (Solvers)
//////////////////////////////////////////////

template <template <typename> class SolverType, class Field, class SPar>
class SolverModule
    : public Parametrized<SPar>,
      public HMCModuleBase<OperatorFunction<Field> > {
 public:
  typedef HMCModuleBase< OperatorFunction<Field> > Base;
  typedef typename Base::Product Product;

  std::unique_ptr< SolverType<Field> > SolverPtr;

  SolverModule(SPar Par) : Parametrized<SPar>(Par) {}

  template <class ReaderClass>
  SolverModule(Reader<ReaderClass>& Reader) : Parametrized<SPar>(Reader){};

  virtual void print_parameters(){
    std::cout << this->Par_ << std::endl;
  }

  Product* getPtr() {
    if (!SolverPtr) initialize();

    return SolverPtr.get();
  }

 private:
  virtual void initialize() = 0;
};


// Factory
template <char const *str, class Field, class ReaderClass >
class HMC_SolverModuleFactory
    : public Factory < HMCModuleBase<OperatorFunction<Field> > ,  Reader<ReaderClass> > {
 public:
  // use SINGLETON FUNCTOR MACRO HERE
  typedef Reader<ReaderClass> TheReader; 

  HMC_SolverModuleFactory(const HMC_SolverModuleFactory& e) = delete;
  void operator=(const HMC_SolverModuleFactory& e) = delete;
  static HMC_SolverModuleFactory& getInstance(void) {
    static HMC_SolverModuleFactory e;
    return e;
  }

 private:
  HMC_SolverModuleFactory(void) = default;
    std::string obj_type() const {
        return std::string(str);
  }
};



class SolverParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(SolverParameters,
    RealD, tolerance,
    RealD, max_iterations);
  // add error on no convergence?
};


class SolverObjName: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(SolverObjName, 
  std::string, name,
  SolverParameters, parameters);

};



template <class Field >
class ConjugateGradientModule: public SolverModule<ConjugateGradient, Field, SolverParameters> {
  typedef SolverModule<ConjugateGradient, Field, SolverParameters> SolverBase;
  using SolverBase::SolverBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->SolverPtr.reset(new ConjugateGradient<Field>(this->Par_.tolerance, this->Par_.max_iterations, true));
  }
};

template <class Field >
class ConjugateResidualModule: public SolverModule<ConjugateResidual, Field, SolverParameters> {
  typedef SolverModule<ConjugateResidual, Field, SolverParameters> SolverBase;
  using SolverBase::SolverBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->SolverPtr.reset(new ConjugateResidual<Field>(this->Par_.tolerance, this->Par_.max_iterations));
  }

};

extern char solver_string[];
} // Grid


#endif //SOLVER_MODULES_H