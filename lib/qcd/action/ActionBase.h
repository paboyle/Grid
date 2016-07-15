/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/ActionBase.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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
#ifndef QCD_ACTION_BASE
#define QCD_ACTION_BASE
namespace Grid {
namespace QCD {

template <class GaugeField>
class Action {
 public:
  bool is_smeared = false;
  // Boundary conditions? // Heatbath?
  virtual void refresh(const GaugeField& U,
                       GridParallelRNG& pRNG) = 0;  // refresh pseudofermions
  virtual RealD S(const GaugeField& U) = 0;         // evaluate the action
  virtual void deriv(const GaugeField& U,
                     GaugeField& dSdU) = 0;  // evaluate the action derivative
  virtual ~Action(){};
};

// Could derive PseudoFermion action with a PF field, FermionField, and a Grid;
// implement refresh
/*
template<class GaugeField, class FermionField>
class PseudoFermionAction : public Action<GaugeField> {
 public:
  FermionField Phi;
  GridParallelRNG &pRNG;
  GridBase &Grid;

  PseudoFermionAction(GridBase &_Grid,GridParallelRNG &_pRNG) : Grid(_Grid),
Phi(&_Grid), pRNG(_pRNG) {
  };

  virtual void refresh(const GaugeField &gauge) {
    gaussian(Phi,pRNG);
  };

};
*/

template <class GaugeField>
struct ActionLevel {
 public:
  typedef Action<GaugeField>*
      ActPtr;  // now force the same colours as the rest of the code

  //Add supported representations here


  unsigned int multiplier;

  std::vector<ActPtr> actions;

  ActionLevel(unsigned int mul = 1) : actions(0), multiplier(mul) {
    assert(mul >= 1);
  };

  void push_back(ActPtr ptr) { actions.push_back(ptr); }
};


template <class GaugeField, class Repr>
struct ActionLevelHirep {
 public:
  unsigned int multiplier; 

  // Fundamental repr actions separated because of the smearing
  typedef Action<GaugeField>* ActPtr;
  //std::vector<ActPtr> actions;
  // construct a tuple of vectors of the actions for the corresponding higher
  // representation fields
  typename AccessTypes<Action, Repr>::VectorCollection actions_hirep;
  typedef typename  AccessTypes<Action, Repr>::ClassCollection actions_hirep_ptrs_type;

  std::vector<ActPtr>& actions;

  // Temporary conversion between ActionLevel and ActionLevelHirep
  ActionLevelHirep(ActionLevel<GaugeField>& AL ):actions(AL.actions), multiplier(AL.multiplier){}



  ActionLevelHirep(unsigned int mul = 1) : actions(std::get<0>(actions_hirep)), multiplier(mul) {
    // initialize the hirep vectors to zero.
    //apply(&ActionLevelHirep::resize, actions_hirep, 0); //need a working resize
    assert(mul >= 1);
  };

  void push_back(ActPtr ptr) { actions.push_back(ptr); }

// SFINAE construct, check
  template <class actionpointer, size_t N>
  void push_back(actionpointer ptr, decltype(std::tuple_element<N, actions_hirep_ptrs_type>::value)* = 0) {
    //insert only in the correct vector
    std::get<N>(actions_hirep).push_back(ptr);
  };

  template < class ActPtr>
  static void resize(ActPtr ap, unsigned int n){
    ap->resize(n);

  }

  

  // Loop on tuple for a callable function
  template <std::size_t I = 0, class Tuple, typename Callable, typename ...Args>
  inline typename std::enable_if<(I == std::tuple_size<Tuple>::value), void>::type apply(
      Callable&, Tuple& , Args...) {}

  template <std::size_t I = 0, class Tuple, typename Callable, typename ...Args>
  inline typename std::enable_if<(I < std::tuple_size<Tuple>::value), void>::type apply(
      Callable& fn,  Tuple& T, Args... arguments) {
    fn(std::get<I>(T), arguments...);
    apply<I + 1>(T, fn, arguments...);
  }  

};


template <class GaugeField>
using ActionSet = std::vector<ActionLevel<GaugeField> >;

template <class GaugeField, class R>
using ActionSetHirep = std::vector<ActionLevelHirep<GaugeField, R> >;

}
}
#endif
