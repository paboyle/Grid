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

// Indexing of tuple types
template <class T, class Tuple>
struct Index;

template <class T, class... Types>
struct Index<T, std::tuple<T, Types...>> {
  static const std::size_t value = 0;
};

template <class T, class U, class... Types>
struct Index<T, std::tuple<U, Types...>> {
  static const std::size_t value = 1 + Index<T, std::tuple<Types...>>::value;
};

/*
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
*/

template <class GaugeField, class Repr = NoHirep >
struct ActionLevel {
 public:
  unsigned int multiplier; 

  // Fundamental repr actions separated because of the smearing
  typedef Action<GaugeField>* ActPtr;

  // construct a tuple of vectors of the actions for the corresponding higher
  // representation fields
  typedef typename AccessTypes<Action, Repr>::VectorCollection action_collection;
  action_collection actions_hirep;
  typedef typename  AccessTypes<Action, Repr>::FieldTypeCollection action_hirep_types;

  std::vector<ActPtr>& actions;

  // Temporary conversion between ActionLevel and ActionLevelHirep
  //ActionLevelHirep(ActionLevel<GaugeField>& AL ):actions(AL.actions), multiplier(AL.multiplier){}

  ActionLevel(unsigned int mul = 1) : actions(std::get<0>(actions_hirep)), multiplier(mul) {
    // initialize the hirep vectors to zero.
    //apply(this->resize, actions_hirep, 0); //need a working resize
    assert(mul >= 1);
  };

  //void push_back(ActPtr ptr) { actions.push_back(ptr); }



  template < class Field >
  void push_back(Action<Field>* ptr) {
    // insert only in the correct vector
    std::get< Index < Field, action_hirep_types>::value >(actions_hirep).push_back(ptr);
  };

 

  template < class ActPtr>
  static void resize(ActPtr ap, unsigned int n){
    ap->resize(n);

  }

  //template <std::size_t I>
  //auto getRepresentation(Repr& R)->decltype(std::get<I>(R).U)  {return std::get<I>(R).U;}

  // Loop on tuple for a callable function
  template <std::size_t I = 1, typename Callable, typename ...Args>
  inline typename std::enable_if<I == std::tuple_size<action_collection>::value, void>::type apply(
      Callable, Repr& R,Args&...) const {}

  template <std::size_t I = 1, typename Callable, typename ...Args>
  inline typename std::enable_if<I < std::tuple_size<action_collection>::value, void>::type apply(
      Callable fn, Repr& R, Args&... arguments) const {
    fn(std::get<I>(actions_hirep), std::get<I>(R.rep), arguments...);
    apply<I + 1>(fn, R, arguments...);
  }  

};


//template <class GaugeField>
//using ActionSet = std::vector<ActionLevel<GaugeField> >;

template <class GaugeField, class R>
using ActionSet = std::vector<ActionLevel<GaugeField, R> >;

}
}
#endif
