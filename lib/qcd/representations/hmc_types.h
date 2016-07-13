#ifndef HMC_TYPES_H
#define HMC_TYPES_H

#include <tuple>
#include <utility>
#include <qcd/representations/adjoint.h>
#include <qcd/representations/fundamental.h>

namespace Grid {
namespace QCD {


// Supported types
//enum {Fundamental, Adjoint} repr_type;

// Utility to add support to the HMC for representations other than the fundamental
template<class... Reptypes>
class Representations{
public:
  typedef std::tuple<Reptypes...> Representation_type;

  // To access the Reptypes (FundamentalRepresentation, AdjointRepresentation)
  template <std::size_t N>
  using repr_type = typename std::tuple_element<N, Representation_type >::type;
  // in order to get the typename of the field use
  // type repr_type::LatticeField

  Representation_type rep;

  // Multiple types constructor
  explicit Representations(GridBase *grid):rep(Reptypes(grid)...){};

  int size(){
    return std::tuple_size< Representation_type >::value;
  }

  // update the fields
  template <std::size_t I = 0>
  inline typename std::enable_if< I == sizeof...(Reptypes), void >::type update(LatticeGaugeField& U) {}

  template <std::size_t I = 0>
      inline typename std::enable_if <
      I<sizeof...(Reptypes), void>::type update(LatticeGaugeField& U) {
    std::get<I>(rep).update_representation(U);
    update<I + 1>(U);
  }  
};


typedef Representations<FundamentalRepresentation> JustTheFundamental;


}
}



#endif



