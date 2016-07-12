#ifndef HMC_TYPES_H
#define HMC_TYPES_H

#include <tuple>
#include <utility>
#include <qcd/representations/adjoint.h>

namespace Grid {
namespace QCD {

// Utility to add support for representations other than the fundamental

template<class... Reptypes>
class Representations{
public:
  typedef std::tuple<Reptypes...> Representation_type;
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
      I<sizeof...(Reptypes), void >::type update(LatticeGaugeField& U) {
    std::get<I>(rep).update_representation(U);
    update<I + 1>(U);
  }  
};

}
}



#endif



