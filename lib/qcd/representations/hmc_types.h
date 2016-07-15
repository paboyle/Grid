#ifndef HMC_TYPES_H
#define HMC_TYPES_H

#include <qcd/representations/adjoint.h>
#include <qcd/representations/fundamental.h>
#include <tuple>
#include <utility>

namespace Grid {
namespace QCD {

// Supported types
// enum {Fundamental, Adjoint} repr_type;

// Utility to add support to the HMC for representations other than the
// fundamental
template <class... Reptypes>
class Representations {
 public:
  typedef std::tuple<Reptypes...> Representation_type;

  // Size of the tuple, known at compile time
  static const int tuple_size = sizeof...(Reptypes);
  // The collection of types for the gauge fields
  typedef std::tuple<typename Reptypes::LatticeField...> Representation_Fields;

  // To access the Reptypes (FundamentalRepresentation, AdjointRepresentation)
  template <std::size_t N>
  using repr_type = typename std::tuple_element<N, Representation_type>::type;
  // in order to get the typename of the field use
  // type repr_type::LatticeField

  Representation_type rep;

  // Multiple types constructor
  explicit Representations(GridBase* grid) : rep(Reptypes(grid)...){};

  int size() { return tuple_size; }

  // update the fields
  template <std::size_t I = 0>
  inline typename std::enable_if<(I == tuple_size), void>::type update(
      LatticeGaugeField& U) {}

  template <std::size_t I = 0>
  inline typename std::enable_if<(I < tuple_size), void>::type update(
      LatticeGaugeField& U) {
    std::get<I>(rep).update_representation(U);
    update<I + 1>(U);
  }



};

typedef Representations<FundamentalRepresentation> NoHirep;

// Helper classes to access the elements
// Strips the first N parameters from the tuple
// sequence of classes to obtain the S sequence
// Creates a type that is a tuple of vectors of the template type A
template <template <typename> class A, class TupleClass,
          size_t N = TupleClass::tuple_size, size_t... S>
struct AccessTypes : AccessTypes<A, TupleClass, N - 1, N - 1, S...> {};

template <template <typename> class A, class TupleClass, size_t... S>
struct AccessTypes<A, TupleClass, 0, S...> {
 public:
  typedef typename TupleClass::Representation_Fields Rfields;

  template <std::size_t N>
  using elem = typename std::tuple_element<N, Rfields>::type;  // fields types

  typedef std::tuple<std::vector< A< elem<S> >* > ... > VectorCollection;
  typedef std::tuple< A< elem<S> >* ... > ClassCollection;

  // Debug
  void return_size() {
    std::cout << GridLogMessage
              << "Access:" << std::tuple_size<std::tuple<elem<S>...> >::value
              << "\n";
    std::cout << GridLogMessage
              << "Access vectors:" << std::tuple_size<VectorCollection>::value
              << "\n";
  }
};
}
}

#endif
