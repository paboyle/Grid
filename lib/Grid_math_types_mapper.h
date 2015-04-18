#ifndef GRID_MATH_TYPE_MAPPER_H
#define GRID_MATH_TYPE_MAPPER_H

namespace Grid {


//////////////////////////////////////////////////////////////////////////////////
// Want to recurse: GridTypeMapper<Matrix<vComplexD> >::scalar_type == ComplexD.
// Use of a helper class like this allows us to template specialise and "dress"
// other classes such as RealD == double, ComplexD == std::complex<double> with these
// traits.
//
// It is possible that we could do this more elegantly if I introduced a 
// queryable trait in iScalar, iMatrix and iVector and used the query on vtype in 
// place of the type mapper?
//
// Not sure how to do this, but probably could be done with a research effort
// to study C++11's type_traits.h file. (std::enable_if<isGridTensorType<vtype> >)
//
//////////////////////////////////////////////////////////////////////////////////
  
  template <class T> class GridTypeMapper {
  public:
    typedef typename T::scalar_type scalar_type;
    typedef typename T::vector_type vector_type;
    typedef typename T::tensor_reduced tensor_reduced;
    enum { TensorLevel = T::TensorLevel };
  };

//////////////////////////////////////////////////////////////////////////////////
// Recursion stops with these template specialisations
//////////////////////////////////////////////////////////////////////////////////
  template<> class GridTypeMapper<RealF> {
  public:
    typedef RealF scalar_type;
    typedef RealF vector_type;
    typedef RealF tensor_reduced ;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<RealD> {
  public:
    typedef RealD scalar_type;
    typedef RealD vector_type;
    typedef RealD tensor_reduced;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<ComplexF> {
  public:
    typedef ComplexF scalar_type;
    typedef ComplexF vector_type;
    typedef ComplexF tensor_reduced;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<ComplexD> {
  public:
    typedef ComplexD scalar_type;
    typedef ComplexD vector_type;
    typedef ComplexD tensor_reduced;
    enum { TensorLevel = 0 };
  };

  template<> class GridTypeMapper<vRealF> {
  public:
    typedef RealF  scalar_type;
    typedef vRealF vector_type;
    typedef vRealF tensor_reduced;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vRealD> {
  public:
    typedef RealD  scalar_type;
    typedef vRealD vector_type;
    typedef vRealD tensor_reduced;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vComplexF> {
  public:
    typedef ComplexF  scalar_type;
    typedef vComplexF vector_type;
    typedef vComplexF tensor_reduced;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vComplexD> {
  public:
    typedef ComplexD  scalar_type;
    typedef vComplexD vector_type;
    typedef vComplexD tensor_reduced;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vInteger> {
  public:
    typedef Integer  scalar_type;
    typedef vInteger vector_type;
    typedef vInteger tensor_reduced;
    enum { TensorLevel = 0 };
  };

}

#endif
