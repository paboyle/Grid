#ifndef GRID_MATH_TYPE_MAPPER_H
#define GRID_MATH_TYPE_MAPPER_H
namespace Grid {
//////////////////////////////////////////////////////////////////////////////////
// Want to recurse: GridTypeMapper<Matrix<vComplexD> >::scalar_type == ComplexD.
//////////////////////////////////////////////////////////////////////////////////
  
  template <class T> class GridTypeMapper {
  public:
    typedef typename T::scalar_type scalar_type;
    typedef typename T::vector_type vector_type;
    typedef typename T::tensor_reduced tensor_reduced;
  };

//////////////////////////////////////////////////////////////////////////////////
// Recursion stops with these template specialisations
//////////////////////////////////////////////////////////////////////////////////
  template<> class GridTypeMapper<RealF> {
  public:
    typedef RealF scalar_type;
    typedef RealF vector_type;
    typedef RealF tensor_reduced ;
  };
  template<> class GridTypeMapper<RealD> {
  public:
    typedef RealD scalar_type;
    typedef RealD vector_type;
    typedef RealD tensor_reduced;
  };
  template<> class GridTypeMapper<ComplexF> {
  public:
    typedef ComplexF scalar_type;
    typedef ComplexF vector_type;
    typedef ComplexF tensor_reduced;
  };
  template<> class GridTypeMapper<ComplexD> {
  public:
    typedef ComplexD scalar_type;
    typedef ComplexD vector_type;
    typedef ComplexD tensor_reduced;
  };

  template<> class GridTypeMapper<vRealF> {
  public:
    typedef RealF  scalar_type;
    typedef vRealF vector_type;
    typedef vRealF tensor_reduced;
  };
  template<> class GridTypeMapper<vRealD> {
  public:
    typedef RealD  scalar_type;
    typedef vRealD vector_type;
    typedef vRealD tensor_reduced;
  };
  template<> class GridTypeMapper<vComplexF> {
  public:
    typedef ComplexF  scalar_type;
    typedef vComplexF vector_type;
    typedef vComplexF tensor_reduced;
  };
  template<> class GridTypeMapper<vComplexD> {
  public:
    typedef ComplexD  scalar_type;
    typedef vComplexD vector_type;
    typedef vComplexD tensor_reduced;
  };
  template<> class GridTypeMapper<vInteger> {
  public:
    typedef Integer  scalar_type;
    typedef vInteger vector_type;
    typedef vInteger tensor_reduced;
  };

  // Again terminate the recursion.
  inline vRealD    TensorRemove(vRealD    arg){ return arg;}
  inline vRealF    TensorRemove(vRealF    arg){ return arg;}
  inline vComplexF TensorRemove(vComplexF arg){ return arg;}
  inline vComplexD TensorRemove(vComplexD arg){ return arg;}
  inline vInteger  TensorRemove(vInteger  arg){ return arg;}

}

#endif
