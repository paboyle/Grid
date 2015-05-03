#ifndef GRID_MATH_TENSORS_H
#define GRID_MATH_TENSORS_H

namespace Grid {

///////////////////////////////////////////////////
// Scalar, Vector, Matrix objects.
// These can be composed to form tensor products of internal indices.
///////////////////////////////////////////////////

// It is useful to NOT have any constructors
// so that these classes assert "is_pod<class> == true"
// because then the standard C++ valarray container eliminates fill overhead on new allocation and 
// non-move copying.
//
// However note that doing this eliminates some syntactical sugar such as 
// calling the constructor explicitly or implicitly
//
#undef TENSOR_IS_POD

template<class vtype> class iScalar
{
public:
  vtype _internal;

  typedef typename GridTypeMapper<vtype>::scalar_type   scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;
  typedef iScalar<tensor_reduced_v> tensor_reduced;
  typedef typename GridTypeMapper<vtype>::scalar_object recurse_scalar_object;
  typedef iScalar<recurse_scalar_object> scalar_object;

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};

  // Scalar no action
  //  template<int Level> using tensor_reduce_level = typename iScalar<GridTypeMapper<vtype>::tensor_reduce_level<Level> >;

#ifndef TENSOR_IS_POD
  iScalar()=default;
  iScalar(scalar_type s) : _internal(s) {};// recurse down and hit the constructor for vector_type
  iScalar(const Zero &z){ *this = zero; };
#endif

    iScalar<vtype> & operator= (const Zero &hero){
      zeroit(*this);
      return *this;
    }
    iScalar<vtype> & operator= (const scalar_type s){
      _internal=s;
      return *this;
    }


    friend void zeroit(iScalar<vtype> &that){
        zeroit(that._internal);
    }
    friend void permute(iScalar<vtype> &out,const iScalar<vtype> &in,int permutetype){
      permute(out._internal,in._internal,permutetype);
    }

    // Unary negation
    friend inline iScalar<vtype> operator -(const iScalar<vtype> &r) {
        iScalar<vtype> ret;
        ret._internal= -r._internal;
        return ret;
    }
    // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
    inline iScalar<vtype> &operator *=(const iScalar<vtype> &r) {
        *this = (*this)*r;
        return *this;
    }
    inline iScalar<vtype> &operator -=(const iScalar<vtype> &r) {
        *this = (*this)-r;
        return *this;
    }
    inline iScalar<vtype> &operator +=(const iScalar<vtype> &r) {
        *this = (*this)+r;
        return *this;
    }
    
    inline vtype & operator ()(void) {
      return _internal;
    }

    inline const vtype & operator ()(void) const {
      return _internal;
    }
    //    inline vtype && operator ()(void) {
    //      return _internal;
    //    }

    operator ComplexD () const { return(TensorRemove(_internal)); };
    operator RealD () const { return(real(TensorRemove(_internal))); }


    template<class T,typename std::enable_if<isGridTensor<T>::notvalue, T>::type* = nullptr > inline auto operator = (T arg) -> iScalar<vtype>
    { 
      _internal = arg;
      return *this;
    }

};
///////////////////////////////////////////////////////////
// Allows to turn scalar<scalar<scalar<double>>>> back to double.
///////////////////////////////////////////////////////////
template<class T>     inline typename std::enable_if<isGridTensor<T>::notvalue, T>::type TensorRemove(T arg) { return arg;}
template<class vtype> inline auto TensorRemove(iScalar<vtype> arg) -> decltype(TensorRemove(arg._internal))
{
  return TensorRemove(arg._internal);
}
    
template<class vtype,int N> class iVector
{
public:
  vtype _internal[N];

  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;
  typedef typename GridTypeMapper<vtype>::scalar_object recurse_scalar_object;
  typedef iScalar<tensor_reduced_v> tensor_reduced;
  typedef iVector<recurse_scalar_object,N> scalar_object;


  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};

#ifndef TENSOR_IS_POD
  iVector(const Zero &z){ *this = zero; };
  iVector() =default;
#endif

    iVector<vtype,N> & operator= (const Zero &hero){
        zeroit(*this);
        return *this;
    }
    friend void zeroit(iVector<vtype,N> &that){
        for(int i=0;i<N;i++){
            zeroit(that._internal[i]);
        }
    }
    friend void permute(iVector<vtype,N> &out,const iVector<vtype,N> &in,int permutetype){
      for(int i=0;i<N;i++){
	permute(out._internal[i],in._internal[i],permutetype);
      }
    }
    // Unary negation
    friend inline iVector<vtype,N> operator -(const iVector<vtype,N> &r) {
        iVector<vtype,N> ret;
        for(int i=0;i<N;i++) ret._internal[i]= -r._internal[i];
        return ret;
    }
    // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
    inline iVector<vtype,N> &operator *=(const iScalar<vtype> &r) {
        *this = (*this)*r;
        return *this;
    }
    inline iVector<vtype,N> &operator -=(const iVector<vtype,N> &r) {
        *this = (*this)-r;
        return *this;
    }
    inline iVector<vtype,N> &operator +=(const iVector<vtype,N> &r) {
        *this = (*this)+r;
        return *this;
    }
    inline vtype & operator ()(int i) {
      return _internal[i];
    }
    inline const vtype & operator ()(int i) const {
      return _internal[i];
    }
    //    inline vtype && operator ()(int i) {
    //      return _internal[i];
    //    }
};
    
template<class vtype,int N> class iMatrix
{
public:
  vtype _internal[N][N];

  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;
  typedef typename GridTypeMapper<vtype>::scalar_object recurse_scalar_object;
  typedef iScalar<tensor_reduced_v> tensor_reduced;
  typedef iMatrix<recurse_scalar_object,N> scalar_object;

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};

#ifndef TENSOR_IS_POD
  iMatrix(const Zero &z){ *this = zero; };
  iMatrix() =default;
#endif

  iMatrix<vtype,N> & operator= (const Zero &hero){
    zeroit(*this);
    return *this;
  }
  template<class T,typename std::enable_if<isGridTensor<T>::notvalue, T>::type* = nullptr > inline auto operator = (T arg) -> iMatrix<vtype,N>
    { 
      zeroit(*this);
      for(int i=0;i<N;i++)
	_internal[i][i] = arg;
      return *this;
    }

  friend void zeroit(iMatrix<vtype,N> &that){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	zeroit(that._internal[i][j]);
    }}
  }
  friend void permute(iMatrix<vtype,N> &out,const iMatrix<vtype,N> &in,int permutetype){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	permute(out._internal[i][j],in._internal[i][j],permutetype);
    }}
  }
  // Unary negation
  friend inline iMatrix<vtype,N> operator -(const iMatrix<vtype,N> &r) {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	ret._internal[i][j]= -r._internal[i][j];
    }}
    return ret;
  }
  // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
  template<class T>
  inline iMatrix<vtype,N> &operator *=(const T &r) {
    *this = (*this)*r;
    return *this;
  }
  template<class T>
  inline iMatrix<vtype,N> &operator -=(const T &r) {
    *this = (*this)-r;
    return *this;
  }
  template<class T>
  inline iMatrix<vtype,N> &operator +=(const T &r) {
    *this = (*this)+r;
    return *this;
  }

  // returns an lvalue reference
  inline vtype & operator ()(int i,int j) {
    return _internal[i][j];
  }
  inline const vtype & operator ()(int i,int j) const {
    return _internal[i][j];
  }

  //  inline vtype && operator ()(int i,int j) {
  //    return _internal[i][j];
  //  }

};



}
#endif
