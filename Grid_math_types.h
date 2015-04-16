#ifndef GRID_MATH_TYPES_H
#define GRID_MATH_TYPES_H

#include <Grid_math_type_mapper.h>
#include <type_traits>

namespace Grid {

  // First some of my own traits
  template<typename T> struct isGridTensor {
    static const bool value = true;
    static const bool notvalue = false;
  };

  template<> struct isGridTensor<RealD > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<RealF > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<ComplexD > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<ComplexF > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<Integer > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vRealD > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vRealF > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vComplexD > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vComplexF > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vInteger > {
    static const bool value = false;
    static const bool notvalue = true;
  };

  // Match the index
  template<typename T,int Level> struct matchGridTensorIndex {
    static const bool value = (Level==T::TensorLevel);
    static const bool notvalue = (Level!=T::TensorLevel);
  };


///////////////////////////////////////////////////
// Scalar, Vector, Matrix objects.
// These can be composed to form tensor products of internal indices.
///////////////////////////////////////////////////

  // Terminates the recursion for temoval of all Grids tensors
  inline vRealD    TensorRemove(vRealD    arg){ return arg;}
  inline vRealF    TensorRemove(vRealF    arg){ return arg;}
  inline vComplexF TensorRemove(vComplexF arg){ return arg;}
  inline vComplexD TensorRemove(vComplexD arg){ return arg;}
  inline vInteger  TensorRemove(vInteger  arg){ return arg;}

template<class vtype> class iScalar
{
public:
  vtype _internal;

  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;
  typedef iScalar<tensor_reduced_v> tensor_reduced;

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};

  // Scalar no action
  //  template<int Level> using tensor_reduce_level = typename iScalar<GridTypeMapper<vtype>::tensor_reduce_level<Level> >;

    iScalar(){};
    
    iScalar(scalar_type s) : _internal(s) {};// recurse down and hit the constructor for vector_type

    iScalar(Zero &z){ *this = zero; };

    iScalar<vtype> & operator= (const Zero &hero){
        zeroit(*this);
        return *this;
    }
    friend void zeroit(iScalar<vtype> &that){
        zeroit(that._internal);
    }
    friend void permute(iScalar<vtype> &out,const iScalar<vtype> &in,int permutetype){
      permute(out._internal,in._internal,permutetype);
    }
    friend void extract(const iScalar<vtype> &in,std::vector<scalar_type *> &out){
      extract(in._internal,out); // extract advances the pointers in out
    }
    friend void merge(iScalar<vtype> &in,std::vector<scalar_type *> &out){
      merge(in._internal,out); // extract advances the pointers in out
    }
    friend inline iScalar<vtype>::vector_type TensorRemove(iScalar<vtype> arg)
    {
      return TensorRemove(arg._internal);
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
};
    
template<class vtype,int N> class iVector
{
public:
  vtype _internal[N];

  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};
  typedef iScalar<tensor_reduced_v> tensor_reduced;

    iVector(Zero &z){ *this = zero; };
    iVector() {};
    iVector<vtype,N> & operator= (Zero &hero){
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
    friend void extract(const iVector<vtype,N> &in,std::vector<scalar_type *> &out){
      for(int i=0;i<N;i++){
	extract(in._internal[i],out);// extract advances pointers in out
      }
    }
    friend void merge(iVector<vtype,N> &in,std::vector<scalar_type *> &out){
      for(int i=0;i<N;i++){
	merge(in._internal[i],out);// extract advances pointers in out
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
};
    
template<class vtype,int N> class iMatrix
{
public:
  vtype _internal[N][N];

  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;

  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};
  typedef iScalar<tensor_reduced_v> tensor_reduced;

  iMatrix(Zero &z){ *this = zero; };
  iMatrix() {};
  iMatrix<vtype,N> & operator= (Zero &hero){
    zeroit(*this);
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
  friend void extract(const iMatrix<vtype,N> &in,std::vector<scalar_type *> &out){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	extract(in._internal[i][j],out);// extract advances pointers in out
    }}
  }
  friend void merge(iMatrix<vtype,N> &in,std::vector<scalar_type *> &out){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	merge(in._internal[i][j],out);// extract advances pointers in out
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

};

    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// ADD         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    

// ADD is simple for now; cannot mix types and straightforward template
// Scalar +/- Scalar
// Vector +/- Vector
// Matrix +/- Matrix
template<class vtype,class ltype,class rtype> inline void add(iScalar<vtype> * __restrict__ ret,
                                                              const iScalar<ltype> * __restrict__ lhs,
                                                              const iScalar<rtype> * __restrict__ rhs)
{
    add(&ret->_internal,&lhs->_internal,&rhs->_internal);
}
template<class vtype,class ltype,class rtype,int N> inline void add(iVector<vtype,N> * __restrict__ ret,
                                                                    const iVector<ltype,N> * __restrict__ lhs,
                                                                    const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c=0;c<N;c++){
        ret->_internal[c]=lhs->_internal[c]+rhs->_internal[c];
    }
    return;
}
template<class vtype,class ltype,class rtype, int N> inline  void add(iMatrix<vtype,N> * __restrict__ ret,
                                                                      const iMatrix<ltype,N> * __restrict__ lhs,
                                                                      const iMatrix<rtype,N> * __restrict__ rhs)
{
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        add(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> inline  void add(iMatrix<vtype,N> * __restrict__ ret,
                                                                      const iScalar<ltype>   * __restrict__ lhs,
                                                                      const iMatrix<rtype,N> * __restrict__ rhs)
{
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        add(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> inline  void add(iMatrix<vtype,N> * __restrict__ ret,
                                                                      const iMatrix<ltype,N> * __restrict__ lhs,
                                                                      const iScalar<rtype>   * __restrict__ rhs)
{
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        if ( c1==c2)
            add(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
        else
            ret->_internal[c1][c2]=lhs->_internal[c1][c2];
    }}
    return;
}
// Need to figure multi-precision.
template<class Mytype>  Mytype timesI(Mytype &r)
{
    iScalar<Complex> i;
    i._internal = Complex(0,1);
    return r*i;
}

                // + operator for scalar, vector, matrix
template<class ltype,class rtype>
//inline auto operator + (iScalar<ltype>& lhs,iScalar<rtype>&& rhs) -> iScalar<decltype(lhs._internal + rhs._internal)>
inline auto operator + (const iScalar<ltype>& lhs,const iScalar<rtype>& rhs) -> iScalar<decltype(lhs._internal + rhs._internal)>
{
    typedef iScalar<decltype(lhs._internal+rhs._internal)> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator + (const iVector<ltype,N>& lhs,const iVector<rtype,N>& rhs) ->iVector<decltype(lhs._internal[0]+rhs._internal[0]),N>
{
    typedef iVector<decltype(lhs._internal[0]+rhs._internal[0]),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator + (const iMatrix<ltype,N>& lhs,const iMatrix<rtype,N>& rhs) ->iMatrix<decltype(lhs._internal[0][0]+rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]+rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator + (const iScalar<ltype>& lhs,const iMatrix<rtype,N>& rhs)->iMatrix<decltype(lhs._internal+rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal+rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator + (const iMatrix<ltype,N>& lhs,const iScalar<rtype>& rhs)->iMatrix<decltype(lhs._internal[0][0]+rhs._internal),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]+rhs._internal),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
}


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// SUB         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    

// SUB is simple for now; cannot mix types and straightforward template
// Scalar +/- Scalar
// Vector +/- Vector
// Matrix +/- Matrix
// Matrix /- scalar
template<class vtype,class ltype,class rtype> inline void sub(iScalar<vtype> * __restrict__ ret,
                                                              const iScalar<ltype> * __restrict__ lhs,
                                                              const iScalar<rtype> * __restrict__ rhs)
{
    sub(&ret->_internal,&lhs->_internal,&rhs->_internal);
}

template<class vtype,class ltype,class rtype,int N> inline void sub(iVector<vtype,N> * __restrict__ ret,
                                                                    const iVector<ltype,N> * __restrict__ lhs,
                                                                    const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c=0;c<N;c++){
        ret->_internal[c]=lhs->_internal[c]-rhs->_internal[c];
    }
    return;
}
template<class vtype,class ltype,class rtype, int N> inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iMatrix<ltype,N> * __restrict__ lhs,
                                                                     const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        sub(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iScalar<ltype> * __restrict__ lhs,
                                                                     const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        if ( c1!=c2) {
            sub(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
        } else {
            // Fails -- need unary minus. Catalogue other unops?
            ret->_internal[c1][c2]=zero;
            ret->_internal[c1][c2]=ret->_internal[c1][c2]-rhs->_internal[c1][c2];

        }
    }}
    return;
}
template<class vtype,class ltype,class rtype, int N> inline void sub(iMatrix<vtype,N> * __restrict__ ret,
                                                                     const iMatrix<ltype,N> * __restrict__ lhs,
                                                                     const iScalar<rtype> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        if ( c1!=c2)
            sub(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
        else
            ret->_internal[c1][c2]=lhs->_internal[c1][c2];
    }}
    return;
}

template<class v> void vprefetch(const iScalar<v> &vv)
{
  vprefetch(vv._internal);
}
template<class v,int N> void vprefetch(const iVector<v,N> &vv)
{
  for(int i=0;i<N;i++){
    vprefetch(vv._internal[i]);
  }
}
template<class v,int N> void vprefetch(const iMatrix<v,N> &vv)
{
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    vprefetch(vv._internal[i][j]);
  }}
}

    // - operator for scalar, vector, matrix
template<class ltype,class rtype> inline auto
operator - (const iScalar<ltype>& lhs, const iScalar<rtype>& rhs) -> iScalar<decltype(lhs._internal - rhs._internal)>
{
    typedef iScalar<decltype(lhs._internal-rhs._internal)> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator - (const iVector<ltype,N>& lhs,const iVector<rtype,N>& rhs) ->iVector<decltype(lhs._internal[0]-rhs._internal[0]),N>
{
    typedef iVector<decltype(lhs._internal[0]-rhs._internal[0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator - (const iMatrix<ltype,N>& lhs,const iMatrix<rtype,N>& rhs) ->iMatrix<decltype(lhs._internal[0][0]-rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]-rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator - (const iScalar<ltype>& lhs,const iMatrix<rtype,N>& rhs)->iMatrix<decltype(lhs._internal-rhs._internal[0][0]),N>
{
    typedef iMatrix<decltype(lhs._internal-rhs._internal[0][0]),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}
template<class ltype,class rtype,int N>
inline auto operator - (const iMatrix<ltype,N>& lhs,const iScalar<rtype>& rhs)->iMatrix<decltype(lhs._internal[0][0]-rhs._internal),N>
{
    typedef iMatrix<decltype(lhs._internal[0][0]-rhs._internal),N> ret_t;
    ret_t ret;
    sub(&ret,&lhs,&rhs);
    return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// MAC         ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////
    // Legal multiplication table
    ///////////////////////////
    // scal x scal = scal
    // mat x  mat  = mat
    // mat  x scal = mat
    // scal x mat  = mat
    // mat  x vec  = vec
    // vec  x scal = vec
    // scal x vec  = vec
    ///////////////////////////
template<class rtype,class vtype,class mtype>
inline  void mac(iScalar<rtype> * __restrict__ ret,const iScalar<vtype> * __restrict__ lhs,const iScalar<mtype> * __restrict__ rhs)
{
    mac(&ret->_internal,&lhs->_internal,&rhs->_internal);
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iMatrix<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
    for(int c3=0;c3<N;c3++){
        mac(&ret->_internal[c1][c2],&lhs->_internal[c1][c3],&rhs->_internal[c3][c2]);
    }}}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iMatrix<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iScalar<rtype> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mac(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
    }}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iMatrix<rrtype,N> * __restrict__ ret,const iScalar<ltype> * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mac(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
    }}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iVector<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mac(&ret->_internal[c1],&lhs->_internal[c1][c2],&rhs->_internal[c2]);
    }}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iVector<rrtype,N> * __restrict__ ret,const iScalar<ltype> * __restrict__ lhs,const iVector<rtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mac(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mac(iVector<rrtype,N> * __restrict__ ret,const iVector<ltype,N> * __restrict__ lhs,const iScalar<rtype> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mac(&ret->_internal[c1],&lhs->_internal[c1],&rhs->_internal);
    }
    return;
}

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// MUL         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    
template<class rtype,class vtype,class mtype>
inline void mult(iScalar<rtype> * __restrict__ ret,const iScalar<mtype> * __restrict__ lhs,const iScalar<vtype> * __restrict__ rhs){
    mult(&ret->_internal,&lhs->_internal,&rhs->_internal);
}

template<class rrtype,class ltype,class rtype,int N>
inline void mult(iMatrix<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1][c2],&lhs->_internal[c1][0],&rhs->_internal[0][c2]);
        for(int c3=1;c3<N;c3++){
            mac(&ret->_internal[c1][c2],&lhs->_internal[c1][c3],&rhs->_internal[c3][c2]);
        }
    }}
    return;
}
template<class rrtype,class ltype,class rtype,int N>
inline void mult(iMatrix<rrtype,N> * __restrict__ ret,const iMatrix<ltype,N> * __restrict__ lhs,const iScalar<rtype> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal);
    }}
    return;
}

template<class rrtype,class ltype,class rtype, int N>
inline void mult(iMatrix<rrtype,N> * __restrict__ ret,const iScalar<ltype>   * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs){
    for(int c2=0;c2<N;c2++){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1][c2],&lhs->_internal,&rhs->_internal[c1][c2]);
    }}
    return;
}
// Matrix left multiplies vector
template<class rtype,class vtype,class mtype,int N>
inline void mult(iVector<rtype,N> * __restrict__ ret,const iMatrix<mtype,N> * __restrict__ lhs,const iVector<vtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1],&lhs->_internal[c1][0],&rhs->_internal[0]);
        for(int c2=1;c2<N;c2++){
            mac(&ret->_internal[c1],&lhs->_internal[c1][c2],&rhs->_internal[c2]);
        }
    }
    return;
}
template<class rtype,class vtype,class mtype,int N>
inline void mult(iVector<rtype,N> * __restrict__ ret,
                 const iScalar<mtype>   * __restrict__ lhs,
                 const iVector<vtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
}
template<class rtype,class vtype,class mtype,int N>
inline void mult(iVector<rtype,N> * __restrict__ ret,
                 const iVector<vtype,N> * __restrict__ rhs,
                 const iScalar<mtype> * __restrict__ lhs){
    mult(ret,lhs,rhs);
}
    


template<class rtype,class vtype,class mtype,int N> inline
iVector<rtype,N> operator * (const iMatrix<mtype,N>& lhs,const iVector<vtype,N>& rhs)
{
    iVector<rtype,N> ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}

template<class rtype,class vtype,class mtype,int N> inline
iVector<rtype,N> operator * (const iScalar<mtype>& lhs,const iVector<vtype,N>& rhs)
{
    iVector<rtype,N> ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}

template<class rtype,class vtype,class mtype,int N> inline
iVector<rtype,N> operator * (const iVector<mtype,N>& lhs,const iScalar<vtype>& rhs)
{
    iVector<rtype,N> ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}
    
    //////////////////////////////////////////////////////////////////
    // Glue operators to mult routines. Must resolve return type cleverly from typeof(internal)
    // since nesting matrix<scalar> x matrix<matrix>-> matrix<matrix>
    // while         matrix<scalar> x matrix<scalar>-> matrix<scalar>
    // so return type depends on argument types in nasty way.
    //////////////////////////////////////////////////////////////////
    // scal x scal = scal
    // mat x  mat  = mat
    // mat  x scal = mat
    // scal x mat  = mat
    // mat  x vec  = vec
    // vec  x scal = vec
    // scal x vec  = vec
    //
    // We can special case scalar_type ??
template<class l,class r>
inline auto operator * (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(lhs._internal * rhs._internal)>
{
    typedef iScalar<decltype(lhs._internal*rhs._internal)> ret_t;
    ret_t ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iMatrix<decltype(lhs._internal[0][0]*rhs._internal[0][0]),N>
{
    typedef decltype(lhs._internal[0][0]*rhs._internal[0][0]) ret_t;
    iMatrix<ret_t,N> ret;
    mult(&ret,&lhs,&rhs);
    return ret;
}
template<class l,class r, int N> inline
auto operator * (const iMatrix<r,N>& lhs,const iScalar<l>& rhs) -> iMatrix<decltype(lhs._internal[0][0]*rhs._internal),N>
{
    typedef decltype(lhs._internal[0][0]*rhs._internal) ret_t;
        
    iMatrix<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mult(&ret._internal[c1][c2],&lhs._internal[c1][c2],&rhs._internal);
    }}
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iScalar<l>& lhs,const iMatrix<r,N>& rhs) -> iMatrix<decltype(lhs._internal*rhs._internal[0][0]),N>
{
    typedef decltype(lhs._internal*rhs._internal[0][0]) ret_t;
    iMatrix<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mult(&ret._internal[c1][c2],&lhs._internal,&rhs._internal[c1][c2]);
    }}
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iMatrix<l,N>& lhs,const iVector<r,N>& rhs) -> iVector<decltype(lhs._internal[0][0]*rhs._internal[0]),N>
{
    typedef decltype(lhs._internal[0][0]*rhs._internal[0]) ret_t;
    iVector<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal[c1][0],&rhs._internal[0]);
        for(int c2=1;c2<N;c2++){
            mac(&ret._internal[c1],&lhs._internal[c1][c2],&rhs._internal[c2]);
        }
    }
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iScalar<l>& lhs,const iVector<r,N>& rhs) -> iVector<decltype(lhs._internal*rhs._internal[0]),N>
{
    typedef decltype(lhs._internal*rhs._internal[0]) ret_t;
    iVector<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal,&rhs._internal[c1]);
    }
    return ret;
}
template<class l,class r,int N> inline
auto operator * (const iVector<l,N>& lhs,const iScalar<r>& rhs) -> iVector<decltype(lhs._internal[0]*rhs._internal),N>
{
    typedef decltype(lhs._internal[0]*rhs._internal) ret_t;
    iVector<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal[c1],&rhs._internal);
    }
    return ret;
}
//////////////////////////////////////////////////////////////////////////////////////////
// Must support native C++ types Integer, Complex, Real
//////////////////////////////////////////////////////////////////////////////////////////

// multiplication by fundamental scalar type
template<class l,int N> inline iScalar<l> operator * (const iScalar<l>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iScalar<l>::tensor_reduced srhs(rhs);
  return lhs*srhs;
}
template<class l,int N> inline iScalar<l> operator * (const typename iScalar<l>::scalar_type lhs,const iScalar<l>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iVector<l,N> operator * (const iVector<l,N>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iVector<l,N>::tensor_reduced srhs(rhs);
  return lhs*srhs;
}
template<class l,int N> inline iVector<l,N> operator * (const typename iScalar<l>::scalar_type lhs,const iVector<l,N>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iMatrix<l,N> operator * (const iMatrix<l,N>& lhs,const typename iScalar<l>::scalar_type &rhs) 
{
  typename iMatrix<l,N>::tensor_reduced srhs(rhs);
  return lhs*srhs;
}
template<class l,int N> inline iMatrix<l,N> operator * (const typename iScalar<l>::scalar_type & lhs,const iMatrix<l,N>& rhs) {  return rhs*lhs; }

////////////////////////////////////////////////////////////////////
// Double support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator * (const iScalar<l>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l> inline iScalar<l> operator * (double lhs,const iScalar<l>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iVector<l,N> operator * (const iVector<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iVector<l,N> operator * (double lhs,const iVector<l,N>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iMatrix<l,N> operator * (const iMatrix<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iMatrix<l,N> operator * (double lhs,const iMatrix<l,N>& rhs) {  return rhs*lhs; }

////////////////////////////////////////////////////////////////////
// Integer support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator * (const iScalar<l>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l> inline iScalar<l> operator * (Integer lhs,const iScalar<l>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iVector<l,N> operator * (const iVector<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iVector<l,N> operator * (Integer lhs,const iVector<l,N>& rhs) {  return rhs*lhs; }

template<class l,int N> inline iMatrix<l,N> operator * (const iMatrix<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs*srhs;
}
template<class l,int N> inline iMatrix<l,N> operator * (Integer lhs,const iMatrix<l,N>& rhs) {  return rhs*lhs; }



///////////////////////////////////////////////////////////////////////////////////////////////
// addition by fundamental scalar type applies to matrix(down diag) and scalar
///////////////////////////////////////////////////////////////////////////////////////////////
template<class l,int N> inline iScalar<l> operator + (const iScalar<l>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iScalar<l>::tensor_reduced srhs(rhs);
  return lhs+srhs;
}
template<class l,int N> inline iScalar<l> operator + (const typename iScalar<l>::scalar_type lhs,const iScalar<l>& rhs) {  return rhs+lhs; }

template<class l,int N> inline iMatrix<l,N> operator + (const iMatrix<l,N>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iMatrix<l,N>::tensor_reduced srhs(rhs);
  return lhs+srhs;
}
template<class l,int N> inline iMatrix<l,N> operator + (const typename iScalar<l>::scalar_type lhs,const iMatrix<l,N>& rhs) {  return rhs+lhs; }

////////////////////////////////////////////////////////////////////
// Double support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator + (const iScalar<l>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs+srhs;
}
template<class l> inline iScalar<l> operator + (double lhs,const iScalar<l>& rhs) {  return rhs+lhs; }

template<class l,int N> inline iMatrix<l,N> operator + (const iMatrix<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs+srhs;
}
template<class l,int N> inline iMatrix<l,N> operator + (double lhs,const iMatrix<l,N>& rhs) {  return rhs+lhs; }

////////////////////////////////////////////////////////////////////
// Integer support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator + (const iScalar<l>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs+srhs;
}
template<class l> inline iScalar<l> operator + (Integer lhs,const iScalar<l>& rhs) {  return rhs+lhs; }

template<class l,int N> inline iMatrix<l,N> operator + (const iMatrix<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs+srhs;
}
template<class l,int N> inline iMatrix<l,N> operator + (Integer lhs,const iMatrix<l,N>& rhs) {  return rhs+lhs; }


///////////////////////////////////////////////////////////////////////////////////////////////
// subtraction of fundamental scalar type applies to matrix(down diag) and scalar
///////////////////////////////////////////////////////////////////////////////////////////////
template<class l,int N> inline iScalar<l> operator - (const iScalar<l>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iScalar<l>::tensor_reduced srhs(rhs);
  return lhs-srhs;
}
template<class l,int N> inline iScalar<l> operator - (const typename iScalar<l>::scalar_type lhs,const iScalar<l>& rhs) 
{
  typename iScalar<l>::tensor_reduced slhs(lhs);
  return slhs-rhs;
}

template<class l,int N> inline iMatrix<l,N> operator - (const iMatrix<l,N>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iScalar<l>::tensor_reduced srhs(rhs);
  return lhs-srhs;
}
template<class l,int N> inline iMatrix<l,N> operator - (const typename iScalar<l>::scalar_type lhs,const iMatrix<l,N>& rhs) 
{
  typename iScalar<l>::tensor_reduced slhs(lhs);
  return slhs-rhs;
}

////////////////////////////////////////////////////////////////////
// Double support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator - (const iScalar<l>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs-srhs;
}
template<class l> inline iScalar<l> operator - (double lhs,const iScalar<l>& rhs) 
{
  typename iScalar<l>::scalar_type t(lhs);
  typename iScalar<l>::tensor_reduced slhs(t);
  return slhs-rhs;
}

template<class l,int N> inline iMatrix<l,N> operator - (const iMatrix<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs-srhs;
}
template<class l,int N> inline iMatrix<l,N> operator - (double lhs,const iMatrix<l,N>& rhs) 
{
  typename iScalar<l>::scalar_type t(lhs);
  typename iScalar<l>::tensor_reduced slhs(t);
  return slhs-rhs;
}

////////////////////////////////////////////////////////////////////
// Integer support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////
template<class l> inline iScalar<l> operator - (const iScalar<l>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs-srhs;
}
template<class l> inline iScalar<l> operator - (Integer lhs,const iScalar<l>& rhs) 
{
  typename iScalar<l>::scalar_type t(lhs);
  typename iScalar<l>::tensor_reduced slhs(t);
  return slhs-rhs;
}
template<class l,int N> inline iMatrix<l,N> operator - (const iMatrix<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t(rhs);
  typename iScalar<l>::tensor_reduced srhs(t);
  return lhs-srhs;
}
template<class l,int N> inline iMatrix<l,N> operator - (Integer lhs,const iMatrix<l,N>& rhs) 
{
  typename iScalar<l>::scalar_type t(lhs);
  typename iScalar<l>::tensor_reduced slhs(t);
  return slhs-rhs;
}




    ///////////////////////////////////////////////////////////////////////////////////////
    // innerProduct Scalar x Scalar -> Scalar
    // innerProduct Vector x Vector -> Scalar
    // innerProduct Matrix x Matrix -> Scalar
    ///////////////////////////////////////////////////////////////////////////////////////
    template<class l,class r,int N> inline
    auto innerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0],rhs._internal[0]))>
    {
        typedef decltype(innerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
        iScalar<ret_t> ret=zero;
        for(int c1=0;c1<N;c1++){
            ret._internal += innerProduct(lhs._internal[c1],rhs._internal[c1]);
        }
        return ret;
    }
    template<class l,class r,int N> inline
    auto innerProduct (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0]))>
    {
        typedef decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
        iScalar<ret_t> ret=zero;
        iScalar<ret_t> tmp;
        for(int c1=0;c1<N;c1++){
        for(int c2=0;c2<N;c2++){
	  ret._internal+=innerProduct(lhs._internal[c1][c2],rhs._internal[c1][c2]);
        }}
        return ret;
    }
    template<class l,class r> inline
    auto innerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(innerProduct(lhs._internal,rhs._internal))>
    {
        typedef decltype(innerProduct(lhs._internal,rhs._internal)) ret_t;
        iScalar<ret_t> ret;
        ret._internal = innerProduct(lhs._internal,rhs._internal);
        return ret;
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    // outerProduct Scalar x Scalar -> Scalar
    //              Vector x Vector -> Matrix
    ///////////////////////////////////////////////////////////////////////////////////////

template<class l,class r,int N> inline
auto outerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iMatrix<decltype(outerProduct(lhs._internal[0],rhs._internal[0])),N>
{
    typedef decltype(outerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
    iMatrix<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = outerProduct(lhs._internal[c1],rhs._internal[c2]);
    }}
    return ret;
}
template<class l,class r> inline
auto outerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(outerProduct(lhs._internal,rhs._internal))>
{
    typedef decltype(outerProduct(lhs._internal,rhs._internal)) ret_t;
    iScalar<ret_t> ret;
    ret._internal = outerProduct(lhs._internal,rhs._internal);
    return ret;
}

inline ComplexF outerProduct(const ComplexF &l, const ComplexF& r)
{
  return l*r;
}
inline ComplexD outerProduct(const ComplexD &l, const ComplexD& r)
{
  return l*r;
}
inline RealF outerProduct(const RealF &l, const RealF& r)
{
  return l*r;
}
inline RealD outerProduct(const RealD &l, const RealD& r)
{
  return l*r;
}
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// CONJ         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
 
// Conj function for scalar, vector, matrix
template<class vtype> inline iScalar<vtype> conj(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = conj(r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> conj(const iVector<vtype,N>&r)
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i] = conj(r._internal[i]);
  }
  return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> conj(const iMatrix<vtype,N>&r)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    ret._internal[i][j] = conj(r._internal[i][j]);
  }}
  return ret;
}

// Adj function for scalar, vector, matrix
template<class vtype> inline iScalar<vtype> adj(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = adj(r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> adj(const iVector<vtype,N>&r)
{
    iVector<vtype,N> ret;
    for(int i=0;i<N;i++){
        ret._internal[i] = adj(r._internal[i]);
    }
    return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> adj(const iMatrix<vtype,N> &arg)
{
    iMatrix<vtype,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2]=adj(arg._internal[c2][c1]);
    }}
    return ret;
}
/////////////////////////////////////////////////////////////////
// Transpose all indices
/////////////////////////////////////////////////////////////////

inline ComplexD transpose(ComplexD &rhs){  return rhs;}
inline ComplexF transpose(ComplexF &rhs){  return rhs;}
inline RealD transpose(RealD &rhs){  return rhs;}
inline RealF transpose(RealF &rhs){  return rhs;}

template<class vtype,int N>
  inline typename std::enable_if<isGridTensor<vtype>::value, iMatrix<vtype,N> >::type 
  transpose(iMatrix<vtype,N> arg)
  {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	ret._internal[i][j] = transpose(arg._internal[j][i]); // NB recurses
      }}
    return ret;
  }
template<class vtype,int N>
  inline typename std::enable_if<isGridTensor<vtype>::notvalue, iMatrix<vtype,N> >::type 
  transpose(iMatrix<vtype,N> arg)
  {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	ret._internal[i][j] = arg._internal[j][i]; // Stop recursion if not a tensor type
      }}
    return ret;
  }

template<class vtype,int N>
  inline typename std::enable_if<isGridTensor<vtype>::value, iScalar<vtype> >::type 
  transpose(iScalar<vtype> arg)
  {
    iScalar<vtype> ret;
    ret._internal = transpose(arg._internal); // NB recurses
    return ret;
  }

template<class vtype,int N>
  inline typename std::enable_if<isGridTensor<vtype>::notvalue, iScalar<vtype> >::type 
  transpose(iScalar<vtype> arg)
  {
    iScalar<vtype> ret;
    ret._internal = arg._internal; // NB recursion stops
    return ret;
  }

////////////////////////////////////////////////////////////////////////////////////////////
// Transpose a specific index; instructive to compare this style of recursion termination
// to that of adj; which is easiers?
////////////////////////////////////////////////////////////////////////////////////////////
template<int Level,class vtype,int N> inline 
  typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::value, iMatrix<vtype,N> >::type 
transposeIndex (const iMatrix<vtype,N> &arg)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = arg._internal[j][i]; 
  }}
  return ret;
}
// or not
template<int Level,class vtype,int N> inline 
typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue, iMatrix<vtype,N> >::type 
transposeIndex (const iMatrix<vtype,N> &arg)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = transposeIndex<Level>(arg._internal[i][j]); 
  }}
  return ret;
}
template<int Level,class vtype> inline 
typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue, iScalar<vtype> >::type 
transposeIndex (const iScalar<vtype> &arg)
{
  return transposeIndex<Level>(arg._internal);
}

//////////////////////////////////////////////////////////////////
// Traces: both all indices and a specific index 
/////////////////////////////////////////////////////////////////

inline ComplexF trace( const ComplexF &arg){    return arg;}
inline ComplexD trace( const ComplexD &arg){    return arg;}
inline RealF trace( const RealF &arg){    return arg;}
inline RealD trace( const RealD &arg){    return arg;}

template<int Level> inline ComplexF traceIndex(const ComplexF arg) { return arg;}
template<int Level> inline ComplexD traceIndex(const ComplexD arg) { return arg;}
template<int Level> inline RealF traceIndex(const RealF arg) { return arg;}
template<int Level> inline RealD traceIndex(const RealD arg) { return arg;}

template<class vtype,int N>
inline auto trace(const iMatrix<vtype,N> &arg) -> iScalar<decltype(trace(arg._internal[0][0]))>
{
    iScalar<decltype( trace(arg._internal[0][0] )) > ret;
    zeroit(ret._internal);
    for(int i=0;i<N;i++){
        ret._internal=ret._internal+trace(arg._internal[i][i]);
    }
    return ret;
}
template<class vtype>
inline auto trace(const iScalar<vtype> &arg) -> iScalar<decltype(trace(arg._internal))>
{
    iScalar<decltype(trace(arg._internal))> ret;
    ret._internal=trace(arg._internal);
    return ret;
}

// Specific indices.
template<int Level,class vtype> inline 
auto traceIndex(const iScalar<vtype> &arg) -> iScalar<decltype(traceIndex<Level>(arg._internal)) >
{
  iScalar<decltype(traceIndex<Level>(arg._internal))> ret;
  ret._internal = traceIndex<Level>(arg._internal);
  return ret;
}

// If we hit the right index, return scalar and trace it with no further recursion
template<int Level,class vtype,int N> inline 
auto traceIndex(const iMatrix<vtype,N> &arg) ->
  typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::value,  // Index matches
                                                    iScalar<vtype> >::type                              // return scalar
{
  iScalar<vtype> ret;
  zeroit(ret._internal);
  for(int i=0;i<N;i++){
    ret._internal = ret._internal + arg._internal[i][i];
  }
  return ret;
}

// not this level, so recurse
template<int Level,class vtype,int N> inline 
auto traceIndex(const iMatrix<vtype,N> &arg) ->
  typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue,// No index match
         iMatrix<decltype(traceIndex<Level>(arg._internal[0][0])),N> >::type     // return matrix
{
  iMatrix<decltype(traceIndex<Level>(arg._internal[0][0])),N> ret;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    ret._internal[i][j] = traceIndex<Level>(arg._internal[i][j]);
  }}
  return ret;
}

/////////////////////////////////////////////////////////////////
// Can only take the real/imag part of scalar objects, since
// lattice objects of different complex nature are non-conformable.
/////////////////////////////////////////////////////////////////
template<class itype> inline auto real(const iScalar<itype> &z) -> iScalar<decltype(real(z._internal))>
{
    iScalar<decltype(real(z._internal))> ret;
    ret._internal = real(z._internal);
    return ret;
}
template<class itype,int N> inline auto real(const iMatrix<itype,N> &z) -> iMatrix<decltype(real(z._internal[0][0])),N>
{
    iMatrix<decltype(real(z._internal[0][0])),N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = real(z._internal[c1][c2]);
    }}
    return ret;
}
template<class itype,int N> inline auto real(const iVector<itype,N> &z) -> iVector<decltype(real(z._internal[0])),N>
{
    iVector<decltype(real(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = real(z._internal[c1]);
    }
    return ret;
}
    
template<class itype> inline auto imag(const iScalar<itype> &z) -> iScalar<decltype(imag(z._internal))>
{
    iScalar<decltype(imag(z._internal))> ret;
    ret._internal = imag(z._internal);
    return ret;
}
template<class itype,int N> inline auto imag(const iMatrix<itype,N> &z) -> iMatrix<decltype(imag(z._internal[0][0])),N>
{
    iMatrix<decltype(imag(z._internal[0][0])),N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = imag(z._internal[c1][c2]);
    }}
    return ret;
}
template<class itype,int N> inline auto imag(const iVector<itype,N> &z) -> iVector<decltype(imag(z._internal[0])),N>
{
    iVector<decltype(imag(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = imag(z._internal[c1]);
    }
    return ret;
}

};
    
#endif
