#ifndef GRID_MATH_TYPES_H
#define GRID_MATH_TYPES_H
namespace dpo {



//////////////////////////////////////////////////////////////////////////////////
// Want to recurse: GridTypeMapper<Matrix<vComplexD> >::scalar_type == ComplexD.
//////////////////////////////////////////////////////////////////////////////////

  template <class T> class GridTypeMapper {
  public:
    typedef typename T::scalar_type scalar_type;
    typedef typename T::vector_type vector_type;
  };

  template<> class GridTypeMapper<RealF> {
  public:
    typedef RealF scalar_type;
    typedef RealF vector_type;
  };
  template<> class GridTypeMapper<RealD> {
  public:
    typedef RealD scalar_type;
    typedef RealD vector_type;
  };
  template<> class GridTypeMapper<ComplexF> {
  public:
    typedef ComplexF scalar_type;
    typedef ComplexF vector_type;
  };
  template<> class GridTypeMapper<ComplexD> {
  public:
    typedef ComplexD scalar_type;
    typedef ComplexD vector_type;
  };

  template<> class GridTypeMapper<vRealF> {
  public:
    typedef RealF  scalar_type;
    typedef vRealF vector_type;
  };
  template<> class GridTypeMapper<vRealD> {
  public:
    typedef RealD  scalar_type;
    typedef vRealD vector_type;
  };
  template<> class GridTypeMapper<vComplexF> {
  public:
    typedef ComplexF  scalar_type;
    typedef vComplexF vector_type;
  };
  template<> class GridTypeMapper<vComplexD> {
  public:
    typedef ComplexD  scalar_type;
    typedef vComplexD vector_type;
  };



///////////////////////////////////////////////////
// Scalar, Vector, Matrix objects.
// These can be composed to form tensor products of internal indices.
///////////////////////////////////////////////////
    
template<class vtype> class iScalar
{
public:
  SIMDalign vtype _internal;

  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;

    iScalar(){};
    iScalar(Zero &z){ *this = zero; };
    iScalar<vtype> & operator= (const Zero &hero){
        zeroit(*this);
        return *this;
    }
    friend void zeroit(iScalar<vtype> &that){
        zeroit(that._internal);
    }
    friend void permute(iScalar<vtype> &out,iScalar<vtype> &in,int permutetype){
      permute(out._internal,in._internal,permutetype);
    }
    friend void extract(iScalar<vtype> &in,std::vector<scalar_type *> &out){
      extract(in._internal,out); // extract advances the pointers in out
    }
    friend void merge(iScalar<vtype> &in,std::vector<scalar_type *> &out){
      merge(in._internal,out); // extract advances the pointers in out
    }
    // Unary negation
    friend inline iScalar<vtype> operator -(const iScalar<vtype> &r) {
        iScalar<vtype> ret;
        ret._internal= -r._internal;
        return ret;
    }
    // *=,+=,-= operators
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
  SIMDalign vtype _internal[N];

  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;

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
    friend void permute(iVector<vtype,N> &out,iVector<vtype,N> &in,int permutetype){
      for(int i=0;i<N;i++){
	permute(out._internal[i],in._internal[i],permutetype);
      }
    }
    friend void extract(iVector<vtype,N> &in,std::vector<scalar_type *> &out){
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
    // *=,+=,-= operators
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
  SIMDalign    vtype _internal[N][N];

  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;

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
    friend void permute(iMatrix<vtype,N> &out,iMatrix<vtype,N> &in,int permutetype){
      for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	permute(out._internal[i][j],in._internal[i][j],permutetype);
      }}
    }
    friend void extract(iMatrix<vtype,N> &in,std::vector<scalar_type *> &out){
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
    // *=,+=,-= operators
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
    ///////////////////////////////////////////////////////////////////////////////////////
    // localInnerProduct Scalar x Scalar -> Scalar
    // localInnerProduct Vector x Vector -> Scalar
    // localInnerProduct Matrix x Matrix -> Scalar
    ///////////////////////////////////////////////////////////////////////////////////////
    template<class l,class r,int N> inline
    auto localInnerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(localInnerProduct(lhs._internal[0],rhs._internal[0]))>
    {
        typedef decltype(localInnerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
        iScalar<ret_t> ret=zero;
        for(int c1=0;c1<N;c1++){
            ret._internal += localInnerProduct(lhs._internal[c1],rhs._internal[c1]);
        }
        return ret;
    }
    template<class l,class r,int N> inline
    auto localInnerProduct (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(localInnerProduct(lhs._internal[0][0],rhs._internal[0][0]))>
    {
        typedef decltype(localInnerProduct(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
        iScalar<ret_t> ret=zero;
        for(int c1=0;c1<N;c1++){
        for(int c2=0;c2<N;c2++){
            ret._internal += localInnerProduct(lhs._internal[c1][c2],rhs._internal[c1][c2]);
        }}
        return ret;
    }
    template<class l,class r> inline
    auto localInnerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(localInnerProduct(lhs._internal,rhs._internal))>
    {
        typedef decltype(localInnerProduct(lhs._internal,rhs._internal)) ret_t;
        iScalar<ret_t> ret;
        ret._internal = localInnerProduct(lhs._internal,rhs._internal);
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
// Can only take the real/imag part of scalar objects, since
// lattice objects of different complexity are non-conformable.
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

    /////////////////////////////////
    // Trace of scalar and matrix
    /////////////////////////////////

inline Complex trace( const Complex &arg){
    return arg;
}
//inline vComplex trace(const vComplex &arg){
//    return arg;
//}
template<class vtype,int N>
inline auto trace(const iMatrix<vtype,N> &arg) -> iScalar<decltype(trace(arg._internal[0][0]))>
{
    iScalar<decltype( trace(arg._internal[0][0] )) > ret;
    ZeroIt(ret._internal);
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
};
/////////////////////////////////////////////////////////////////////////
// Generic routine to promote object<complex> -> object<vcomplex>
// Supports the array reordering transformation that gives me SIMD utilisation
/////////////////////////////////////////////////////////////////////////
/*
template<template<class> class object>
inline object<vComplex> splat(object<Complex >s){
    object<vComplex> ret;
    vComplex * v_ptr = (vComplex *)& ret;
    Complex * s_ptr = (Complex *) &s;
    for(int i=0;i<sizeof(ret);i+=sizeof(vComplex)){
        vsplat(*(v_ptr++),*(s_ptr++));
    }
    return ret;
}
*/
    
#endif
