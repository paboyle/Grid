//---------------------------------------------------------------------------
/*! @file Grid_vector_types.h
  @brief Defines templated class Grid_simd to deal with inner vector types
*/
// Time-stamp: <2015-05-27 12:04:06 neo>
//---------------------------------------------------------------------------
#ifndef GRID_VECTOR_TYPES
#define GRID_VECTOR_TYPES

#ifdef SSE4
#include "Grid_sse4.h"
#endif
#if defined (AVX1)|| defined (AVX2)
#include "Grid_avx.h"
#endif
#if defined AVX512
#include "Grid_avx512.h"
#endif
#if defined QPX
#include "Grid_qpx.h"
#endif

namespace Grid {

  // To take the floating point type of real/complex type
  template <typename T> struct RealPart {
    typedef T type;
  };
  template <typename T> struct RealPart< std::complex<T> >{
    typedef T type;
  };

  // type alias used to simplify the syntax of std::enable_if
  template <typename T> using Invoke                                  =  typename T::type;
  template <typename Condition, typename ReturnType> using EnableIf   =  Invoke<std::enable_if<Condition::value, ReturnType> >;
  template <typename Condition, typename ReturnType> using NotEnableIf=  Invoke<std::enable_if<!Condition::value, ReturnType> >;


  ////////////////////////////////////////////////////////
  // Check for complexity with type traits
  template <typename T>   struct is_complex : public std::false_type {};
  template <> struct is_complex<std::complex<double> >: public std::true_type {};
  template <> struct is_complex<std::complex<float> > : public std::true_type {};

  template <typename T> using IfReal    = Invoke<std::enable_if<std::is_floating_point<T>::value,int> > ;
  template <typename T> using IfComplex = Invoke<std::enable_if<is_complex<T>::value,int> > ;
  template <typename T> using IfInteger = Invoke<std::enable_if<std::is_integral<T>::value,int> > ;

  template <typename T> using IfNotReal    = Invoke<std::enable_if<!std::is_floating_point<T>::value,int> > ;
  template <typename T> using IfNotComplex = Invoke<std::enable_if<!is_complex<T>::value,int> > ;
  template <typename T> using IfNotInteger = Invoke<std::enable_if<!std::is_integral<T>::value,int> > ;

  ////////////////////////////////////////////////////////
  // Define the operation templates functors
  // general forms to allow for vsplat syntax
  // need explicit declaration of types when used since
  // clang cannot automatically determine the output type sometimes
  // use decltype?
  template < class Out, class Input1, class Input2, class Operation > 
    Out binary(Input1 src_1, Input2 src_2, Operation op){
    return op(src_1, src_2);
  } 

  template < class Out, class Input, class Operation > 
    Out unary(Input src, Operation op){
    return op(src);
  } 
  ///////////////////////////////////////////////



  /*
    @brief Grid_simd class for the SIMD vector type operations
   */
  template < class Scalar_type, class Vector_type > 
    class Grid_simd {
    
  public:
    typedef typename RealPart < Scalar_type >::type Real; 
    typedef Vector_type     vector_type;
    typedef Scalar_type     scalar_type;
   
    Vector_type v;
    
    static inline int Nsimd(void) { return sizeof(Vector_type)/sizeof(Scalar_type);}

    
    Grid_simd& operator=(const Grid_simd&& rhs){v=rhs.v;return *this;};
    Grid_simd& operator=(const Grid_simd& rhs){v=rhs.v;return *this;}; //faster than not declaring it and leaving to the compiler
    Grid_simd()=default; 
    Grid_simd(const Grid_simd& rhs) :v(rhs.v){};    //compiles in movaps
    Grid_simd(const Grid_simd&& rhs):v(rhs.v){};  

    /////////////////////////////
    // Constructors
    /////////////////////////////
    Grid_simd & operator = ( Zero & z){
      vzero(*this);
      return (*this);
    }
  
    //Enable if complex type
    template < typename S = Scalar_type > 
    Grid_simd(const typename std::enable_if< is_complex < S >::value, S>::type a){
      vsplat(*this,a);
    };

    Grid_simd(const Real a){
      vsplat(*this,Scalar_type(a));
    };
       
    ///////////////////////////////////////////////
    // mac, mult, sub, add, adj
    ///////////////////////////////////////////////
    friend inline void mac (Grid_simd * __restrict__ y,const Grid_simd * __restrict__ a,const Grid_simd *__restrict__ x){ *y = (*a)*(*x)+(*y); };
    friend inline void mult(Grid_simd * __restrict__ y,const Grid_simd * __restrict__ l,const Grid_simd *__restrict__ r){ *y = (*l) * (*r); }
    friend inline void sub (Grid_simd * __restrict__ y,const Grid_simd * __restrict__ l,const Grid_simd *__restrict__ r){ *y = (*l) - (*r); }
    friend inline void add (Grid_simd * __restrict__ y,const Grid_simd * __restrict__ l,const Grid_simd *__restrict__ r){ *y = (*l) + (*r); }

    friend inline void mac (Grid_simd *__restrict__ y,const Scalar_type *__restrict__ a,const Grid_simd   *__restrict__ x){ *y = (*a)*(*x)+(*y); };
    friend inline void mult(Grid_simd *__restrict__ y,const Scalar_type *__restrict__ l,const Grid_simd   *__restrict__ r){ *y = (*l) * (*r); }
    friend inline void sub (Grid_simd *__restrict__ y,const Scalar_type *__restrict__ l,const Grid_simd   *__restrict__ r){ *y = (*l) - (*r); }
    friend inline void add (Grid_simd *__restrict__ y,const Scalar_type *__restrict__ l,const Grid_simd   *__restrict__ r){ *y = (*l) + (*r); }

    friend inline void mac (Grid_simd *__restrict__ y,const Grid_simd   *__restrict__ a,const Scalar_type *__restrict__ x){ *y = (*a)*(*x)+(*y); };
    friend inline void mult(Grid_simd *__restrict__ y,const Grid_simd   *__restrict__ l,const Scalar_type *__restrict__ r){ *y = (*l) * (*r); }
    friend inline void sub (Grid_simd *__restrict__ y,const Grid_simd   *__restrict__ l,const Scalar_type *__restrict__ r){ *y = (*l) - (*r); }
    friend inline void add (Grid_simd *__restrict__ y,const Grid_simd   *__restrict__ l,const Scalar_type *__restrict__ r){ *y = (*l) + (*r); }

    ////////////////////////////////////////////////////////////////////////
    // FIXME:  gonna remove these load/store, get, set, prefetch
    ////////////////////////////////////////////////////////////////////////
    friend inline void vset(Grid_simd &ret, Scalar_type *a){
      ret.v = unary<Vector_type>(a, VsetSIMD());
    }
        
    ///////////////////////
    // Vstore
    ///////////////////////
    friend inline void vstore(const Grid_simd &ret, Scalar_type *a){
      binary<void>(ret.v, (Real*)a, VstoreSIMD());
    }

    ///////////////////////
    // Vprefetch
    ///////////////////////
    friend inline void vprefetch(const Grid_simd &v)
    {
      _mm_prefetch((const char*)&v.v,_MM_HINT_T0);
    }

    ///////////////////////
    // Reduce
    ///////////////////////
    friend inline Scalar_type Reduce(const Grid_simd & in)
    {
      return unary<Scalar_type>(in.v, ReduceSIMD<Scalar_type, Vector_type>());
    }

    ////////////////////////////
    // opreator scalar * simd
    ////////////////////////////
    friend inline Grid_simd operator * (const Scalar_type &a, Grid_simd b){
      Grid_simd va;
      vsplat(va,a);
      return va*b;
    }
    friend inline Grid_simd operator * (Grid_simd b,const Scalar_type &a){
      return a*b;
    }

    ///////////////////////
    // Unary negation
    ///////////////////////
    friend inline Grid_simd operator -(const Grid_simd &r) {
      Grid_simd ret;
      vzero(ret);
      ret = ret - r;
      return ret;
    }
    // *=,+=,-= operators
    inline Grid_simd &operator *=(const Grid_simd &r) {
      *this = (*this)*r;
      return *this;
      // return (*this)*r; ?
    }
    inline Grid_simd &operator +=(const Grid_simd &r) {
      *this = *this+r;
      return *this;
    }
    inline Grid_simd &operator -=(const Grid_simd &r) {
      *this = *this-r;
      return *this;
    }

    ////////////////////////////////////////////////////////////////////
    // General permute; assumes vector length is same across 
    // all subtypes; may not be a good assumption, but could
    // add the vector width as a template param for BG/Q for example
    ////////////////////////////////////////////////////////////////////
    friend inline void permute(Grid_simd &y,Grid_simd b,int perm)
    {
      Gpermute<Grid_simd>(y,b,perm);
    }
    
  };// end of Grid_simd class definition 


  ///////////////////////
  // Splat
  ///////////////////////
  
  // this is only for the complex version
  template <class S, class V, IfComplex<S> =0, class ABtype> 
    inline void vsplat(Grid_simd<S,V> &ret,ABtype a, ABtype b){
    ret.v = binary<V>(a, b, VsplatSIMD());
  }    

  // overload if complex
  template <class S,class V> inline void vsplat(Grid_simd<S,V> &ret, EnableIf<is_complex < S >, S> c) {
    Real a = real(c);
    Real b = imag(c);
    vsplat(ret,a,b);
  }

  //if real fill with a, if complex fill with a in the real part (first function above)
  template <class S,class V>
    inline void vsplat(Grid_simd<S,V> &ret,NotEnableIf<is_complex< S>,S> a){
    ret.v = unary<V>(a, VsplatSIMD());
  }    
  //////////////////////////

  ///////////////////////////////////////////////
  // Initialise to 1,0,i for the correct types
  ///////////////////////////////////////////////
  // For complex types
  template <class S,class V, IfComplex<S> = 0 > inline void vone(Grid_simd<S,V>  &ret)     { vsplat(ret,S(1.0,0.0)); }
  template <class S,class V, IfComplex<S> = 0 > inline void vzero(Grid_simd<S,V> &ret)     { vsplat(ret,S(0.0,0.0)); }// use xor?
  template <class S,class V, IfComplex<S> = 0 > inline void vcomplex_i(Grid_simd<S,V> &ret){ vsplat(ret,S(0.0,1.0));} 

  // if not complex overload here 
  template <class S,class V, IfReal<S> = 0 > inline void vone (Grid_simd<S,V> &ret){ vsplat(ret,1.0); }
  template <class S,class V, IfReal<S> = 0 > inline void vzero(Grid_simd<S,V> &ret)     { vsplat(ret,0.0); }
   
  // For integral types
  template <class S,class V,IfInteger<S> = 0 > inline void vone(Grid_simd<S,V> &ret)  {vsplat(ret,1); }
  template <class S,class V,IfInteger<S> = 0 > inline void vzero(Grid_simd<S,V> &ret) {vsplat(ret,0); }
  template <class S,class V,IfInteger<S> = 0 > inline void vtrue (Grid_simd<S,V> &ret){vsplat(ret,0xFFFFFFFF);}
  template <class S,class V,IfInteger<S> = 0 > inline void vfalse(Grid_simd<S,V> &ret){vsplat(ret,0);}

  template<class S,class V> inline void zeroit(Grid_simd<S,V> &z){ vzero(z);}

  ///////////////////////
  // Vstream
  ///////////////////////
  template <class S,class V, IfNotInteger<S> = 0 > 
    inline void vstream(Grid_simd<S,V> &out,const Grid_simd<S,V> &in){
      binary<void>((Real*)&out.v, in.v, VstreamSIMD());
    }

  template <class S,class V, IfInteger<S> = 0 > 
    inline void vstream(Grid_simd<S,V> &out,const Grid_simd<S,V> &in){
    out=in;
  }

  ////////////////////////////////////
  // Arithmetic operator overloads +,-,*
  ////////////////////////////////////
  template<class S,class V> inline Grid_simd<S,V> operator + (Grid_simd<S,V> a, Grid_simd<S,V> b) {
    Grid_simd<S,V> ret;
    ret.v = binary<V>(a.v, b.v, SumSIMD());
    return ret;
  };
        
  template<class S,class V> inline Grid_simd<S,V> operator - (Grid_simd<S,V> a, Grid_simd<S,V> b) {
    Grid_simd<S,V> ret;
    ret.v = binary<V>(a.v, b.v, SubSIMD());
    return ret;
  };
        
  // Distinguish between complex types and others
  template<class S,class V, IfComplex<S> = 0 > inline Grid_simd<S,V> operator * (Grid_simd<S,V> a, Grid_simd<S,V> b) {
    Grid_simd<S,V> ret;
    ret.v = binary<V>(a.v,b.v, MultComplexSIMD());
    return ret;
  };

    // Real/Integer types
  template<class S,class V, IfNotComplex<S> = 0 > inline Grid_simd<S,V> operator * (Grid_simd<S,V> a, Grid_simd<S,V> b) {
    Grid_simd<S,V> ret;
    ret.v = binary<V>(a.v,b.v, MultSIMD());
    return ret;
  };
  

    ///////////////////////
    // Conjugate
    ///////////////////////
  template <class S,class V, IfComplex<S> = 0 > 
    inline Grid_simd<S,V> conjugate(const Grid_simd<S,V>  &in){
    Grid_simd<S,V>  ret ; 
    ret.v = unary<V>(in.v, ConjSIMD());
    return ret;
  }
  template <class S,class V, IfNotComplex<S> = 0 > inline Grid_simd<S,V> conjugate(const Grid_simd<S,V>  &in){
    return in; // for real objects
  }

  //Suppress adj for integer types... // odd; why conjugate above but not adj??
  template < class S, class V, IfNotInteger<S> = 0 > 
    inline Grid_simd<S,V> adj(const Grid_simd<S,V> &in){ return conjugate(in); }
  
  ///////////////////////
  // timesMinusI
  ///////////////////////
  template<class S,class V,IfComplex<S> = 0 > 
    inline void timesMinusI( Grid_simd<S,V> &ret,const Grid_simd<S,V> &in){
    ret.v = binary<V>(in.v, ret.v, TimesMinusISIMD());
  }

  template<class S,class V,IfComplex<S> = 0 > 
    inline Grid_simd<S,V> timesMinusI(const Grid_simd<S,V> &in){
    Grid_simd<S,V> ret; 
    timesMinusI(ret,in);
    return ret;
  }

  template<class S,class V,IfNotComplex<S> = 0 > 
    inline Grid_simd<S,V> timesMinusI(const Grid_simd<S,V> &in){
    return in;
  }

    ///////////////////////
    // timesI
    ///////////////////////
  template<class S,class V,IfComplex<S> = 0 > 
    inline void timesI(Grid_simd<S,V> &ret,const Grid_simd<S,V> &in){
    ret.v =   binary<V>(in.v, ret.v, TimesISIMD());     
  }
        
  template<class S,class V,IfComplex<S> = 0 > 
    inline Grid_simd<S,V> timesI(const Grid_simd<S,V> &in){
    Grid_simd<S,V> ret; 
    timesI(ret,in);
    return ret;
  }

  template<class S,class V,IfNotComplex<S> = 0 > 
    inline Grid_simd<S,V> timesI(const Grid_simd<S,V> &in){
    return in;
  }


  /////////////////////
  // Inner, outer
  /////////////////////

  template<class S, class V > 
    inline Grid_simd< S, V>  innerProduct(const Grid_simd< S, V> & l, const Grid_simd< S, V> & r) 
  {
    return conjugate(l)*r; 
  }

  template<class S, class V >
  inline Grid_simd< S, V> outerProduct(const Grid_simd< S, V> &l, const Grid_simd< S, V> & r)
  {
    return l*r;
  }

  template<class S, class V >
  inline Grid_simd< S, V> trace(const Grid_simd< S, V> &arg){
    return arg;
  }

  ///////////////////////////////
  // Define available types
  ///////////////////////////////
  typedef Grid_simd< float                 , SIMD_Ftype > vRealF;
  typedef Grid_simd< double                , SIMD_Dtype > vRealD;
  typedef Grid_simd< std::complex< float > , SIMD_Ftype > vComplexF;
  typedef Grid_simd< std::complex< double >, SIMD_Dtype > vComplexD;
  typedef Grid_simd< Integer               , SIMD_Itype > vInteger;
}

#endif
