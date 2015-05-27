//---------------------------------------------------------------------------
/*! @file Grid_vector_types.h
  @brief Defines templated class Grid_simd to deal with inner vector types
*/
// Time-stamp: <2015-05-26 14:08:13 neo>
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
#include "Grid_knc.h"
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
  template <typename Condition, typename ReturnType> using EnableIf   =    Invoke<std::enable_if<Condition::value, ReturnType>>;
  template <typename Condition, typename ReturnType> using NotEnableIf=    Invoke<std::enable_if<!Condition::value, ReturnType>>;
  
  ////////////////////////////////////////////////////////
  // Check for complexity with type traits
  template <typename T>     struct is_complex : std::false_type {};
  template < typename T >   struct is_complex< std::complex<T> >: std::true_type {};
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


  // Move to the simd files
//////////////////////////////////////////////////////////
// Permute
// Permute 0 every ABCDEFGH -> BA DC FE HG
// Permute 1 every ABCDEFGH -> CD AB GH EF
// Permute 2 every ABCDEFGH -> EFGH ABCD
// Permute 3 possible on longer iVector lengths (512bit = 8 double = 16 single)
// Permute 4 possible on half precision @512bit vectors.
//////////////////////////////////////////////////////////
template<class vsimd>
inline void Gpermute(vsimd &y,const vsimd &b,int perm){
	union { 
	  SIMD_Ftype f;
	  decltype(vsimd::v) v;
	} conv;
	conv.v = b.v;
      switch (perm){
#if defined(AVX1)||defined(AVX2)
      // 8x32 bits=>3 permutes
      case 2: 
	conv.f = _mm256_shuffle_ps(conv.f,conv.f,_MM_SHUFFLE(2,3,0,1)); 
	break;
      case 1: conv.f = _mm256_shuffle_ps(conv.f,conv.f,_MM_SHUFFLE(1,0,3,2)); break;
      case 0: conv.f = _mm256_permute2f128_ps(conv.f,conv.f,0x01); break;
#endif
#ifdef SSE4
      case 1: conv.f = _mm_shuffle_ps(conv.f,conv.f,_MM_SHUFFLE(2,3,0,1)); break;
      case 0: conv.f = _mm_shuffle_ps(conv.f,conv.f,_MM_SHUFFLE(1,0,3,2));break;
#endif
#ifdef AVX512
	// 16 floats=> permutes
        // Permute 0 every abcd efgh ijkl mnop -> badc fehg jilk nmpo 
        // Permute 1 every abcd efgh ijkl mnop -> cdab ghef jkij opmn 
        // Permute 2 every abcd efgh ijkl mnop -> efgh abcd mnop ijkl
        // Permute 3 every abcd efgh ijkl mnop -> ijkl mnop abcd efgh
      case 3: conv.f = _mm512_swizzle_ps(conv.f,_MM_SWIZ_REG_CDAB); break;
      case 2: conv.f = _mm512_swizzle_ps(conv.f,_MM_SWIZ_REG_BADC); break;
      case 1: conv.f = _mm512_permute4f128_ps(conv.f,(_MM_PERM_ENUM)_MM_SHUFFLE(2,3,0,1)); break;
      case 0: conv.f = _mm512_permute4f128_ps(conv.f,(_MM_PERM_ENUM)_MM_SHUFFLE(1,0,3,2)); break;
#endif
#ifdef QPX
#error not implemented
#endif
      default: assert(0); break;
      }
      y.v=conv.v;

 };

///////////////////////////////////////



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
    
    // Constructors
    Grid_simd & operator = ( Zero & z){
      vzero(*this);
      return (*this);
    }
    
    Grid_simd& operator=(const Grid_simd&& rhs){v=rhs.v;return *this;};
    Grid_simd& operator=(const Grid_simd& rhs){v=rhs.v;return *this;}; //faster than not declaring it and leaving to the compiler
    Grid_simd()=default; 
    Grid_simd(const Grid_simd& rhs):v(rhs.v){};    //compiles in movaps
    Grid_simd(const Grid_simd&& rhs):v(rhs.v){};  
  
    //Enable if complex type
    template < class S = Scalar_type > 
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



    //not for integer types... 
    template <  class S = Scalar_type, NotEnableIf<std::is_integral < S >, int> = 0 > 
    friend inline Grid_simd adj(const Grid_simd &in){ return conjugate(in); }
        
    ///////////////////////////////////////////////
    // Initialise to 1,0,i for the correct types
    ///////////////////////////////////////////////
    // For complex types
    template <  class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 
      friend inline void vone(Grid_simd &ret)      { vsplat(ret,1.0,0.0); }
    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 
      friend inline void vzero(Grid_simd &ret)     { vsplat(ret,0.0,0.0); }// use xor?
    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 
      friend inline void vcomplex_i(Grid_simd &ret){ vsplat(ret,0.0,1.0);} 

    // if not complex overload here 
    template <  class S = Scalar_type, EnableIf<std::is_floating_point < S >,int> = 0 > 
      friend inline void vone(Grid_simd &ret)      { vsplat(ret,1.0); }
    template <  class S = Scalar_type, EnableIf<std::is_floating_point < S >,int> = 0 > 
      friend inline void vzero(Grid_simd &ret)     { vsplat(ret,0.0); }
    

   
    // For integral types
    template <  class S = Scalar_type, EnableIf<std::is_integral < S >, int> = 0 > 
      friend inline void vone(Grid_simd &ret)      { vsplat(ret,1); }
    template <  class S = Scalar_type, EnableIf<std::is_integral < S >, int> = 0 > 
      friend inline void vzero(Grid_simd &ret)      { vsplat(ret,0); }
    template <  class S = Scalar_type, EnableIf<std::is_integral < S >, int> = 0 > 
      friend inline void vtrue (Grid_simd &ret){vsplat(ret,0xFFFFFFFF);}
    template <  class S = Scalar_type, EnableIf<std::is_integral < S >, int> = 0 > 
      friend inline void vfalse(Grid_simd &ret){vsplat(ret,0);}
   
    ////////////////////////////////////
    // Arithmetic operator overloads +,-,*
    ////////////////////////////////////
    friend inline Grid_simd operator + (Grid_simd a, Grid_simd b)
    {
      Grid_simd ret;
      ret.v = binary<Vector_type>(a.v, b.v, SumSIMD());
      return ret;
    };
        
    friend inline Grid_simd operator - (Grid_simd a, Grid_simd b)
    {
      Grid_simd ret;
      ret.v = binary<Vector_type>(a.v, b.v, SubSIMD());
      return ret;
    };
        
    // Distinguish between complex types and others
    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 >
      friend inline Grid_simd operator * (Grid_simd a, Grid_simd b)
      {
	Grid_simd ret;
	ret.v = binary<Vector_type>(a.v,b.v, MultComplexSIMD());
	return ret;
      };

    // Real/Integer types
    template <  class S = Scalar_type, NotEnableIf<is_complex < S >, int> = 0 > 
    friend inline Grid_simd operator * (Grid_simd a, Grid_simd b)
      {
	Grid_simd ret;
	ret.v = binary<Vector_type>(a.v,b.v, MultSIMD());
	return ret;
      };

    ////////////////////////////////////////////////////////////////////////
    // FIXME:  gonna remove these load/store, get, set, prefetch
    ////////////////////////////////////////////////////////////////////////
    friend inline void vset(Grid_simd &ret, Scalar_type *a){
      ret.v = unary<Vector_type>(a, VsetSIMD());
    }
        
    ///////////////////////
    // Splat
    ///////////////////////
    // overload if complex
    template < class S = Scalar_type > 
    friend inline void vsplat(Grid_simd &ret, EnableIf<is_complex < S >, S> c){
      Real a = real(c);
      Real b = imag(c);
      vsplat(ret,a,b);
    }

    // this is only for the complex version
    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 
    friend inline void vsplat(Grid_simd &ret,Real a, Real b){
      ret.v = binary<Vector_type>(a, b, VsplatSIMD());
    }    

    //if real fill with a, if complex fill with a in the real part (first function above)
    friend inline void vsplat(Grid_simd &ret,Real a){
      ret.v = unary<Vector_type>(a, VsplatSIMD());
    }    

    ///////////////////////
    // Vstore
    ///////////////////////
    friend inline void vstore(const Grid_simd &ret, Scalar_type *a){
      binary<void>(ret.v, (Real*)a, VstoreSIMD());
    }

    ///////////////////////
    // Vstream
    ///////////////////////
    template <  class S = Scalar_type, NotEnableIf<std::is_integral < S >, int> = 0 > 
    friend inline void vstream(Grid_simd &out,const Grid_simd &in){
      binary<void>((Real*)&out.v, in.v, VstreamSIMD());
    }

    template <  class S = Scalar_type, EnableIf<std::is_integral < S >, int> = 0 > 
      friend inline void vstream(Grid_simd &out,const Grid_simd &in){
      out=in;
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
    // Conjugate
    ///////////////////////
    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 							     
    friend inline Grid_simd  conjugate(const Grid_simd  &in){
      Grid_simd  ret ; 
      ret.v = unary<Vector_type>(in.v, ConjSIMD());
      return ret;
    }
    template < class S = Scalar_type, NotEnableIf<is_complex < S >, int> = 0 > 
    friend inline Grid_simd  conjugate(const Grid_simd  &in){
      return in; // for real objects
    }


    ///////////////////////
    // timesMinusI
    ///////////////////////
    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 
    friend inline void timesMinusI( Grid_simd &ret,const Grid_simd &in){
      ret.v = binary<Vector_type>(in.v, ret.v, TimesMinusISIMD());
    }

    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 
    friend inline Grid_simd timesMinusI(const Grid_simd &in){
      Grid_simd ret; 
      timesMinusI(ret,in);
      return ret;
    }

    template < class S = Scalar_type, NotEnableIf<is_complex < S >, int> = 0 > 
    friend inline Grid_simd timesMinusI(const Grid_simd &in){
      return in;
    }


    ///////////////////////
    // timesI
    ///////////////////////
    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 
    friend inline void timesI(Grid_simd &ret,const Grid_simd &in){
      ret.v =   binary<Vector_type>(in.v, ret.v, TimesISIMD());     
    }
        
    template < class S = Scalar_type, EnableIf<is_complex < S >, int> = 0 > 
    friend inline Grid_simd timesI(const Grid_simd &in){
      Grid_simd ret; 
      timesI(ret,in);
      return ret;
    }

    template < class S = Scalar_type, NotEnableIf<is_complex < S >, int> = 0 > 
    friend inline Grid_simd timesI(const Grid_simd &in){
      return in;
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



  template<class scalar_type, class vector_type > 
    inline Grid_simd< scalar_type, vector_type>  innerProduct(const Grid_simd< scalar_type, vector_type> & l, const Grid_simd< scalar_type, vector_type> & r) 
  {
    return conjugate(l)*r; 
  }

  template<class scalar_type, class vector_type >
    inline void zeroit(Grid_simd< scalar_type, vector_type> &z){ vzero(z);}


  template<class scalar_type, class vector_type >
  inline Grid_simd< scalar_type, vector_type> outerProduct(const Grid_simd< scalar_type, vector_type> &l, const Grid_simd< scalar_type, vector_type>& r)
  {
    return l*r;
  }


  template<class scalar_type, class vector_type >
  inline Grid_simd< scalar_type, vector_type> trace(const Grid_simd< scalar_type, vector_type> &arg){
    return arg;
  }


  // Define available types

  typedef Grid_simd< float                 , SIMD_Ftype > vRealF;
  typedef Grid_simd< double                , SIMD_Dtype > vRealD;
  typedef Grid_simd< std::complex< float > , SIMD_Ftype > vComplexF;
  typedef Grid_simd< std::complex< double >, SIMD_Dtype > vComplexD;
  typedef Grid_simd< Integer               , SIMD_Itype > vInteger;







}

#endif
