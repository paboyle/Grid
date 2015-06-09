//----------------------------------------------------------------------
/*! @file Grid_sse4.h
  @brief Optimization libraries for NEON (ARM) instructions set ARMv7

  Experimental - Using intrinsics - DEVELOPING! 
*/
// Time-stamp: <2015-06-09 15:25:40 neo>
//----------------------------------------------------------------------

#include <arm_neon.h>

namespace Optimization {

  template<class vtype>
  union uconv {
    float32x4_t f;
    vtype v;
  };

  union u128f {
    float32x4_t v;
    float f[4];
  };
  union u128d {
    float32x4_t v;
    float f[4];
  };
  
  struct Vsplat{
    //Complex float
    inline float32x4_t operator()(float a, float b){
      float32x4_t foo;
      return foo;
    }
    // Real float
    inline float32x4_t operator()(float a){
      float32x4_t foo;
      return foo;
    }
    //Complex double
    inline float32x4_t operator()(double a, double b){
      float32x4_t foo;
      return foo;
    }
    //Real double
    inline float32x4_t operator()(double a){
      float32x4_t foo;
      return foo;
    }
    //Integer
    inline uint32x4_t operator()(Integer a){
      uint32x4_t foo;
      return foo;
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(float32x4_t a, float* F){
      
    }
    //Double
    inline void operator()(float32x4_t a, double* D){
      
    }
    //Integer
    inline void operator()(uint32x4_t a, Integer* I){
     
    }

  };

  struct Vstream{
    //Float
    inline void operator()(float * a, float32x4_t b){
    
    }
    //Double
    inline void operator()(double * a, float32x4_t b){
  
    }


  };

  struct Vset{
    // Complex float 
    inline float32x4_t operator()(Grid::ComplexF *a){
      float32x4_t foo;
      return foo;
    }
    // Complex double 
    inline float32x4_t operator()(Grid::ComplexD *a){
      float32x4_t foo;
      return foo;
    }
    // Real float 
    inline float32x4_t operator()(float *a){
      float32x4_t foo;
      return foo;
    }
    // Real double
    inline float32x4_t operator()(double *a){
      float32x4_t foo;
      return foo;
    }
    // Integer
    inline uint32x4_t operator()(Integer *a){
      uint32x4_t foo;
      return foo;
    }


  };

  template <typename Out_type, typename In_type>
  struct Reduce{
    //Need templated class to overload output type
    //General form must generate error if compiled
    inline Out_type operator()(In_type in){
      printf("Error, using wrong Reduce function\n");
      exit(1);
      return 0;
    }
  };

  /////////////////////////////////////////////////////
  // Arithmetic operations
  /////////////////////////////////////////////////////
  struct Sum{
    //Complex/Real float
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      float32x4_t foo;
      return foo;
    }
    //Complex/Real double
    //inline float32x4_t operator()(float32x4_t a, float32x4_t b){
    //  float32x4_t foo;
    //  return foo;
    //}
    //Integer
    inline uint32x4_t operator()(uint32x4_t a, uint32x4_t b){
      uint32x4_t foo;
      return foo;
    }
  };

  struct Sub{
    //Complex/Real float
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      float32x4_t foo;
      return foo;
    }
    //Complex/Real double
    //inline float32x4_t operator()(float32x4_t a, float32x4_t b){
    //  float32x4_t foo;
    //  return foo;
    //}
    //Integer
    inline uint32x4_t operator()(uint32x4_t a, uint32x4_t b){
      uint32x4_t foo;
      return foo;
    }
  };

  struct MultComplex{
    // Complex float
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      float32x4_t foo;
      return foo;
    }
    // Complex double
    //inline float32x4_t operator()(float32x4_t a, float32x4_t b){
    //  float32x4_t foo;
    //  return foo;
    //}
  };

  struct Mult{
    // Real float
    inline float32x4_t operator()(float32x4_t a, float32x4_t b){
      return a;
    }
    // Real double
    //inline float32x4_t operator()(float32x4_t a, float32x4_t b){
    //  return 0;
    //}
    // Integer
    inline uint32x4_t operator()(uint32x4_t a, uint32x4_t b){
      return a;
    }
  };

  struct Conj{
    // Complex single
    inline float32x4_t operator()(float32x4_t in){
      return in;
    }
    // Complex double
    //inline float32x4_t operator()(float32x4_t in){
    // return 0;
    //}
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
      return in;
    }
    //Complex double
    //inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
    //  return in;
    //}


  };

  struct TimesI{
    //Complex single
    inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
      return in;
    }
    //Complex double
    //inline float32x4_t operator()(float32x4_t in, float32x4_t ret){
    //  return 0;
    //}
  };

  //////////////////////////////////////////////
  // Some Template specialization
  template < typename vtype > 
    void permute(vtype &a, vtype b, int perm) {

  }; 

  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, float32x4_t>::operator()(float32x4_t in){
    return 0;
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, float32x4_t>::operator()(float32x4_t in){
    return 0;
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, float32x4_t>::operator()(float32x4_t in){
    return 0;
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, float32x4_t>::operator()(float32x4_t in){
    return 0;
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, uint32x4_t>::operator()(uint32x4_t in){
    // FIXME unimplemented
   printf("Reduce : Missing integer implementation -> FIX\n");
    assert(0);
  }
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 
namespace Grid {

  typedef float32x4_t  SIMD_Ftype; // Single precision type
  typedef float32x4_t  SIMD_Dtype; // Double precision type - no double on ARMv7
  typedef uint32x4_t   SIMD_Itype; // Integer type

  inline void v_prefetch0(int size, const char *ptr){};  // prefetch utilities
  inline void prefetch_HINT_T0(const char *ptr){};


  // Gpermute function
  template < typename VectorSIMD > 
    inline void Gpermute(VectorSIMD &y,const VectorSIMD &b, int perm ) {
    Optimization::permute(y.v,b.v,perm);
  }


  // Function name aliases
  typedef Optimization::Vsplat   VsplatSIMD;
  typedef Optimization::Vstore   VstoreSIMD;
  typedef Optimization::Vset     VsetSIMD;
  typedef Optimization::Vstream  VstreamSIMD;
  template <typename S, typename T> using ReduceSIMD = Optimization::Reduce<S,T>;

 


  // Arithmetic operations
  typedef Optimization::Sum         SumSIMD;
  typedef Optimization::Sub         SubSIMD;
  typedef Optimization::Mult        MultSIMD;
  typedef Optimization::MultComplex MultComplexSIMD;
  typedef Optimization::Conj        ConjSIMD;
  typedef Optimization::TimesMinusI TimesMinusISIMD;
  typedef Optimization::TimesI      TimesISIMD;

}
