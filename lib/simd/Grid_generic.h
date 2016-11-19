    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_generic.h

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */

static_assert(GEN_SIMD_WIDTH % 16u == 0, "SIMD vector size is not an integer multiple of 16 bytes");

//#define VECTOR_LOOPS

// playing with compiler pragmas
#ifdef VECTOR_LOOPS
#ifdef __clang__
#define VECTOR_FOR(i, w, inc)\
_Pragma("clang loop unroll(full) vectorize(enable) interleave(enable) vectorize_width(w)")\
for (unsigned int i = 0; i < w; i += inc)
#elif defined __INTEL_COMPILER
#define VECTOR_FOR(i, w, inc)\
_Pragma("simd vectorlength(w*8)")\
for (unsigned int i = 0; i < w; i += inc)
#else
#define VECTOR_FOR(i, w, inc)\
for (unsigned int i = 0; i < w; i += inc)
#endif
#else
#define VECTOR_FOR(i, w, inc)\
for (unsigned int i = 0; i < w; i += inc)
#endif

namespace Grid {
namespace Optimization {

  // type traits giving the number of elements for each vector type
  template <typename T> struct W;
  template <> struct W<double> {
    constexpr static unsigned int c = GEN_SIMD_WIDTH/16u;
    constexpr static unsigned int r = GEN_SIMD_WIDTH/8u;
  };
  template <> struct W<float> {
    constexpr static unsigned int c = GEN_SIMD_WIDTH/8u;
    constexpr static unsigned int r = GEN_SIMD_WIDTH/4u;
  };
  
  // SIMD vector types
  template <typename T>
  struct vec {
    alignas(GEN_SIMD_WIDTH) T v[W<T>::r];
  };
  
  typedef vec<float>   vecf;
  typedef vec<double>  vecd;
  
  struct Vsplat{
    // Complex
    template <typename T>
    inline vec<T> operator()(T a, T b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::r, 2)
      {
        out.v[i]   = a;
        out.v[i+1] = b;
      }

      return out;
    }
    
    // Real
    template <typename T>
    inline vec<T> operator()(T a){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::r, 1)
      {
        out.v[i] = a;
      }
      
      return out;
    }
    
    // Integer
    inline int operator()(Integer a){
      return a;
    }
  };

  struct Vstore{
    // Real
    template <typename T>
    inline void operator()(vec<T> a, T *D){
      *((vec<T> *)D) = a;
    }
    //Integer
    inline void operator()(int a, Integer *I){
      *I = a;
    }

  };

  struct Vstream{
    // Real
    template <typename T>
    inline void operator()(T * a, vec<T> b){
      *((vec<T> *)a) = b;
    }
  };

  struct Vset{
    // Complex
    template <typename T>
    inline vec<T> operator()(std::complex<T> *a){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::c, 1)
      {
        out.v[2*i]   = a[i].real();
        out.v[2*i+1] = a[i].imag();
      }
      
      return out;
    }
    
    // Real
    template <typename T>
    inline vec<T> operator()(T *a){
      vec<T> out;
      
      out = *((vec<T> *)a);
      
      return out;
    }

    // Integer
    inline int operator()(Integer *a){
      return *a;
    }
  };

  /////////////////////////////////////////////////////
  // Arithmetic operations
  /////////////////////////////////////////////////////
  struct Sum{
    // Complex/Real
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::r, 1)
      {
        out.v[i] = a.v[i] + b.v[i];
      }
      
      return out;
    }
    
    //I nteger
    inline int operator()(int a, int b){
      return a + b;
    }
  };

  struct Sub{
    // Complex/Real
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::r, 1)
      {
        out.v[i] = a.v[i] - b.v[i];
      }
      
      return out;
    }
    
    //Integer
    inline int operator()(int a, int b){
      return a-b;
    }
  };

  struct Mult{
    // Real
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::r, 1)
      {
        out.v[i] = a.v[i]*b.v[i];
      }
      
      return out;
    }
    
    // Integer
    inline int operator()(int a, int b){
      return a*b;
    }
  };
  
  #define cmul(a, b, c, i)\
  c[i]   = a[i]*b[i]   - a[i+1]*b[i+1];\
  c[i+1] = a[i]*b[i+1] + a[i+1]*b[i];
  
  struct MultComplex{
    // Complex
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::c, 1)
      {
        cmul(a.v, b.v, out.v, 2*i);
      }      
      
      return out;
    }
  };
  
  #undef cmul

  struct Div{
    // Real
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::r, 1)
      {
        out.v[i] = a.v[i]/b.v[i];
      }
      
      return out;
    }
  };
  
  #define conj(a, b, i)\
  b[i]   = a[i];\
  b[i+1] = -a[i+1];
  
  struct Conj{
    // Complex
    template <typename T>
    inline vec<T> operator()(vec<T> a){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::c, 1)
      {
        conj(a.v, out.v, 2*i);
      }
      
      return out;
    }
  };
  
  #undef conj

  #define timesmi(a, b, i)\
  b[i]   = a[i+1];\
  b[i+1] = -a[i];
  
  struct TimesMinusI{
    // Complex
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::c, 1)
      {
        timesmi(a.v, out.v, 2*i);
      }
      
      return out;
    }
  };

  #undef timesmi
  
  #define timesi(a, b, i)\
  b[i]   = -a[i+1];\
  b[i+1] = a[i];
  
  struct TimesI{
    // Complex
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::c, 1)
      {
        timesi(a.v, out.v, 2*i);
      }
      
      return out;
    }
  };
  
  #undef timesi

  //////////////////////////////////////////////
  // Some Template specialization
  #define perm(a, b, n, w)\
  unsigned int _mask = w >> (n + 1);\
  VECTOR_FOR(i, w, 1)\
  {\
    b[i] = a[i^_mask];\
  }
  
  #define DECL_PERMUTE_N(n)\
  template <typename T>\
  static inline vec<T> Permute##n(vec<T> in) {\
    vec<T> out;\
    perm(in.v, out.v, n, W<T>::r);\
    return out;\
  }
  
  struct Permute{
    DECL_PERMUTE_N(0);
    DECL_PERMUTE_N(1);
    DECL_PERMUTE_N(2);
    DECL_PERMUTE_N(3);
  };
  
  #undef perm
  #undef DECL_PERMUTE_N
  
  #define rot(a, b, n, w)\
  VECTOR_FOR(i, w, 1)\
  {\
    b[i] = a[(i + n)%w];\
  }
  
  struct Rotate{
    template <typename T>
    static inline vec<T> rotate(vec<T> in, int n){
      vec<T> out;
      
      rot(in.v, out.v, n, W<T>::r);
      
      return out;
    }
  };

  #undef rot
  
  #define acc(v, a, off, step, n)\
  for (unsigned int i = off; i < n; i += step)\
  {\
    a += v[i];\
  }
  
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
  
  //Complex float Reduce
  template <>
  inline Grid::ComplexF Reduce<Grid::ComplexF, vecf>::operator()(vecf in){
    float a = 0.f, b = 0.f;
    
    acc(in.v, a, 0, 2, W<float>::r);
    acc(in.v, b, 1, 2, W<float>::r);
    
    return Grid::ComplexF(a, b);
  }
  
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, vecf>::operator()(vecf in){
    float a = 0.;
    
    acc(in.v, a, 0, 1, W<float>::r);
    
    return a;
  }
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, vecd>::operator()(vecd in){
    double a = 0., b = 0.;
    
    acc(in.v, a, 0, 2, W<double>::r);
    acc(in.v, b, 1, 2, W<double>::r);
    
    return Grid::ComplexD(a, b);
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, vecd>::operator()(vecd in){
    double a = 0.f;
    
    acc(in.v, a, 0, 1, W<double>::r);
    
    return a;
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, int>::operator()(int in){
    return in;
  }
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 

  typedef Optimization::vecf SIMD_Ftype; // Single precision type
  typedef Optimization::vecd SIMD_Dtype; // Double precision type
  typedef int SIMD_Itype; // Integer type

  // prefetch utilities
  inline void v_prefetch0(int size, const char *ptr){};
  inline void prefetch_HINT_T0(const char *ptr){};

  // Function name aliases
  typedef Optimization::Vsplat   VsplatSIMD;
  typedef Optimization::Vstore   VstoreSIMD;
  typedef Optimization::Vset     VsetSIMD;
  typedef Optimization::Vstream  VstreamSIMD;
  template <typename S, typename T> using ReduceSIMD = Optimization::Reduce<S,T>;

  // Arithmetic operations
  typedef Optimization::Sum         SumSIMD;
  typedef Optimization::Sub         SubSIMD;
  typedef Optimization::Div         DivSIMD;
  typedef Optimization::Mult        MultSIMD;
  typedef Optimization::MultComplex MultComplexSIMD;
  typedef Optimization::Conj        ConjSIMD;
  typedef Optimization::TimesMinusI TimesMinusISIMD;
  typedef Optimization::TimesI      TimesISIMD;
}
