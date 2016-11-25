/*******************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/simd/Grid_qpx.h
 
 Copyright (C) 2016
 
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
 ******************************************************************************/

namespace Grid {
namespace Optimization {
  typedef struct 
  {
    float v0,v1,v2,v3;
  } vector4float;

  inline std::ostream & operator<<(std::ostream& stream, const vector4double a)
  {
    stream << "{"<<vec_extract(a,0)<<","<<vec_extract(a,1)<<","<<vec_extract(a,2)<<","<<vec_extract(a,3)<<"}";
    return stream;
  };

  inline std::ostream & operator<<(std::ostream& stream, const vector4float a)
  {
    stream << "{"<< a.v0 <<","<< a.v1 <<","<< a.v2 <<","<< a.v3 <<"}";
    return stream;
  };
  
  struct Vsplat{
    //Complex float
    inline vector4float operator()(float a, float b){
      return (vector4float){a, b, a, b};
    }
    // Real float
    inline vector4float operator()(float a){
      return (vector4float){a, a, a, a};
    }
    //Complex double
    inline vector4double operator()(double a, double b){
      return (vector4double){a, b, a, b};
    }
    //Real double
    inline vector4double operator()(double a){
      return (vector4double){a, a, a, a};
    }
    //Integer
    inline int operator()(Integer a){
      return a;
    }
  };
  
  struct Vstore{
    //Float
    inline void operator()(vector4double a, float *f){
      vec_st(a, 0, f);
    }

    inline void operator()(vector4double a, vector4float &f){
      vec_st(a, 0, (float *)(&f));
    }

    inline void operator()(vector4float a, float *f){
      f[0] = a.v0;
      f[1] = a.v1;
      f[2] = a.v2;
      f[3] = a.v3;
    }

    //Double
    inline void operator()(vector4double a, double *d){
      vec_st(a, 0, d);
    }
    //Integer
    inline void operator()(int a, Integer *i){
      i[0] = a;
    }
  };
  
  struct Vstream{
    //Float
    inline void operator()(float *f, vector4double a){
      vec_st(a, 0, f);
    }

    inline void operator()(vector4float f, vector4double a){
      vec_st(a, 0, (float *)(&f));
    }

    inline void operator()(float *f, vector4float a){
      f[0] = a.v0;
      f[1] = a.v1;
      f[2] = a.v2;
      f[3] = a.v3;
    }

    //Double
    inline void operator()(double *d, vector4double a){
      vec_st(a, 0, d);
    }

  };
  
  struct Vset{
    // Complex float
    inline vector4float operator()(Grid::ComplexF *a){
      return (vector4float){a[0].real(), a[0].imag(), a[1].real(), a[1].imag()};
    }
    // Complex double
    inline vector4double operator()(Grid::ComplexD *a){
      return vec_ld(0, (double *)a);
    }

    // Real float
    inline vector4float operator()(float *a){
      return (vector4float){a[0], a[1], a[2], a[3]};
    }

    inline vector4double operator()(vector4float a){
      return vec_ld(0, (float *)(&a));
    }

    // Real double
    inline vector4double operator()(double *a){
      return vec_ld(0, a);
    }
    // Integer
    inline int operator()(Integer *a){
      return a[0];
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
  #define FLOAT_WRAP_2(fn, pref)\
  pref vector4float fn(vector4float a, vector4float b)\
  {\
    vector4double ad, bd, rd;\
    vector4float  r;\
    \
    ad = Vset()(a);\
    bd = Vset()(b);\
    rd = fn(ad, bd);\
    Vstore()(rd, r);\
    \
    return r;\
  }

  #define FLOAT_WRAP_1(fn, pref)\
  pref vector4float fn(vector4float a)\
  {\
    vector4double ad, rd;\
    vector4float  r;\
    \
    ad = Vset()(a);\
    rd = fn(ad);\
    Vstore()(rd, r);\
    \
    return r;\
  }

  struct Sum{
    //Complex/Real double
    inline vector4double operator()(vector4double a, vector4double b){
      return vec_add(a, b);
    }

    //Complex/Real float
    FLOAT_WRAP_2(operator(), inline)

    //Integer
    inline int operator()(int a, int b){
      return a + b;
    }
  };
  
  struct Sub{
    //Complex/Real double
    inline vector4double operator()(vector4double a, vector4double b){
      return vec_sub(a, b);
    }

    //Complex/Real float
    FLOAT_WRAP_2(operator(), inline)

    //Integer
    inline int operator()(int a, int b){
      return a - b;
    }
  };
  
  struct MultComplex{
    // Complex double
    inline vector4double operator()(vector4double a, vector4double b){
      return vec_xxnpmadd(a, b, vec_xmul(b, a));
    }

    // Complex float
    FLOAT_WRAP_2(operator(), inline)
  };
  
  struct Mult{
    // Real double
    inline vector4double operator()(vector4double a, vector4double b){
      return vec_mul(a, b);
    }

    // Real float
    FLOAT_WRAP_2(operator(), inline)

    // Integer
    inline int operator()(int a, int b){
      return a*b;
    }
  };

  struct Div{
    // Real double
    inline vector4double operator()(vector4double a, vector4double b){
      return vec_swdiv(a, b);
    }

    // Real float
    FLOAT_WRAP_2(operator(), inline)

    // Integer
    inline int operator()(int a, int b){
      return a/b;
    }
  };

  struct Conj{
    // Complex double
    inline vector4double operator()(vector4double v){
      return vec_mul(v, (vector4double){1., -1., 1., -1.});
    }

    // Complex float
    FLOAT_WRAP_1(operator(), inline)
  };
  
  struct TimesMinusI{
    //Complex double
    inline vector4double operator()(vector4double v, vector4double ret){
      return vec_xxcpnmadd(v, (vector4double){1., 1., 1., 1.},
                               (vector4double){0., 0., 0., 0.});
    }

    // Complex float
    FLOAT_WRAP_2(operator(), inline)
  };
  
  struct TimesI{
    //Complex double
    inline vector4double operator()(vector4double v, vector4double ret){
      return vec_xxcpnmadd(v, (vector4double){-1., -1., -1., -1.},
                              (vector4double){0., 0., 0., 0.});
    }

    // Complex float
    FLOAT_WRAP_2(operator(), inline)
  };

  struct Permute{
    //Complex double
    static inline vector4double Permute0(vector4double v){ //0123 -> 2301
      return vec_perm(v, v, vec_gpci(02301));
    };
    static inline vector4double Permute1(vector4double v){ //0123 -> 1032
      return vec_perm(v, v, vec_gpci(01032));
    };
    static inline vector4double Permute2(vector4double v){
      return v;
    };
    static inline vector4double Permute3(vector4double v){
      return v;
    };

    // Complex float
    FLOAT_WRAP_1(Permute0, static inline)
    FLOAT_WRAP_1(Permute1, static inline)
    FLOAT_WRAP_1(Permute2, static inline)
    FLOAT_WRAP_1(Permute3, static inline)
  };
  
  struct Rotate{
    static inline vector4double rotate(vector4double v, int n){
      switch(n){
        case 0:
          return v;
          break;
        case 1:
          return vec_perm(v, v, vec_gpci(01230));
          break;
        case 2:
          return vec_perm(v, v, vec_gpci(02301));
          break;
        case 3:
          return vec_perm(v, v, vec_gpci(03012));
          break;
        default: assert(0);
      }
    }

    static inline vector4float rotate(vector4float v, int n){
      vector4double vd, rd;
      vector4float  r;

      vd = Vset()(v);
      rd = rotate(vd, n);
      Vstore()(rd, r);

      return r;
    }
  };
  
  //Complex float Reduce
  template<>
  inline Grid::ComplexF
  Reduce<Grid::ComplexF, vector4float>::operator()(vector4float v) { //2 complex
    vector4float v1,v2;
    
    v1 = Optimization::Permute::Permute0(v);
    v1 = Optimization::Sum()(v1, v);
    
    return Grid::ComplexF(v1.v0, v1.v1);
  }
  //Real float Reduce
  template<>
  inline Grid::RealF
  Reduce<Grid::RealF, vector4float>::operator()(vector4float v){ //4 floats
    vector4float v1,v2;
    
    v1 = Optimization::Permute::Permute0(v);
    v1 = Optimization::Sum()(v1, v);
    v2 = Optimization::Permute::Permute1(v1);
    v1 = Optimization::Sum()(v1, v2);
    
    return v1.v0;
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD
  Reduce<Grid::ComplexD, vector4double>::operator()(vector4double v){ //2 complex
    vector4double v1;
    
    v1 = Optimization::Permute::Permute0(v);
    v1 = vec_add(v1, v);
    
    return Grid::ComplexD(vec_extract(v1, 0), vec_extract(v1, 1));
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD
  Reduce<Grid::RealD, vector4double>::operator()(vector4double v){ //4 doubles
    vector4double v1,v2;
    
    v1 = Optimization::Permute::Permute0(v);
    v1 = vec_add(v1, v);
    v2 = Optimization::Permute::Permute1(v1);
    v1 = vec_add(v1, v2);

    return vec_extract(v1, 0);
  }
  
  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, int>::operator()(int in){
    // FIXME unimplemented
    printf("Reduce : Missing integer implementation -> FIX\n");
    assert(0);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Here assign types
typedef Optimization::vector4float SIMD_Ftype;  // Single precision type
typedef vector4double              SIMD_Dtype; // Double precision type
typedef int                        SIMD_Itype; // Integer type

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
typedef Optimization::Mult        MultSIMD;
typedef Optimization::Div         DivSIMD;
typedef Optimization::MultComplex MultComplexSIMD;
typedef Optimization::Conj        ConjSIMD;
typedef Optimization::TimesMinusI TimesMinusISIMD;
typedef Optimization::TimesI      TimesISIMD;
  
}
