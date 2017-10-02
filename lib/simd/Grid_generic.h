    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_generic.h

    Copyright (C) 2015
    Copyright (C) 2017

Author: Antonin Portelli <antonin.portelli@me.com>
        Andrew Lawson    <andrew.lawson1991@gmail.com>

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

#include "Grid_generic_types.h"

namespace Grid {
namespace Optimization {
  
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
  };

  struct Vstore{
    // Real
    template <typename T>
    inline void operator()(vec<T> a, T *D){
      *((vec<T> *)D) = a;
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
  };
  
  #define cmul(a, b, c, i)\
  c[i]   = a[i]*b[i]   - a[i+1]*b[i+1];\
  c[i+1] = a[i]*b[i+1] + a[i+1]*b[i];

  struct MultRealPart{
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::c, 1)
      {
         out.v[2*i]   = a.v[2*i]*b.v[2*i];
         out.v[2*i+1] = a.v[2*i]*b.v[2*i+1];
      }      
      return out;
    }
  };

  struct MaddRealPart{
    template <typename T>
    inline vec<T> operator()(vec<T> a, vec<T> b, vec<T> c){
      vec<T> out;
      
      VECTOR_FOR(i, W<T>::c, 1)
      {
         out.v[2*i]   = a.v[2*i]*b.v[2*i] + c.v[2*i];
         out.v[2*i+1] = a.v[2*i]*b.v[2*i+1] + c.v[2*i+1];
      }      
      return out;
    }
  };
  
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

  struct PrecisionChange {
    static inline vech StoH (const vecf &a,const vecf &b) {
      vech ret;
#ifdef USE_FP16
      vech *ha = (vech *)&a;
      vech *hb = (vech *)&b;
      const int nf = W<float>::r;
      //      VECTOR_FOR(i, nf,1){ ret.v[i]    = ( (uint16_t *) &a.v[i])[1] ; }
      //      VECTOR_FOR(i, nf,1){ ret.v[i+nf] = ( (uint16_t *) &b.v[i])[1] ; }
      VECTOR_FOR(i, nf,1){ ret.v[i]    = ha->v[2*i+1]; }
      VECTOR_FOR(i, nf,1){ ret.v[i+nf] = hb->v[2*i+1]; }
#else
      assert(0);
#endif
      return ret;
    }
    static inline void  HtoS (vech h,vecf &sa,vecf &sb) {
#ifdef USE_FP16
      const int nf = W<float>::r;
      const int nh = W<uint16_t>::r;
      vech *ha = (vech *)&sa;
      vech *hb = (vech *)&sb;
      VECTOR_FOR(i, nf, 1){ sb.v[i]= sa.v[i] = 0; }
      //      VECTOR_FOR(i, nf, 1){ ( (uint16_t *) (&sa.v[i]))[1] = h.v[i];}
      //      VECTOR_FOR(i, nf, 1){ ( (uint16_t *) (&sb.v[i]))[1] = h.v[i+nf];}
      VECTOR_FOR(i, nf, 1){ ha->v[2*i+1]=h.v[i]; }
      VECTOR_FOR(i, nf, 1){ hb->v[2*i+1]=h.v[i+nf]; }
#else
      assert(0);
#endif
    }
    static inline vecf DtoS (vecd a,vecd b) {
      const int nd = W<double>::r;
      const int nf = W<float>::r;
      vecf ret;
      VECTOR_FOR(i, nd,1){ ret.v[i]    = a.v[i] ; }
      VECTOR_FOR(i, nd,1){ ret.v[i+nd] = b.v[i] ; }
      return ret;
    }
    static inline void StoD (vecf s,vecd &a,vecd &b) {
      const int nd = W<double>::r;
      VECTOR_FOR(i, nd,1){ a.v[i] = s.v[i] ; }
      VECTOR_FOR(i, nd,1){ b.v[i] = s.v[i+nd] ; }
    }
    static inline vech DtoH (vecd a,vecd b,vecd c,vecd d) {
      vecf sa,sb;
      sa = DtoS(a,b);
      sb = DtoS(c,d);
      return StoH(sa,sb);
    }
    static inline void HtoD (vech h,vecd &a,vecd &b,vecd &c,vecd &d) {
      vecf sa,sb;
      HtoS(h,sa,sb);
      StoD(sa,a,b);
      StoD(sb,c,d);
    }
  };

  //////////////////////////////////////////////
  // Exchange support
  struct Exchange{

    template <typename T,int n>
    static inline void ExchangeN(vec<T> &out1,vec<T> &out2,vec<T> &in1,vec<T> &in2){
      const int w = W<T>::r;
      unsigned int mask = w >> (n + 1);
      //      std::cout << " Exchange "<<n<<" nsimd "<<w<<" mask 0x" <<std::hex<<mask<<std::dec<<std::endl;
      VECTOR_FOR(i, w, 1) {	
	int j1 = i&(~mask);
	if  ( (i&mask) == 0 ) { out1.v[i]=in1.v[j1];}
	else                  { out1.v[i]=in2.v[j1];}
	int j2 = i|mask;
	if  ( (i&mask) == 0 ) { out2.v[i]=in1.v[j2];}
	else                  { out2.v[i]=in2.v[j2];}
      }      
    }
    template <typename T>
    static inline void Exchange0(vec<T> &out1,vec<T> &out2,vec<T> &in1,vec<T> &in2){
      ExchangeN<T,0>(out1,out2,in1,in2);
    };
    template <typename T>
    static inline void Exchange1(vec<T> &out1,vec<T> &out2,vec<T> &in1,vec<T> &in2){
      ExchangeN<T,1>(out1,out2,in1,in2);
    };
    template <typename T>
    static inline void Exchange2(vec<T> &out1,vec<T> &out2,vec<T> &in1,vec<T> &in2){
      ExchangeN<T,2>(out1,out2,in1,in2);
    };
    template <typename T>
    static inline void Exchange3(vec<T> &out1,vec<T> &out2,vec<T> &in1,vec<T> &in2){
      ExchangeN<T,3>(out1,out2,in1,in2);
    };
  };


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
      
    template <int n, typename T> static inline vec<T> tRotate(vec<T> in){
      return rotate(in, n);
    }
    
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
  inline Integer Reduce<Integer, veci>::operator()(veci in){
    Integer a = 0;
    
    acc(in.v, a, 0, 1, W<Integer>::r);
    
    return a;
  }

  #undef acc  // EIGEN compatibility
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 

  typedef Optimization::vech SIMD_Htype; // Reduced precision type
  typedef Optimization::vecf SIMD_Ftype; // Single precision type
  typedef Optimization::vecd SIMD_Dtype; // Double precision type
  typedef Optimization::veci SIMD_Itype; // Integer type

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
  typedef Optimization::MultRealPart MultRealPartSIMD;
  typedef Optimization::MaddRealPart MaddRealPartSIMD;
  typedef Optimization::Conj        ConjSIMD;
  typedef Optimization::TimesMinusI TimesMinusISIMD;
  typedef Optimization::TimesI      TimesISIMD;
}
