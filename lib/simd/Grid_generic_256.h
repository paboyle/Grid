    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_generic.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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

#ifndef GEN_SIMD_WIDTH
#define GEN_SIMD_DCOMPLEX_WIDTH 2
#endif

#include "Grid_generic.h"

namespace Grid {
namespace Optimization {

  constexpr unsigned int dcw = GEN_SIMD_DCOMPLEX_WIDTH;
  constexpr unsigned int fcw = 2*dcw;
  constexpr unsigned int dw  = 2*dcw;
  constexpr unsigned int fw  = 2*fcw;

  struct vecf {
    float v[fw];
  };
  
  struct vecd {
    double v[dw];
  };
  
  struct Vsplat{
    //Complex float
    inline vecf operator()(float a, float b){
      vecf out;
      
      for (unsigned int i = 0; i < fw; i += 2)
      {
        out.v[i]   = a;
        out.v[i+1] = b;
      }

      return out;
    }
    
    // Real float
    inline vecf operator()(float a){
      vecf out;
      
      for (unsigned int i = 0; i < fw; ++i)
      {
        out.v[i] = a;
      }
      
      return out;
    }
    
    //Complex double
    inline vecd operator()(double a, double b){
      vecd out;
      
      for (unsigned int i = 0; i < dw; i += 2)
      {
        out.v[i]   = a;
        out.v[i+1] = b;
      }
      
      return out;
    }
    
    //Real double
    inline vecd operator()(double a){
      vecd out;
      
      for (unsigned int i = 0; i < dw; ++i)
      {
        out.v[i] = a;
      }
      
      return out;
    }
    
    //Integer
    inline int operator()(Integer a){
      return a;
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(vecf a, float* F){
      memcpy(F,a.v,fw*sizeof(float));
    }
    //Double
    inline void operator()(vecd a, double* D){
      memcpy(D,a.v,dw*sizeof(double));
    }
    //Integer
    inline void operator()(int a, Integer* I){
      I[0] = a;
    }

  };

  struct Vstream{
    //Float
    inline void operator()(float * a, vecf b){
      memcpy(a,b.v,fw*sizeof(float));
    }
    //Double
    inline void operator()(double * a, vecd b){
      memcpy(a,b.v,dw*sizeof(double));
    }


  };

  struct Vset{
    // Complex float 
    inline vecf operator()(Grid::ComplexF *a){
      vecf out;
      
      for (unsigned int i = 0; i < fcw; ++i)
      {
        out.v[2*i]   = a[i].real();
        out.v[2*i+1] = a[i].imag();
      }
      
      return out;
    }
    
    // Complex double 
    inline vecd operator()(Grid::ComplexD *a){
      vecd out;
      
      for (unsigned int i = 0; i < dcw; ++i)
      {
        out.v[2*i]   = a[i].real();
        out.v[2*i+1] = a[i].imag();
      }
      
      return out;
    }
    
    // Real float 
    inline vecf operator()(float *a){
      vecf out;
      
      memcpy(out.v,a,fw*sizeof(float));
      
      return out;
    }
    // Real double
    inline vecd operator()(double *a){
      vecd out; 
      
      memcpy(out.v,a,dw*sizeof(float));
      
      return out;
    }
    // Integer
    inline int operator()(Integer *a){
      return a[0];
    }


  };

  /////////////////////////////////////////////////////
  // Arithmetic operations
  /////////////////////////////////////////////////////
  struct Sum{
    //Complex/Real float
    inline vecf operator()(vecf a, vecf b){
      vecf out;
      
      for (unsigned int i = 0; i < fw; ++i)
      {
        out.v[i] = a.v[i] + b.v[i];
      }
      
      return out;
    }
    
    //Complex/Real double
    inline vecd operator()(vecd a, vecd b){
      vecd out;
      
      for (unsigned int i = 0; i < dw; ++i)
      {
        out.v[i] = a.v[i] + b.v[i];
      }

      return out;
    }
    
    //Integer
    inline int operator()(int a, int b){
      return a + b;
    }
  };

  struct Sub{
    //Complex/Real float
    inline vecf operator()(vecf a, vecf b){
      vecf out;
      
      for (unsigned int i = 0; i < fw; ++i)
      {
        out.v[i] = a.v[i] - b.v[i];
      }
      
      return out;
    }
    
    //Complex/Real double
    inline vecd operator()(vecd a, vecd b){
      vecd out;
      
      for (unsigned int i = 0; i < dw; ++i)
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

  #define cmul(a, b, c, i)\
  c[i]   = a[i]*b[i]   - a[i+1]*b[i+1];\
  c[i+1] = a[i]*b[i+1] + a[i+1]*b[i];
  
  struct MultComplex{
    // Complex float
    inline vecf operator()(vecf a, vecf b){
      vecf out;
      
      for (unsigned int i = 0; i < fcw; ++i)
      {
        cmul(a.v, b.v, out.v, 2*i);
      }
      
      return out;
    }
    
    // Complex double
    inline vecd operator()(vecd a, vecd b){
      vecd out;
      
      for (unsigned int i = 0; i < dcw; ++i)
      {
        cmul(a.v, b.v, out.v, 2*i);
      }
      
      return out;
    }
  };
  
  #undef cmul

  struct Mult{
    // Real float
    inline vecf operator()(vecf a, vecf b){
      vecf out;
      
      for (unsigned int i = 0; i < fw; ++i)
      {
        out.v[i] = a.v[i]*b.v[i];
      }
      
      return out;
    }
    
    // Real double
    inline vecd operator()(vecd a, vecd b){
      vecd out;
      
      for (unsigned int i = 0; i < dw; ++i)
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

  struct Div{
    // Real float
    inline vecf operator()(vecf a, vecf b){
      vecf out;
      
      for (unsigned int i = 0; i < fw; ++i)
      {
        out.v[i] = a.v[i]/b.v[i];
      }
      
      return out;
    }
    // Real double
    inline vecd operator()(vecd a, vecd b){
      vecd out;
      
      for (unsigned int i = 0; i < dw; ++i)
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
    // Complex single
    inline vecf operator()(vecf in){
      vecf out;
      
      for (unsigned int i = 0; i < fcw; ++i)
      {
        conj(in.v, out.v, 2*i);
      }
      
      return out;
    }
    
    // Complex double
    inline vecd operator()(vecd in){
      vecd out;
      
      for (unsigned int i = 0; i < dcw; ++i)
      {
        conj(in.v, out.v, 2*i);
      }
      
      return out;
    }
  };
  
  #undef conj

  #define timesmi(a, b, i)\
  b[i]   = a[i+1];\
  b[i+1] = -a[i];
  
  struct TimesMinusI{
    // Complex single
    inline vecf operator()(vecf in, vecf ret){
      vecf out;
      
      for (unsigned int i = 0; i < fcw; ++i)
      {
        timesmi(in.v, out.v, 2*i);
      }
      
      return out;
    }
    
    // Complex double
    inline vecd operator()(vecd in, vecd ret){
      vecd out;
      
      for (unsigned int i = 0; i < dcw; ++i)
      {
        timesmi(in.v, out.v, 2*i);
      }
      
      return out;
    }
  };

  #undef timesmi
  
  #define timespi(a, b, i)\
  b[i]   = -a[i+1];\
  b[i+1] = a[i];
  
  struct TimesI{
    // Complex single
    inline vecf operator()(vecf in, vecf ret){
      vecf out;
      
      for (unsigned int i = 0; i < fcw; ++i)
      {
        timespi(in.v, out.v, 2*i);
      }
      
      return out;
    }
    
    // Complex double
    inline vecd operator()(vecd in, vecd ret){
      vecd out;
      
      for (unsigned int i = 0; i < dcw; ++i)
      {
        timespi(in.v, out.v, 2*i);
      }
      
      return out;
    }
  };
  
  #undef timespi

  //////////////////////////////////////////////
  // Some Template specialization
  struct Permute{
    static inline vecf Permute0(vecf in){ //AB CD -> CD AB
      vecf out;
      
      out.v[0] = in.v[4];
      out.v[1] = in.v[5];
      out.v[2] = in.v[6];
      out.v[3] = in.v[7];
      out.v[4] = in.v[0];
      out.v[5] = in.v[1];
      out.v[6] = in.v[2];
      out.v[7] = in.v[3];
      
      return out;
    };
    
    static inline vecf Permute1(vecf in){ //AB CD -> BA DC
      vecf out;
      
      out.v[0] = in.v[2];
      out.v[1] = in.v[3];
      out.v[2] = in.v[0];
      out.v[3] = in.v[1];
      out.v[4] = in.v[6];
      out.v[5] = in.v[7];
      out.v[6] = in.v[4];
      out.v[7] = in.v[5];
      
      return out;
    };
    
    static inline vecf Permute2(vecf in){
      vecf out;
      
      out.v[0] = in.v[1];
      out.v[1] = in.v[0];
      out.v[2] = in.v[3];
      out.v[3] = in.v[2];
      out.v[4] = in.v[5];
      out.v[5] = in.v[4];
      out.v[6] = in.v[7];
      out.v[7] = in.v[6];
      
      return out;
    };

    static inline vecf Permute3(vecf in){
      return in;
    };

    static inline vecd Permute0(vecd in){ //AB -> BA
      vecd out;
      
      out.v[0] = in.v[2];
      out.v[1] = in.v[3];
      out.v[2] = in.v[0];
      out.v[3] = in.v[1];
      
      return out;
    };
    
    static inline vecd Permute1(vecd in){
      vecd out;
      
      out.v[0] = in.v[1];
      out.v[1] = in.v[0];
      out.v[2] = in.v[3];
      out.v[3] = in.v[2];
      
      return out;
    };
    
    static inline vecd Permute2(vecd in){
      return in;
    };
    
    static inline vecd Permute3(vecd in){
      return in;
    };

  };
  
  #define rot(a, b, n, w)\
  for (unsigned int i = 0; i < w; ++i)\
  {\
    b[i] = a[(i + n)%w];\
  }
  
  struct Rotate{

    static inline vecf rotate(vecf in, int n){
      vecf out;
      
      rot(in.v, out.v, n, fw);
      
      return out;
    }
    
    static inline vecd rotate(vecd in,int n){
      vecd out;
      
      rot(in.v, out.v, n, dw);
      
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
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, vecf>::operator()(vecf in){
    float a = 0.f, b = 0.f;
    
    acc(in.v, a, 0, 2, fw);
    acc(in.v, b, 1, 2, fw);
    
    return Grid::ComplexF(a, b);
  }
  
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, vecf>::operator()(vecf in){
    float a = 0.;
    
    acc(in.v, a, 0, 1, fw);
    
    return a;
  }
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, vecd>::operator()(vecd in){
    double a = 0., b = 0.;
    
    acc(in.v, a, 0, 2, dw);
    acc(in.v, b, 1, 2, dw);
    
    return Grid::ComplexD(a, b);
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, vecd>::operator()(vecd in){
    double a = 0.f;
    
    acc(in.v, a, 0, 1, dw);
    
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

  typedef Optimization::vecf SIMD_Ftype;  // Single precision type
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
