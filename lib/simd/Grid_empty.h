    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_empty.h

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
//----------------------------------------------------------------------
/*! @file Grid_sse4.h
  @brief Empty Optimization libraries for debugging

  Using intrinsics
*/
// Time-stamp: <2015-06-09 14:28:02 neo>
//----------------------------------------------------------------------

namespace Grid {
namespace Optimization {

  template<class vtype>
  union uconv {
    float f;
    vtype v;
  };

  union u128f {
    float v;
    float f[4];
  };
  union u128d {
    double v;
    double f[2];
  };
  
  struct Vsplat{
    //Complex float
    inline u128f operator()(float a, float b){
      u128f out; 
      out.f[0] = a;
      out.f[1] = b;
      out.f[2] = a;
      out.f[3] = b;
      return out;
    }
    // Real float
    inline u128f operator()(float a){
      u128f out; 
      out.f[0] = a;
      out.f[1] = a;
      out.f[2] = a;
      out.f[3] = a;
      return out;
    }
    //Complex double
    inline u128d operator()(double a, double b){
      u128d out; 
      out.f[0] = a;
      out.f[1] = b;
      return out;
    }
    //Real double
    inline u128d operator()(double a){
      u128d out; 
      out.f[0] = a;
      out.f[1] = a;
      return out;
    }
    //Integer
    inline int operator()(Integer a){
      return a;
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(u128f a, float* F){
      memcpy(F,a.f,4*sizeof(float));
    }
    //Double
    inline void operator()(u128d a, double* D){
      memcpy(D,a.f,2*sizeof(double));
    }
    //Integer
    inline void operator()(int a, Integer* I){
      I[0] = a;
    }

  };

  struct Vstream{
    //Float
    inline void operator()(float * a, u128f b){
      memcpy(a,b.f,4*sizeof(float));
    }
    //Double
    inline void operator()(double * a, u128d b){
      memcpy(a,b.f,2*sizeof(double));
    }


  };

  struct Vset{
    // Complex float 
    inline u128f operator()(Grid::ComplexF *a){
      u128f out; 
      out.f[0] = a[0].real();
      out.f[1] = a[0].imag();
      out.f[2] = a[1].real();
      out.f[3] = a[1].imag();
      return out;
    }
    // Complex double 
    inline u128d operator()(Grid::ComplexD *a){
      u128d out; 
      out.f[0] = a[0].real();
      out.f[1] = a[0].imag();
      return out;
    }
    // Real float 
    inline u128f operator()(float *a){
      u128f out; 
      out.f[0] = a[0];
      out.f[1] = a[1];
      out.f[2] = a[2];
      out.f[3] = a[3];
      return out;
    }
    // Real double
    inline u128d operator()(double *a){
      u128d out; 
      out.f[0] = a[0];
      out.f[1] = a[1];
      return out;
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
  struct Sum{
    //Complex/Real float
    inline u128f operator()(u128f a, u128f b){
      u128f out;
      out.f[0] = a.f[0] + b.f[0];
      out.f[1] = a.f[1] + b.f[1];
      out.f[2] = a.f[2] + b.f[2];
      out.f[3] = a.f[3] + b.f[3];
      return out;
    }
    //Complex/Real double
    inline u128d operator()(u128d a, u128d b){
      u128d out;
      out.f[0] = a.f[0] + b.f[0];
      out.f[1] = a.f[1] + b.f[1];
      return out;
    }
    //Integer
    inline int operator()(int a, int b){
      return a + b;
    }
  };

  struct Sub{
    //Complex/Real float
    inline u128f operator()(u128f a, u128f b){
      u128f out;
      out.f[0] = a.f[0] - b.f[0];
      out.f[1] = a.f[1] - b.f[1];
      out.f[2] = a.f[2] - b.f[2];
      out.f[3] = a.f[3] - b.f[3];
      return out;
    }
    //Complex/Real double
    inline u128d operator()(u128d a, u128d b){
      u128d out;
      out.f[0] = a.f[0] - b.f[0];
      out.f[1] = a.f[1] - b.f[1];
      return out;
    }
    //Integer
    inline int operator()(int a, int b){
      return a-b;
    }
  };

  struct MultComplex{
    // Complex float
    inline u128f operator()(u128f a, u128f b){
      u128f out;
      out.f[0] = a.f[0]*b.f[0] - a.f[1]*b.f[1];
      out.f[1] = a.f[0]*b.f[1] + a.f[1]*b.f[0];
      out.f[2] = a.f[2]*b.f[2] - a.f[3]*b.f[3];
      out.f[3] = a.f[2]*b.f[3] + a.f[3]*b.f[2];
      return out;
    }
    // Complex double
    inline u128d operator()(u128d a, u128d b){
      u128d out;
      out.f[0] = a.f[0]*b.f[0] - a.f[1]*b.f[1];
      out.f[1] = a.f[0]*b.f[1] + a.f[1]*b.f[0];
      return out;
    }
  };

  struct Mult{
    //CK: Appear unneeded
    // inline float  mac(float a, float b,double c){
    //   return 0;
    // }
    // inline double mac(double a, double b,double c){
    //   return 0;
    // }

    // Real float
    inline u128f operator()(u128f a, u128f b){
      u128f out;
      out.f[0] = a.f[0]*b.f[0];
      out.f[1] = a.f[1]*b.f[1];
      out.f[2] = a.f[2]*b.f[2];
      out.f[3] = a.f[3]*b.f[3];
      return out;
    }
    // Real double
    inline u128d operator()(u128d a, u128d b){
      u128d out;
      out.f[0] = a.f[0]*b.f[0];
      out.f[1] = a.f[1]*b.f[1];
      return out;
    }
    // Integer
    inline int operator()(int a, int b){
      return a*b;
    }
  };

  struct Conj{
    // Complex single
    inline u128f operator()(u128f in){
      u128f out;
      out.f[0] = in.f[0];
      out.f[1] = -in.f[1];
      out.f[2] = in.f[2];
      out.f[3] = -in.f[3];
      return out;
    }
    // Complex double
    inline u128d operator()(u128d in){
      u128d out;
      out.f[0] = in.f[0];
      out.f[1] = -in.f[1];
      return out;
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline u128f operator()(u128f in, u128f ret){ //note ret is ignored
      u128f out;
      out.f[0] = in.f[1];
      out.f[1] = -in.f[0];
      out.f[2] = in.f[3];
      out.f[3] = -in.f[2];
      return out;
    }
    //Complex double
    inline u128d operator()(u128d in, u128d ret){
      u128d out;
      out.f[0] = in.f[1];
      out.f[1] = -in.f[0];
      return out;
    }
  };

  struct TimesI{
    //Complex single
    inline u128f operator()(u128f in, u128f ret){ //note ret is ignored
      u128f out;
      out.f[0] = -in.f[1];
      out.f[1] = in.f[0];
      out.f[2] = -in.f[3];
      out.f[3] = in.f[2];
      return out;
    }
    //Complex double
    inline u128d operator()(u128d in, u128d ret){
      u128d out;
      out.f[0] = -in.f[1];
      out.f[1] = in.f[0];
      return out;
    }
  };

  //////////////////////////////////////////////
  // Some Template specialization
  struct Permute{
    //We just have to mirror the permutes of Grid_sse4.h
    static inline u128f Permute0(u128f in){ //AB CD -> CD AB
      u128f out;
      out.f[0] = in.f[2];
      out.f[1] = in.f[3];
      out.f[2] = in.f[0];
      out.f[3] = in.f[1];
      return out;
    };
    static inline u128f Permute1(u128f in){ //AB CD -> BA DC
      u128f out;
      out.f[0] = in.f[1];
      out.f[1] = in.f[0];
      out.f[2] = in.f[3];
      out.f[3] = in.f[2];
      return out;
    };
    static inline u128f Permute2(u128f in){
      return in;
    };
    static inline u128f Permute3(u128f in){
      return in;
    };

    static inline u128d Permute0(u128d in){ //AB -> BA
      u128d out;
      out.f[0] = in.f[1];
      out.f[1] = in.f[0];
      return out;      
    };
    static inline u128d Permute1(u128d in){
      return in;
    };
    static inline u128d Permute2(u128d in){
      return in;
    };
    static inline u128d Permute3(u128d in){
      return in;
    };

  };
  
  template < typename vtype > 
    void permute(vtype &a, vtype b, int perm) {
   }; 

  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, u128f>::operator()(u128f in){ //2 complex
    return Grid::ComplexF(in.f[0] + in.f[2], in.f[1] + in.f[3]);
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, u128f>::operator()(u128f in){ //4 floats
    return in.f[0] + in.f[1] + in.f[2] + in.f[3];
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, u128d>::operator()(u128d in){ //1 complex
    return Grid::ComplexD(in.f[0],in.f[1]);
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, u128d>::operator()(u128d in){ //2 doubles
    return in.f[0] + in.f[1];
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, int>::operator()(int in){
    // FIXME unimplemented
   printf("Reduce : Missing integer implementation -> FIX\n");
    assert(0);
  }
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 

  typedef Optimization::u128f SIMD_Ftype;  // Single precision type
  typedef Optimization::u128d SIMD_Dtype; // Double precision type
  typedef int SIMD_Itype; // Integer type

  // prefetch utilities
  inline void v_prefetch0(int size, const char *ptr){};
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
