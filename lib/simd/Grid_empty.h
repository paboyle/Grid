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
    inline float operator()(float a, float b){
      return 0;
    }
    // Real float
    inline float operator()(float a){
      return 0;
    }
    //Complex double
    inline double operator()(double a, double b){
      return 0;
    }
    //Real double
    inline double operator()(double a){
      return 0;
    }
    //Integer
    inline int operator()(Integer a){
      return 0;
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(float a, float* F){
      
    }
    //Double
    inline void operator()(double a, double* D){
     
    }
    //Integer
    inline void operator()(int a, Integer* I){
      
    }

  };

  struct Vstream{
    //Float
    inline void operator()(float * a, float b){
     
    }
    //Double
    inline void operator()(double * a, double b){
     
    }


  };

  struct Vset{
    // Complex float 
    inline float operator()(Grid::ComplexF *a){
      return 0;
    }
    // Complex double 
    inline double operator()(Grid::ComplexD *a){
      return 0;
    }
    // Real float 
    inline float operator()(float *a){
      return  0;
    }
    // Real double
    inline double operator()(double *a){
      return 0;
    }
    // Integer
    inline int operator()(Integer *a){
      return 0;
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
    inline float operator()(float a, float b){
      return 0;
    }
    //Complex/Real double
    inline double operator()(double a, double b){
      return 0;
    }
    //Integer
    inline int operator()(int a, int b){
      return 0;
    }
  };

  struct Sub{
    //Complex/Real float
    inline float operator()(float a, float b){
      return 0;
    }
    //Complex/Real double
    inline double operator()(double a, double b){
      return 0;
    }
    //Integer
    inline int operator()(int a, int b){
      return 0;
    }
  };

  struct MultComplex{
    // Complex float
    inline float operator()(float a, float b){
      return 0;
    }
    // Complex double
    inline double operator()(double a, double b){
      return 0;
    }
  };

  struct Mult{
    inline float  mac(float a, float b,double c){
      return 0;
    }
    inline double mac(double a, double b,double c){
      return 0;
    }
    // Real float
    inline float operator()(float a, float b){
      return 0;
    }
    // Real double
    inline double operator()(double a, double b){
      return 0;
    }
    // Integer
    inline int operator()(int a, int b){
      return 0;
    }
  };

  struct Conj{
    // Complex single
    inline float operator()(float in){
      return 0;
    }
    // Complex double
    inline double operator()(double in){
      return 0;
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline float operator()(float in, float ret){
      return 0;
    }
    //Complex double
    inline double operator()(double in, double ret){
      return 0;
    }


  };

  struct TimesI{
    //Complex single
    inline float operator()(float in, float ret){
      return 0;
    }
    //Complex double
    inline double operator()(double in, double ret){
      return 0;
    }
  };

  //////////////////////////////////////////////
  // Some Template specialization
  template < typename vtype > 
    void permute(vtype &a, vtype b, int perm) {
   }; 

  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, float>::operator()(float in){
    return 0;
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, float>::operator()(float in){
    return 0;
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, double>::operator()(double in){
    return 0;
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, double>::operator()(double in){
    return 0;
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
namespace Grid {

  typedef float SIMD_Ftype;  // Single precision type
  typedef double SIMD_Dtype; // Double precision type
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
