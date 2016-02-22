    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_qpx.h

    Copyright (C) 2015

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
/*! @file Grid_qpx.h
  @brief Optimization libraries for QPX instructions set for BG/Q

  Using intrinsics
*/
// Time-stamp: <2015-05-27 11:30:21 neo>
//----------------------------------------------------------------------

// lot of undefined functions

namespace Optimization {
  
  struct Vsplat{
    //Complex float
    inline float operator()(float a, float b){
      return {a,b,a,b};
    }
    // Real float
    inline float operator()(float a){
      return {a,a,a,a};
    }
    //Complex double
    inline vector4double operator()(double a, double b){
      return {a,b,a,b};
    }
    //Real double
    inline vector4double operator()(double a){
      return {a,a,a,a};
    }
    //Integer
    inline int operator()(Integer a){
#error
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(float a, float* F){
      assert(0);
    }
    //Double
    inline void operator()(vector4double a, double* D){
      assert(0);
    }
    //Integer
    inline void operator()(int a, Integer* I){
      assert(0);
    }

  };


  struct Vstream{
    //Float
    inline void operator()(float * a, float b){
      assert(0);
    }
    //Double
    inline void operator()(double * a, vector4double b){
      assert(0);
    }


  };



  struct Vset{
    // Complex float 
    inline float operator()(Grid::ComplexF *a){
      return {a[0].real(),a[0].imag(),a[1].real(),a[1].imag(),a[2].real(),a[2].imag(),a[3].real(),a[3].imag()};
    }
    // Complex double 
    inline vector4double operator()(Grid::ComplexD *a){
      return {a[0].real(),a[0].imag(),a[1].real(),a[1].imag(),a[2].real(),a[2].imag(),a[3].real(),a[3].imag()};
    }
    // Real float 
    inline float operator()(float *a){
      return {a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]};
    }
    // Real double
    inline vector4double operator()(double *a){
      return {a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]};
    }
    // Integer
    inline int operator()(Integer *a){
#error
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
#error
    }
    //Complex/Real double
    inline vector4double operator()(vector4double a, vector4double b){
      return vec_add(a,b);
    }
    //Integer
    inline int operator()(int a, int b){
#error
    }
  };

  struct Sub{
    //Complex/Real float
    inline float operator()(float a, float b){
#error
    }
    //Complex/Real double
    inline vector4double operator()(vector4double a, vector4double b){
#error
    }
    //Integer
    inline floati operator()(int a, int b){
#error
    }
  };


  struct MultComplex{
    // Complex float
    inline float operator()(float a, float b){
#error
    }
    // Complex double
    inline vector4double operator()(vector4double a, vector4double b){
#error
    }
  };

  struct Mult{
    // Real float
    inline float operator()(float a, float b){
#error
    }
    // Real double
    inline vector4double operator()(vector4double a, vector4double b){
#error
    }
    // Integer
    inline int operator()(int a, int b){
#error
    }
  };


  struct Conj{
    // Complex single
    inline float operator()(float in){
      assert(0);
    }
    // Complex double
    inline vector4double operator()(vector4double in){
      assert(0);
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline float operator()(float in, float ret){
      assert(0);
    }
    //Complex double
    inline vector4double operator()(vector4double in, vector4double ret){
      assert(0);
    }


  };

  struct TimesI{
    //Complex single
    inline float operator()(float in, float ret){
  
    }
    //Complex double
    inline vector4double operator()(vector4double in, vector4double ret){
  
    }


  };


  


  //////////////////////////////////////////////
  // Some Template specialization
  
  //Complex float Reduce
  template<>
    inline Grid::ComplexF Reduce<Grid::ComplexF, float>::operator()(float in){
    assert(0);
  }
  //Real float Reduce
  template<>
    inline Grid::RealF Reduce<Grid::RealF, float>::operator()(float in){
    assert(0);
  }
  
  
  //Complex double Reduce
  template<>
    inline Grid::ComplexD Reduce<Grid::ComplexD, vector4double>::operator()(vector4double in){
    assert(0);
  }
  
  //Real double Reduce
  template<>
    inline Grid::RealD Reduce<Grid::RealD, vector4double>::operator()(vector4double in){
    assert(0);
  }

  //Integer Reduce
  template<>
    inline Integer Reduce<Integer, floati>::operator()(float in){
    assert(0);
  }
  
  
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 
namespace Grid {
  typedef float SIMD_Ftype  __attribute__ ((vector_size (16)));         // Single precision type
  typedef vector4double SIMD_Dtype; // Double precision type
  typedef int SIMD_Itype;           // Integer type

  inline void v_prefetch0(int size, const char *ptr){};

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
