    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_gpu.h

    Copyright (C) 2018

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
/*! @file Grid_gpu.h
  @brief Optimization libraries for GPU
  Use float4, double2 
*/
//----------------------------------------------------------------------

#include <cuda_fp16.h>

namespace Grid {

  // re im, re, im, re, im etc..
struct half8 {
  half ax, ay, az, aw, bx, by, bz, bw;
};
accelerator_inline float half2float(half h)
{
  float f;
#ifdef __CUDA_ARCH__
  f = __half2float(h);
#else 
  //f = __half2float(h);
  __half_raw hr(h);
  Grid_half hh; 
  hh.x = hr.x;
  f=  sfw_half_to_float(hh);
#endif
  return f;
}
accelerator_inline half float2half(float f)
{
  half h;
#ifdef __CUDA_ARCH__
  h = __float2half(f);
#else
  Grid_half hh = sfw_float_to_half(f);
  __half_raw hr;  
  hr.x = hh.x;
  h = __half(hr);
#endif
  return h;
}

namespace Optimization {
  
  inline accelerator float4 operator*(float4 a,float4 b) {return make_float4(a.x*b.x,a.y*b.y,a.z*b.z,a.w*b.w);}
  inline accelerator float4 operator+(float4 a,float4 b) {return make_float4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w);}
  inline accelerator float4 operator-(float4 a,float4 b) {return make_float4(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w);}
  inline accelerator float4 operator/(float4 a,float4 b) {return make_float4(a.x/b.x,a.y/b.y,a.z/b.z,a.w/b.w);}

  inline accelerator double2 operator*(double2 a,double2 b) {return make_double2(a.x*b.x,a.y*b.y);}
  inline accelerator double2 operator+(double2 a,double2 b) {return make_double2(a.x+b.x,a.y+b.y);}
  inline accelerator double2 operator-(double2 a,double2 b) {return make_double2(a.x-b.x,a.y-b.y);}
  inline accelerator double2 operator/(double2 a,double2 b) {return make_double2(a.x/b.x,a.y/b.y);}

  inline accelerator int4 operator*(int4 a,int4 b) {return make_int4(a.x*b.x,a.y*b.y,a.z*b.z,a.w*b.w);}
  inline accelerator int4 operator+(int4 a,int4 b) {return make_int4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w);}
  inline accelerator int4 operator-(int4 a,int4 b) {return make_int4(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w);}
  inline accelerator int4 operator/(int4 a,int4 b) {return make_int4(a.x/b.x,a.y/b.y,a.z/b.z,a.w/b.w);}

  struct Vsplat{
    //Complex float
    accelerator_inline float4 operator()(float a, float b){
      float4 ret;
      ret.x=ret.z=a;
      ret.y=ret.w=b;
      return ret;
    }
    // Real float
    accelerator_inline float4 operator()(float a){
      float4 ret;
      ret.x=ret.y=ret.z=ret.w = a;
      return ret;
    }
    //Complex double
    accelerator_inline double2 operator()(double a, double b){
      double2 ret;
      ret.x=a;
      ret.y=b;
      return ret;
    }
    //Real double
    accelerator_inline double2 operator()(double a){
      double2 ret; 
      ret.x = ret.y = a;
      return ret;
    }
    //Integer
    accelerator_inline int4 operator()(Integer a){
      int4 ret;
      ret.x=ret.y=ret.z=ret.w=a;
      return ret;
    }
  };

  struct Vstore{
    //Float 
    accelerator_inline void operator()(float4 a, float* F){
      float4 *F4 = (float4 *)F;
      *F4 = a;
    }
    //Double
    accelerator_inline void operator()(double2 a, double* D){
      double2 *D2 = (double2 *)D;
      *D2 = a;
    }
    //Integer
    accelerator_inline void operator()(int4 a, Integer* I){
      int4 *I4 = (int4 *)I;
      *I4 = a;
    }

  };

  struct Vstream{
    //Float
    accelerator_inline void operator()(float * a, float4 b){
      float4 * a4 = (float4 *)a;
      *a4 = b;
    }
    //Double
    accelerator_inline void operator()(double * a, double2 b){
      double2 * a2 = (double2 *)a;
      *a2 = b;
    }

  };

  struct Vset{
    // Complex float 
    accelerator_inline float4 operator()(Grid::ComplexF *a){
      float4 ret;
      ret.x = a[0].real();
      ret.y = a[0].imag();
      ret.z = a[1].real();
      ret.w = a[1].imag();
      return ret;
    }
    // Complex double 
    accelerator_inline double2 operator()(Grid::ComplexD *a){
      double2 ret;
      ret.x = a[0].real();
      ret.y = a[0].imag();
      return ret;
    }
    // Real float 
    accelerator_inline float4 operator()(float *a){
      float4 ret;
      ret.x = a[0];
      ret.y = a[1];
      ret.z = a[2];
      ret.w = a[3];
      return ret;
    }
    // Real double
    accelerator_inline double2 operator()(double *a){
      double2 ret;
      ret.x = a[0];
      ret.y = a[1];
      return ret;
    }
    // Integer
    accelerator_inline int4 operator()(Integer *a){
      int4 ret;
      ret.x = a[0];
      ret.y = a[1];
      ret.z = a[2];
      ret.w = a[3];
      return ret;
    }

  };

  template <typename Out_type, typename In_type>
  struct Reduce{
    //Need templated class to overload output type
    //General form must generate error if compiled
    accelerator_inline Out_type operator()(In_type in){
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
    accelerator_inline float4 operator()(float4 a, float4 b){
      return a+b;
    }
    //Complex/Real double
    accelerator_inline double2 operator()(double2 a, double2 b){
      return a+b;
    }
    //Integer
    accelerator_inline int4 operator()(int4 a,int4 b){
      return a+b;
    }
  };

  struct Sub{
    //Complex/Real float
    accelerator_inline float4 operator()(float4 a, float4 b){
      return a-b;
    }
    //Complex/Real double
    accelerator_inline double2 operator()(double2 a, double2 b){
      return a-b;
    }
    //Integer
    accelerator_inline int4 operator()(int4 a, int4 b){
      return a-b;
    }
  };

  struct MultRealPart{
    accelerator_inline float4 operator()(float4 a, float4 b){
      float4 ymm0;
      ymm0.x = a.y;
      ymm0.y = a.y;
      ymm0.z = a.w;
      ymm0.w = a.w;
      return  ymm0*b;
      // ymm0 = _mm_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(2,2,0,0)); // ymm0 <- ar ar,
      // return _mm_mul_ps(ymm0,b);                       // ymm0 <- ar bi, ar br
    }
    accelerator_inline double2 operator()(double2 a, double2 b){
      double2 ymm0;
      ymm0.x = a.y;
      ymm0.y = a.y;
      return ymm0*b;
      //      ymm0 = _mm_shuffle_pd(a,a,0x0); // ymm0 <- ar ar, ar,ar b'00,00
      //      return _mm_mul_pd(ymm0,b);      // ymm0 <- ar bi, ar br
    }
  };
  struct MaddRealPart{
    accelerator_inline float4 operator()(float4 a, float4 b, float4 c){
      float4 ymm0; // =  _mm_shuffle_ps(a,a,_MM_SELECT_FOUR_FOUR(2,2,0,0)); // ymm0 <- ar ar,
      ymm0.x = a.y;
      ymm0.y = a.y;
      ymm0.z = a.w;
      ymm0.w = a.w;
      return c+ymm0*b;
    }
    accelerator_inline double2 operator()(double2 a, double2 b, double2 c){
      //      ymm0 = _mm_shuffle_pd( a, a, 0x0 );
      double2 ymm0;
      ymm0.x = a.y;
      ymm0.y = a.y;
      return c+ymm0*b;
    }
  };

  struct MultComplex{
    // Complex float
    accelerator_inline float4 operator()(float4 a, float4 b){
      float4 ymm0;
      ymm0.x = a.x*b.x -  a.y*b.y ; // rr - ii
      ymm0.y = a.x*b.y +  a.y*b.x ; // ir + ri
      ymm0.z = a.z*b.z -  a.w*b.w ; // rr - ii
      ymm0.w = a.w*b.z +  a.z*b.w ; // ir + ri
      return ymm0;
    }
    // Complex double
    accelerator_inline double2 operator()(double2 a, double2 b){
      double2 ymm0;
      ymm0.x = a.x*b.x -  a.y*b.y ; // rr - ii
      ymm0.y = a.x*b.y +  a.y*b.x ; // ir + ri
      return ymm0;
    }
  };

  struct Mult{

    accelerator_inline void mac(float4 &a, float4 b, float4 c){
      a= a+b*c;
    }

    accelerator_inline void mac(double2 &a, double2 b, double2 c){
      a= a+b*c;
    }

    // Real float
    accelerator_inline float4 operator()(float4 a, float4 b){
      return a*b;
    }
    // Real double
    accelerator_inline double2 operator()(double2 a, double2 b){
      return a*b;
    }
    // Integer
    accelerator_inline int4 operator()(int4 a, int4 b){
      return a*b;
    }
  };

  struct Div{
    // Real float
    accelerator_inline float4 operator()(float4 a, float4 b){
      return a/b;
    }
    // Real double
    accelerator_inline double2 operator()(double2 a, double2 b){
      return a/b;
    }
  };


  struct Conj{
    // Complex single
    accelerator_inline float4 operator()(float4 in){
      float4 ret;
      ret.x =   in.x;
      ret.y = - in.y;
      ret.z =   in.z;
      ret.w = -  in.w;
      return ret;
    }
    // Complex double
    accelerator_inline double2 operator()(double2 in){
      double2 ret;
      ret.x =   in.x;
      ret.y = - in.y;
      return ret;
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    accelerator_inline float4 operator()(float4 in, float4 ret){
      float4 tmp;
      tmp.x =   in.y;
      tmp.y = - in.x;
      tmp.z =   in.w;
      tmp.w = - in.z;
      return tmp;
    }
    //Complex double
    accelerator_inline double2 operator()(double2 in, double2 ret){
      double2 tmp;
      tmp.x =   in.y;
      tmp.y = - in.x;
      return tmp;
    }
  };

  struct TimesI{
    //Complex single
    accelerator_inline float4 operator()(float4 in, float4 ret){
      float4 tmp;
      tmp.x = - in.y;
      tmp.y =   in.x;
      tmp.z = - in.w;
      tmp.w =   in.z;
      return tmp;
    }
    //Complex double
    accelerator_inline double2 operator()(double2 in, double2 ret){
      double2 tmp ;
      tmp.x = - in.y;
      tmp.y =   in.x;
      return tmp;
    }
  };

  struct Permute{

    static accelerator_inline float4 Permute0(float4 in){
      float4 tmp;
      tmp.x = in.z;
      tmp.y = in.w;
      tmp.z = in.x;
      tmp.w = in.y;
      return tmp;
    };
    static accelerator_inline float4 Permute1(float4 in){
      float4 tmp;
      tmp.x = in.y;
      tmp.y = in.x;
      tmp.z = in.w;
      tmp.w = in.z;
      return tmp;
    };
    static accelerator_inline float4 Permute2(float4 in){
      return in;
    };
    static accelerator_inline float4 Permute3(float4 in){
      return in;
    };

    static accelerator_inline double2 Permute0(double2 in){ //AB -> BA
      double2 tmp;
      tmp.x = in.y;
      tmp.y = in.x;
      return tmp;
    };
    static accelerator_inline double2 Permute1(double2 in){
      return in;
    };
    static accelerator_inline double2 Permute2(double2 in){
      return in;
    };
    static accelerator_inline double2 Permute3(double2 in){
      return in;
    };
  };

  struct PrecisionChange {
    static accelerator_inline half8 StoH (float4 a,float4 b) {
      half8 h;
      h.ax = float2half(a.x);
      h.ay = float2half(a.y);
      h.az = float2half(a.z);
      h.aw = float2half(a.w);
      h.bx = float2half(b.x);
      h.by = float2half(b.y);
      h.bz = float2half(b.z);
      h.bw = float2half(b.w);
      return h;
    }
    static accelerator_inline void  HtoS (half8 h,float4 &sa,float4 &sb) {
      sa.x = half2float(h.ax);
      sa.y = half2float(h.ay);
      sa.z = half2float(h.az);
      sa.w = half2float(h.aw);
      sb.x = half2float(h.bx);
      sb.y = half2float(h.by);
      sb.z = half2float(h.bz);
      sb.w = half2float(h.bw);
    }
    static accelerator_inline float4 DtoS (double2 a,double2 b) {
      float4 s;
      s.x = a.x;
      s.y = a.y;
      s.z = b.x;
      s.w = b.y;
      return s;
    }
    static accelerator_inline void StoD (float4 s,double2 &a,double2 &b) {
      a.x = s.x;
      a.y = s.y;
      b.x = s.z;
      b.y = s.w;
    }
    static accelerator_inline half8 DtoH (double2 a,double2 b,double2 c,double2 d) {
      float4 sa,sb;
      sa = DtoS(a,b);
      sb = DtoS(c,d);
      return StoH(sa,sb);
    }
    static accelerator_inline void HtoD (half8 h,double2 &a,double2 &b,double2 &c,double2 &d) {
      float4 sa,sb;
      HtoS(h,sa,sb);
      StoD(sa,a,b);
      StoD(sb,c,d);
    }
  };

  struct Exchange{
    // 3210 ordering

    static accelerator_inline void Exchange0(float4 &out1,float4 &out2,float4 in1,float4 in2){
      out1.x = in1.x;      
      out1.y = in1.y;
      out1.z = in2.x;      
      out1.w = in2.y;

      out2.x = in1.z;      
      out2.y = in1.w;
      out2.z = in2.z;      
      out2.w = in2.w;

	//      out1= _mm_shuffle_ps(in1,in2,_MM_SELECT_FOUR_FOUR(1,0,1,0));
	//      out2= _mm_shuffle_ps(in1,in2,_MM_SELECT_FOUR_FOUR(3,2,3,2));
      return;
    };

    static accelerator_inline void Exchange1(float4 &out1,float4 &out2,float4 in1,float4 in2){

      out1.x = in1.x;      
      out1.y = in2.x;
      out1.z = in1.z;      
      out1.w = in2.z;

      out2.x = in1.y;      
      out2.y = in2.y;
      out2.z = in1.w;      
      out2.w = in2.w;

      //      out1= _mm_shuffle_ps(in1,in2,_MM_SELECT_FOUR_FOUR(2,0,2,0)); /*ACEG*/
      //      out2= _mm_shuffle_ps(in1,in2,_MM_SELECT_FOUR_FOUR(3,1,3,1)); /*BDFH*/
      //      out1= _mm_shuffle_ps(out1,out1,_MM_SELECT_FOUR_FOUR(3,1,2,0)); /*AECG*/
      //      out2= _mm_shuffle_ps(out2,out2,_MM_SELECT_FOUR_FOUR(3,1,2,0)); /*AECG*/
    };
    static accelerator_inline void Exchange2(float4 &out1,float4 &out2,float4 in1,float4 in2){
      assert(0);
      return;
    };
    static accelerator_inline void Exchange3(float4 &out1,float4 &out2,float4 in1,float4 in2){
      assert(0);
      return;
    };

    static accelerator_inline void Exchange0(double2 &out1,double2 &out2,double2 in1,double2 in2){
      out1.x = in1.x;      
      out1.y = in2.x;
      out2.x = in1.y;      
      out2.y = in2.y;
      //      out1= _mm_shuffle_pd(in1,in2,0x0);
      //      out2= _mm_shuffle_pd(in1,in2,0x3);
    };
    static accelerator_inline void Exchange1(double2 &out1,double2 &out2,double2 in1,double2 in2){
      assert(0);
      return;
    };
    static accelerator_inline void Exchange2(double2 &out1,double2 &out2,double2 in1,double2 in2){
      assert(0);
      return;
    };
    static accelerator_inline void Exchange3(double2 &out1,double2 &out2,double2 in1,double2 in2){
      assert(0);
      return;
    };
  };

  struct Rotate{

    static accelerator_inline float4 rotate(float4 in,int n){ 
      float4 ret;
      switch(n){
      case 0: ret = in ; break;
      case 1: ret.x = in.y; ret.y = in.z ; ret.z = in.w ; ret.w = in.x; break;
      case 2: ret.x = in.z; ret.y = in.w ; ret.z = in.x ; ret.w = in.y; break;
      case 3: ret.x = in.w; ret.y = in.x ; ret.z = in.y ; ret.w = in.z; break;
      default: break;
      }
      return ret;
    }
    static accelerator_inline double2 rotate(double2 in,int n){ 
      double2 ret;
      switch(n){
      case 0: ret = in; break;
      case 1: ret.x = in.y; ret.y = in.x ; break;
      default: break;
      }
      return ret;
    }

    template<int n> static accelerator_inline float4  tRotate(float4  in){ 
      return rotate(in,n);
    };
    template<int n> static accelerator_inline double2 tRotate(double2 in){ 
      return rotate(in,n);
    };

  };
  //////////////////////////////////////////////
  // Some Template specialization

  //Complex float Reduce
  template<>
  accelerator_inline Grid::ComplexF Reduce<Grid::ComplexF, float4>::operator()(float4 in){
    Grid::ComplexF ret(in.x+in.z,in.y+in.w);
    return ret;
  }
  //Real float Reduce
  template<>
  accelerator_inline Grid::RealF Reduce<Grid::RealF, float4>::operator()(float4 in){
    return in.x+in.y+in.z+in.w;
  }
  
  //Complex double Reduce
  template<>
  accelerator_inline Grid::ComplexD Reduce<Grid::ComplexD, double2>::operator()(double2 in){
    return Grid::ComplexD(in.x,in.y);
  }
  
  //Real double Reduce
  template<>
  accelerator_inline Grid::RealD Reduce<Grid::RealD, double2>::operator()(double2 in){
    return in.x+in.y;
  }

  //Integer Reduce
  template<>
  accelerator_inline Integer Reduce<Integer, int4>::operator()(int4 in){
    return in.x+in.y+in.z+in.w;
  }
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 
//////////////////////////////////////////////////////////////////////////////////////
  typedef half8   SIMD_Htype;  // Single precision type
  typedef float4  SIMD_Ftype;  // Single precision type
  typedef double2 SIMD_Dtype; // Double precision type
  typedef int4    SIMD_Itype; // Integer type

  // prefetch utilities
  accelerator_inline void v_prefetch0(int size, const char *ptr){};
  accelerator_inline void prefetch_HINT_T0(const char *ptr){};

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
