    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_gpu.h

    Copyright (C) 2021

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
/*! @file Grid_gpu_rrii.h*/
//----------------------------------------------------------------------

//////////////////////////////
// fp16
//////////////////////////////
#ifdef GRID_CUDA
#include <cuda_fp16.h>
#endif
#ifdef GRID_HIP
#include <hip/hip_fp16.h>
#endif
#if !defined(GRID_HIP) && !defined(GRID_CUDA) 
namespace Grid {
  typedef struct { uint16_t x;} half;
}
#endif
namespace Grid {
  accelerator_inline float half2float(half h)
  {
    float f;
#if defined(GRID_CUDA) || defined(GRID_HIP)
    f = __half2float(h);
#else 
    Grid_half hh; 
    hh.x = h.x;
    f=  sfw_half_to_float(hh);
#endif
    return f;
  }
  accelerator_inline half float2half(float f)
  {
    half h;
#if defined(GRID_CUDA) || defined(GRID_HIP)
    h = __float2half(f);
#else
    Grid_half hh = sfw_float_to_half(f);
    h.x = hh.x;
#endif
    return h;
  }
}


#define COALESCE_GRANULARITY ( GEN_SIMD_WIDTH )

namespace Grid {

////////////////////////////////////////////////////////////////////////
// Real vector
////////////////////////////////////////////////////////////////////////  
template<int _N, class _datum>
struct GpuVector {
  _datum rrrr[_N];
  static const int N = _N;
  typedef _datum datum;
};
template<int N,class datum>
inline accelerator GpuVector<N,datum> operator*(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.rrrr[i] = l.rrrr[i]*r.rrrr[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuVector<N,datum> operator-(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.rrrr[i] = l.rrrr[i]-r.rrrr[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuVector<N,datum> operator+(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.rrrr[i] = l.rrrr[i]+r.rrrr[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuVector<N,datum> operator/(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.rrrr[i] = l.rrrr[i]/r.rrrr[i];
  }
  return ret;
}

////////////////////////////////////////////////////////////////////////
// Complex vector
////////////////////////////////////////////////////////////////////////  
template<int _N, class _datum>
struct GpuComplexVector {
  _datum rrrr[_N];
  _datum iiii[_N];
  static const int N = _N;
  typedef _datum datum;
};
template<int N,class datum>
inline accelerator GpuComplexVector<N,datum> operator*(const GpuComplexVector<N,datum> l,const GpuComplexVector<N,datum> r) {
  GpuComplexVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.rrrr[i] = l.rrrr[i]*r.rrrr[i] - l.iiii[i]*r.iiii[i];
    ret.iiii[i] = l.rrrr[i]*r.iiii[i] + l.iiii[i]*r.rrrr[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuComplexVector<N,datum> operator-(const GpuComplexVector<N,datum> l,const GpuComplexVector<N,datum> r) {
  GpuComplexVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.rrrr[i] = l.rrrr[i]-r.rrrr[i];
    ret.iiii[i] = l.iiii[i]-r.iiii[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuComplexVector<N,datum> operator+(const GpuComplexVector<N,datum> l,const GpuComplexVector<N,datum> r) {
  GpuComplexVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.rrrr[i] = l.rrrr[i]+r.rrrr[i];
    ret.iiii[i] = l.iiii[i]+r.iiii[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuComplexVector<N,datum> operator/(const GpuComplexVector<N,datum> l,const GpuComplexVector<N,datum> r) {
  GpuComplexVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.rrrr[i] = l.rrrr[i]/r.rrrr[i];
    ret.iiii[i] = l.iiii[i]/r.iiii[i];
  }
  return ret;
}

////////////////////////////////
// SIMD counts
////////////////////////////////

constexpr int NSIMD_RealH    = COALESCE_GRANULARITY / sizeof(half);
constexpr int NSIMD_ComplexH = COALESCE_GRANULARITY / sizeof(half);
constexpr int NSIMD_RealF    = COALESCE_GRANULARITY / sizeof(float);
constexpr int NSIMD_ComplexF = COALESCE_GRANULARITY / sizeof(float);
constexpr int NSIMD_RealD    = COALESCE_GRANULARITY / sizeof(double);
constexpr int NSIMD_ComplexD = COALESCE_GRANULARITY / sizeof(double);
constexpr int NSIMD_Integer  = COALESCE_GRANULARITY / sizeof(Integer);

typedef GpuVector<NSIMD_RealH   , half        > GpuVectorRH;
typedef GpuComplexVector<NSIMD_ComplexH, half > GpuVectorCH;
typedef GpuVector<NSIMD_RealF,    float       > GpuVectorRF;
typedef GpuComplexVector<NSIMD_ComplexF, float> GpuVectorCF;
typedef GpuVector<NSIMD_RealD,    double      > GpuVectorRD;
typedef GpuComplexVector<NSIMD_ComplexD,double> GpuVectorCD;
typedef GpuVector<NSIMD_Integer,  Integer     > GpuVectorI;

namespace Optimization {

  struct Vsplat{
    //Complex float
    accelerator_inline GpuVectorCF operator()(float a, float b){
      GpuVectorCF ret;
      for(int i=0;i<GpuVectorCF::N;i++){
	ret.rrrr[i] = typename GpuVectorCF::datum(a);
	ret.iiii[i] = typename GpuVectorCF::datum(b);
      }
      return ret;
    }
    // Real float
    accelerator_inline GpuVectorRF operator()(float a){
      GpuVectorRF ret;
      for(int i=0;i<GpuVectorRF::N;i++){
	ret.rrrr[i] = typename GpuVectorRF::datum(a);
      }
      return ret;
    }
    //Complex double
    accelerator_inline GpuVectorCD operator()(double a, double b){
      GpuVectorCD ret;
      for(int i=0;i<GpuVectorCD::N;i++){
	ret.rrrr[i] = typename GpuVectorCD::datum(a);
	ret.iiii[i] = typename GpuVectorCD::datum(b);
      }
      return ret;
    }
    //Real double
    accelerator_inline GpuVectorRD operator()(double a){
      GpuVectorRD ret; 
      for(int i=0;i<GpuVectorRD::N;i++){
	ret.rrrr[i] = typename GpuVectorRD::datum(a);
      }
      return ret;
    }
    //Integer
    accelerator_inline GpuVectorI operator()(Integer a){
      GpuVectorI ret;
      for(int i=0;i<GpuVectorI::N;i++){
	ret.rrrr[i] = typename GpuVectorI::datum(a);
      }
      return ret;
    }
  };

  struct Vstore{
    template<int N,class datum,class P>
    accelerator_inline void operator()(GpuVector<N,datum> a, P* Fp){
      GpuVector<N,datum> *vF = (GpuVector<N,datum> *)Fp;
      *vF = a;
    }
    template<int N,class datum,class P>
    accelerator_inline void operator()(GpuComplexVector<N,datum> a, P* Fp){
      GpuComplexVector<N,datum> *vF = (GpuComplexVector<N,datum> *)Fp;
      *vF = a;
    }
  };

  struct Vstream{
    template<int N,class datum, class P>
    accelerator_inline void operator()(P* F,GpuVector<N,datum> a){
      GpuVector<N,datum> *vF = (GpuVector<N,datum> *)F;
      *vF = a;
    }
    template<int N,class datum, class P>
    accelerator_inline void operator()(P* F,GpuComplexVector<N,datum> a){
      GpuComplexVector<N,datum> *vF = (GpuComplexVector<N,datum> *)F;
      *vF = a;
    }
  };

  struct Vset{
    // Complex float 
    accelerator_inline GpuVectorCF operator()(Grid::ComplexF *a){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = vec::datum(a[i].real());
	ret.iiii[i] = vec::datum(a[i].imag());
      }
      return ret;
    }
    // Complex double 
    accelerator_inline GpuVectorCD operator()(Grid::ComplexD *a){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = vec::datum(a[i].real());
	ret.iiii[i] = vec::datum(a[i].imag());
      }
      return ret;
    }
    // Real float 
    accelerator_inline GpuVectorRF operator()(float *a){
      typedef GpuVectorRF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = vec::datum(a[i]);
      }
      return ret;
    }
    // Real double
    accelerator_inline GpuVectorRD operator()(double *a){
      typedef GpuVectorRD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = vec::datum(a[i]);
      }
      return ret;
    }
    // Integer
    accelerator_inline GpuVectorI operator()(Integer *a){
      typedef GpuVectorI vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = vec::datum(a[i]);
      }
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
    //Real float
    accelerator_inline GpuVectorRF operator()(GpuVectorRF a,GpuVectorRF b){
      return a+b;
    }
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a,GpuVectorRD b){
      return a+b;
    }
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      return a+b;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      return a+b;
    }
    accelerator_inline GpuVectorI operator()(GpuVectorI a,GpuVectorI b){
      return a+b;
    }
  };

  struct Sub{
    accelerator_inline GpuVectorRF operator()(GpuVectorRF a,GpuVectorRF b){
      return a-b;
    }
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a,GpuVectorRD b){
      return a-b;
    }
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      return a-b;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      return a-b;
    }
    accelerator_inline GpuVectorI operator()(GpuVectorI a,GpuVectorI b){
      return a-b;
    }
  };

  struct MultRealPart{
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = a.rrrr[i]*b.rrrr[i];
	ret.iiii[i] = a.rrrr[i]*b.iiii[i];
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = a.rrrr[i]*b.rrrr[i];
	ret.iiii[i] = a.rrrr[i]*b.iiii[i];
      }
      return ret;
    }
  };

  struct MaddRealPart{
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b,GpuVectorCF c){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = a.rrrr[i]*b.rrrr[i]+c.rrrr[i];
	ret.iiii[i] = a.rrrr[i]*b.iiii[i]+c.iiii[i];
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b,GpuVectorCD c){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = a.rrrr[i]*b.rrrr[i]+c.rrrr[i];
	ret.iiii[i] = a.rrrr[i]*b.iiii[i]+c.iiii[i];
      }
      return ret;
    }
  };

  struct MultComplex{

    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b){
      return a*b;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      return a*b;
    }
  };

  struct Mult{
    accelerator_inline void mac(GpuVectorRF &a, GpuVectorRF b, GpuVectorRF c){
      a= a+b*c;
    }
    accelerator_inline void mac(GpuVectorRD &a, GpuVectorRD b, GpuVectorRD c){
      a= a+b*c;
    }
    // Real float
    accelerator_inline GpuVectorRF operator()(GpuVectorRF a, GpuVectorRF b){
      return a*b;
    }
    // Real double
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a, GpuVectorRD b){
      return a*b;
    }
    accelerator_inline GpuVectorI operator()(GpuVectorI a, GpuVectorI b){
      return a*b;
    }
  };

  struct Div{
    // Real float
    accelerator_inline GpuVectorRF operator()(GpuVectorRF a, GpuVectorRF b){
      return a/b;
    }
    accelerator_inline GpuVectorRD operator()(GpuVectorRD a, GpuVectorRD b){
      return a/b;
    }
    accelerator_inline GpuVectorI operator()(GpuVectorI a, GpuVectorI b){
      return a/b;
    }

    // Danger -- element wise divide fro complex, not complex div. 
    // See Grid_vector_types.h lines around 735, applied after "toReal"
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a, GpuVectorCF b){
      return a/b;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a, GpuVectorCD b){
      return a/b;
    }
  };


  struct Conj{
    // Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = in.rrrr[i];
	ret.iiii[i] =-in.iiii[i];
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = in.rrrr[i];
	ret.iiii[i] =-in.iiii[i];
      }
      return ret;
    }
  };

  struct TimesMinusI{
    //Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = in.iiii[i];
	ret.iiii[i] =-in.rrrr[i];
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] = in.iiii[i];
	ret.iiii[i] =-in.rrrr[i];
      }
      return ret;
    }
  };

  struct TimesI{
    //Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] =-in.iiii[i];
	ret.iiii[i] = in.rrrr[i];
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.rrrr[i] =-in.iiii[i];
	ret.iiii[i] = in.rrrr[i];
      }
      return ret;
    }
  };

  struct Permute{

    template <int n,int _N, class _datum >
    static accelerator_inline GpuVector<_N,_datum> PermuteN(GpuVector<_N,_datum> &in) {   
      typedef GpuVector<_N,_datum> vec;
      vec out;					
      unsigned int _mask = vec::N >> (n + 1);	
      for(int i=0;i<vec::N;i++) {
	out.rrrr[i] = in.rrrr[i^_mask];
      }
      return out;	
    }
    template <int n,int _N, class _datum >
    static accelerator_inline GpuComplexVector<_N,_datum> PermuteN(GpuComplexVector<_N,_datum> &in) {   
      typedef GpuComplexVector<_N,_datum> vec;
      vec out;					
      unsigned int _mask = vec::N >> (n + 1);	
      for(int i=0;i<vec::N;i++) {
	out.rrrr[i] = in.rrrr[i^_mask];
	out.iiii[i] = in.iiii[i^_mask];
      }
      return out;	
    }
    
    template <typename vec>  static accelerator_inline vec Permute0(vec in) { return PermuteN<0,vec::N,typename vec::datum>(in);  }
    template <typename vec>  static accelerator_inline vec Permute1(vec in) { return PermuteN<1,vec::N,typename vec::datum>(in);  }
    template <typename vec>  static accelerator_inline vec Permute2(vec in) { return PermuteN<2,vec::N,typename vec::datum>(in);  }
    template <typename vec>  static accelerator_inline vec Permute3(vec in) { return PermuteN<3,vec::N,typename vec::datum>(in);  }
    
  };
  
  struct PrecisionChange {

    ////////////////////////////////////////////////////////////////////////////////////
    // Single / Half
    ////////////////////////////////////////////////////////////////////////////////////
     static accelerator_inline GpuVectorCH StoH (GpuVectorCF a,GpuVectorCF b) {
      int N = GpuVectorCF::N;
      GpuVectorCH h;
      for(int i=0;i<N;i++) {
        h.rrrr[i  ] = float2half(a.rrrr[i]);
        h.iiii[i  ] = float2half(a.iiii[i]);
	h.rrrr[i+N] = float2half(b.rrrr[i]);
	h.iiii[i+N] = float2half(b.iiii[i]);
      }
      return h;
    }
    static accelerator_inline void  HtoS (GpuVectorCH h,GpuVectorCF &sa,GpuVectorCF &sb) {
      int N = GpuVectorCF::N;
      for(int i=0;i<N;i++) {
	sa.rrrr[i] = half2float(h.rrrr[i  ]);
	sa.iiii[i] = half2float(h.iiii[i  ]);
	sb.rrrr[i] = half2float(h.rrrr[i+N]);
	sb.iiii[i] = half2float(h.iiii[i+N]);
      }
    }
    static accelerator_inline GpuVectorRH StoH (GpuVectorRF a,GpuVectorRF b) {
      int N = GpuVectorRF::N;
      GpuVectorRH h;
      for(int i=0;i<N;i++) {
        h.rrrr[i  ] = float2half(a.rrrr[i]);
	h.rrrr[i+N] = float2half(b.rrrr[i]);
      }
      return h;
    }
    static accelerator_inline void  HtoS (GpuVectorRH h,GpuVectorRF &sa,GpuVectorRF &sb) {
      int N = GpuVectorRF::N;
      for(int i=0;i<N;i++) {
	sa.rrrr[i] = half2float(h.rrrr[i  ]);
	sb.rrrr[i] = half2float(h.rrrr[i+N]);
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // Double Single
    ////////////////////////////////////////////////////////////////////////////////////
    static accelerator_inline GpuVectorCF DtoS (GpuVectorCD a,GpuVectorCD b) {
      int N = GpuVectorCD::N;
      GpuVectorCF h;
      for(int i=0;i<N;i++) {
        h.rrrr[i  ] = a.rrrr[i];
        h.iiii[i  ] = a.iiii[i];
	h.rrrr[i+N] = b.rrrr[i];
	h.iiii[i+N] = b.iiii[i];
      }
      return h;
    }

    static accelerator_inline void  StoD (GpuVectorCF h,GpuVectorCD &sa,GpuVectorCD &sb) {
      int N = GpuVectorCD::N;
      for(int i=0;i<N;i++) {
	sa.rrrr[i] = h.rrrr[i  ];
	sa.iiii[i] = h.iiii[i  ];
	sb.rrrr[i] = h.rrrr[i+N];
	sb.iiii[i] = h.iiii[i+N];
      }
    }

    static accelerator_inline GpuVectorRF DtoS (GpuVectorRD a,GpuVectorRD b) {
      int N = GpuVectorRD::N;
      GpuVectorRF h;
      for(int i=0;i<N;i++) {
        h.rrrr[i  ] = a.rrrr[i];
	h.rrrr[i+N] = b.rrrr[i];
      }
      return h;
    }

    static accelerator_inline void  StoD (GpuVectorRF h,GpuVectorRD &sa,GpuVectorRD &sb) {
      int N = GpuVectorRD::N;
      for(int i=0;i<N;i++) {
	sa.rrrr[i] = h.rrrr[i  ];
	sb.rrrr[i] = h.rrrr[i+N];
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // Double Half
    ////////////////////////////////////////////////////////////////////////////////////
    static accelerator_inline GpuVectorCH DtoH (GpuVectorCD a,GpuVectorCD b,GpuVectorCD c,GpuVectorCD d) {
      GpuVectorCF sa,sb;
      sa = DtoS(a,b);
      sb = DtoS(c,d);
      return StoH(sa,sb);
    }
    static accelerator_inline void HtoD (GpuVectorCH h,GpuVectorCD &a,GpuVectorCD &b,GpuVectorCD &c,GpuVectorCD &d) {
      GpuVectorCF sa,sb;
      HtoS(h,sa,sb);
      StoD(sa,a,b);
      StoD(sb,c,d);
    }
    static accelerator_inline GpuVectorRH DtoH (GpuVectorRD a,GpuVectorRD b,GpuVectorRD c,GpuVectorRD d) {
      GpuVectorRF sa,sb;
      sa = DtoS(a,b);
      sb = DtoS(c,d);
      return StoH(sa,sb);
    }
    static accelerator_inline void HtoD (GpuVectorRH h,GpuVectorRD &a,GpuVectorRD &b,GpuVectorRD &c,GpuVectorRD &d) {
      GpuVectorRF sa,sb;
      HtoS(h,sa,sb);
      StoD(sa,a,b);
      StoD(sb,c,d);
    }
  };

struct Exchange{

  template <int n,int _N, class _datum >
  static accelerator_inline void ExchangeN(GpuVector<_N,_datum> &out1,
					   GpuVector<_N,_datum> &out2,
					   GpuVector<_N,_datum> &in1,
					   GpuVector<_N,_datum> &in2 )
  {   
    typedef GpuVector<_N,_datum> vec;
    unsigned int mask = vec::N >> (n + 1);
    for(int i=0;i<vec::N;i++) {
      int j1 = i&(~mask);
      if  ( (i&mask) == 0 ) { out1.rrrr[i]=in1.rrrr[j1];}
      else                  { out1.rrrr[i]=in2.rrrr[j1];}
      int j2 = i|mask;
      if  ( (i&mask) == 0 ) { out2.rrrr[i]=in1.rrrr[j2];}
      else                  { out2.rrrr[i]=in2.rrrr[j2];}
    }      
  }
  template <int n,int _N, class _datum >
  static accelerator_inline void ExchangeN(GpuComplexVector<_N,_datum> &out1,
					   GpuComplexVector<_N,_datum> &out2,
					   GpuComplexVector<_N,_datum> &in1,
					   GpuComplexVector<_N,_datum> &in2 )
  {   
    typedef GpuComplexVector<_N,_datum> vec;
    unsigned int mask = vec::N >> (n + 1);
    for(int i=0;i<vec::N;i++) {
      int j1 = i&(~mask);
      if  ( (i&mask) == 0 ) {
	out1.rrrr[i]=in1.rrrr[j1];
	out1.iiii[i]=in1.iiii[j1];
      }
      else                  {
	out1.rrrr[i]=in2.rrrr[j1];
	out1.iiii[i]=in2.iiii[j1];
      }
      int j2 = i|mask;
      if  ( (i&mask) == 0 ) {
	out2.rrrr[i]=in1.rrrr[j2];
	out2.iiii[i]=in1.iiii[j2];
      }
      else                  {
	out2.rrrr[i]=in2.rrrr[j2];
	out2.iiii[i]=in2.iiii[j2];
      }
    }      
  }
  template <typename vec>
  static accelerator_inline void Exchange0(vec &out1,vec &out2,vec &in1,vec &in2){
    ExchangeN<0>(out1,out2,in1,in2);
  };
  template <typename vec>
  static accelerator_inline void Exchange1(vec &out1,vec &out2,vec &in1,vec &in2){
    ExchangeN<1>(out1,out2,in1,in2);
  };
  template <typename vec>
  static accelerator_inline void Exchange2(vec &out1,vec &out2,vec &in1,vec &in2){
    ExchangeN<2>(out1,out2,in1,in2);
  };
  template <typename vec>
  static accelerator_inline void Exchange3(vec &out1,vec &out2,vec &in1,vec &in2){
    ExchangeN<3>(out1,out2,in1,in2);
  };

};

struct Rotate{

  template <int n, typename vec> static accelerator_inline vec tRotate(vec in){
    return rotate(in, n);
  }
    
  template <int _N, class _datum >
  static accelerator_inline GpuComplexVector<_N,_datum> rotate_template(GpuComplexVector<_N,_datum> &in, int n)
  {
    typedef GpuComplexVector<_N,_datum> vec;
    vec out;
    for(int i=0;i<vec::N;i++){
      out.rrrr[i] = in.rrrr[(i + n)%vec::N];
      out.iiii[i] = in.iiii[(i + n)%vec::N];
    }
    return out;
  }

  template <int _N, class _datum >
  static accelerator_inline GpuVector<_N,_datum> rotate_template(GpuVector<_N,_datum> &in, int n)
  {
    typedef GpuVector<_N,_datum> vec;
    vec out;
    for(int i=0;i<vec::N;i++){
      out.rrrr[i] = in.rrrr[(i + n)%vec::N];
    }
    return out;
  }

  typedef GpuVectorRH  SIMD_Htype; // Single precision type
  typedef GpuVectorRF  SIMD_Ftype; // Single precision type
  typedef GpuVectorRD  SIMD_Dtype; // Double precision type
  typedef GpuVectorI   SIMD_Itype; // Integer type

  typedef GpuVectorCH  SIMD_CHtype; // Single precision type
  typedef GpuVectorCF  SIMD_CFtype; // Single precision type
  typedef GpuVectorCD  SIMD_CDtype; // Double precision type

  static accelerator_inline GpuVectorRH rotate(GpuVectorRH in, int n){ return rotate_template(in,n);}
  static accelerator_inline GpuVectorRF rotate(GpuVectorRF in, int n){ return rotate_template(in,n);}
  static accelerator_inline GpuVectorRD rotate(GpuVectorRD in, int n){ return rotate_template(in,n);}
  static accelerator_inline GpuVectorI  rotate(GpuVectorI  in, int n){ return rotate_template(in,n);}
  static accelerator_inline GpuVectorCH rotate(GpuVectorCH in, int n){ return rotate_template(in,n/2);} // Measure in complex not float
  static accelerator_inline GpuVectorCF rotate(GpuVectorCF in, int n){ return rotate_template(in,n/2);}
  static accelerator_inline GpuVectorCD rotate(GpuVectorCD in, int n){ return rotate_template(in,n/2);}

};

//////////////////////////////////////////////
// Some Template specialization

  //Complex float Reduce
  template<>
  accelerator_inline Grid::ComplexF 
  Reduce<Grid::ComplexF, GpuVectorCF>::operator()(GpuVectorCF in)
  {
    Grid::ComplexF greduce(in.rrrr[0],in.iiii[0]);
    for(int i=1;i<GpuVectorCF::N;i++) {
      greduce = greduce+Grid::ComplexF(in.rrrr[i],in.iiii[i]);
    }
    return greduce;
  }

  template<>
  accelerator_inline Grid::ComplexD
  Reduce<Grid::ComplexD, GpuVectorCD>::operator()(GpuVectorCD in)
  {
    Grid::ComplexD greduce(in.rrrr[0],in.iiii[0]);
    for(int i=1;i<GpuVectorCD::N;i++) {
      greduce = greduce+ Grid::ComplexD(in.rrrr[i],in.iiii[i]);
    }
    return greduce;
  }

  // Real
  template<>
  accelerator_inline Grid::RealF 
  Reduce<RealF, GpuVectorRF>::operator()(GpuVectorRF in)
  {
    RealF ret = in.rrrr[0];
    for(int i=1;i<GpuVectorRF::N;i++) {
      ret = ret+in.rrrr[i];
    }
    return ret;
  }

  template<>
  accelerator_inline Grid::RealD 
  Reduce<RealD, GpuVectorRD>::operator()(GpuVectorRD in)
  {
    RealD ret = in.rrrr[0];
    for(int i=1;i<GpuVectorRD::N;i++) {
      ret = ret+in.rrrr[i];
    }
    return ret;
  }

  template<>
  accelerator_inline Integer
  Reduce<Integer, GpuVectorI>::operator()(GpuVectorI in)
  {
    Integer ret = in.rrrr[0];
    for(int i=1;i<GpuVectorI::N;i++) {
      ret = ret+in.rrrr[i];
    }
    return ret;
  }

}// End optimizatoin

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 
//////////////////////////////////////////////////////////////////////////////////////
  typedef GpuVectorRH  SIMD_Htype; // Single precision type
  typedef GpuVectorRF  SIMD_Ftype; // Single precision type
  typedef GpuVectorRD  SIMD_Dtype; // Double precision type
  typedef GpuVectorI   SIMD_Itype; // Integer type

  typedef GpuVectorCH  SIMD_CHtype; // Single precision type
  typedef GpuVectorCF  SIMD_CFtype; // Single precision type
  typedef GpuVectorCD  SIMD_CDtype; // Double precision type

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
