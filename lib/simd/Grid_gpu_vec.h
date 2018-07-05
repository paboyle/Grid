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

#define COALESCE_GRANULARITY (16) // bytes for coalesce granularity of target: Pascal, Volta

template<class pair>
class GpuComplex {
public:
  pair z;
  typedef decltype(z.x) real;
public: 
  accelerator_inline GpuComplex() = default;
  accelerator_inline GpuComplex(real re,real im) { z.x=re; z.y=im; };
  accelerator_inline GpuComplex(const GpuComplex &zz) { z = zz.z;};
  friend accelerator_inline  GpuComplex operator+(const GpuComplex &lhs,const GpuComplex &rhs) { 
    GpuComplex r ; 
    r.z.x = lhs.z.x + rhs.z.x; 
    r.z.y = lhs.z.y + rhs.z.y; 
    return r; 
  }
  friend accelerator_inline GpuComplex operator-(const GpuComplex &lhs,const GpuComplex &rhs) { 
    GpuComplex r ; 
    r.z.x = lhs.z.x - rhs.z.x; 
    r.z.y = lhs.z.y - rhs.z.y; 
    return r; 
  }
  friend accelerator_inline GpuComplex operator*(const GpuComplex &lhs,const GpuComplex &rhs) { 
    GpuComplex r ; 
    r.z.x= lhs.z.x*rhs.z.x - lhs.z.y*rhs.z.y; // rr-ii
    r.z.y= lhs.z.x*rhs.z.y + lhs.z.y*rhs.z.x; // ri+ir
    return r;
  }
  friend accelerator_inline GpuComplex real_mult(const GpuComplex &l,const GpuComplex &r) 
  {
    GpuComplex ret;
    ret.z.x = l.z.x*r.z.x;
    ret.z.y = l.z.x*r.z.y;
    return ret;
  }
  friend std::ostream& operator<< (std::ostream& stream, const GpuComplex o){
    stream << "("<< o.z.x << ","<< o.z.y <<")";
    return stream;
  }
};

template<int _N, class _datum>
struct GpuVector {
  _datum v[_N];
  static const int N = _N;
  typedef _datum datum;
};


template<int N,class datum>
inline accelerator GpuVector<N,datum> operator*(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.v[i] = l.v[i]*r.v[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuVector<N,datum> operator-(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.v[i] = l.v[i]-r.v[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuVector<N,datum> operator+(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.v[i] = l.v[i]+r.v[i];
  }
  return ret;
}
template<int N,class datum>
inline accelerator GpuVector<N,datum> operator/(const GpuVector<N,datum> l,const GpuVector<N,datum> r) {
  GpuVector<N,datum> ret;
  for(int i=0;i<N;i++) { 
    ret.v[i] = l.v[i]/r.v[i];
  }
  return ret;
}

constexpr int NSIMD_RealH    = COALESCE_GRANULARITY / sizeof(half);
constexpr int NSIMD_ComplexH = COALESCE_GRANULARITY / sizeof(half2);
constexpr int NSIMD_RealF    = COALESCE_GRANULARITY / sizeof(float);
constexpr int NSIMD_ComplexF = COALESCE_GRANULARITY / sizeof(float2);
constexpr int NSIMD_RealD    = COALESCE_GRANULARITY / sizeof(double);
constexpr int NSIMD_ComplexD = COALESCE_GRANULARITY / sizeof(double2);
constexpr int NSIMD_Integer  = COALESCE_GRANULARITY / sizeof(Integer);

typedef GpuComplex<half2  > GpuComplexH;
typedef GpuComplex<float2 > GpuComplexF;
typedef GpuComplex<double2> GpuComplexD;

typedef GpuVector<NSIMD_RealH   , half        > GpuVectorRH;
typedef GpuVector<NSIMD_ComplexH, GpuComplexH > GpuVectorCH;
typedef GpuVector<NSIMD_RealF,    float       > GpuVectorRF;
typedef GpuVector<NSIMD_ComplexF, GpuComplexF > GpuVectorCF;
typedef GpuVector<NSIMD_RealD,    double      > GpuVectorRD;
typedef GpuVector<NSIMD_ComplexD, GpuComplexD > GpuVectorCD;
typedef GpuVector<NSIMD_Integer,  Integer     > GpuVectorI;

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

  struct Vsplat{
    //Complex float
    accelerator_inline GpuVectorCF operator()(float a, float b){
      GpuVectorCF ret;
      for(int i=0;i<GpuVectorCF::N;i++){
	ret.v[i] = typename GpuVectorCF::datum(a,b);
      }
      return ret;
    }
    // Real float
    accelerator_inline GpuVectorRF operator()(float a){
      GpuVectorRF ret;
      for(int i=0;i<GpuVectorRF::N;i++){
	ret.v[i] = typename GpuVectorRF::datum(a);
      }
      return ret;
    }
    //Complex double
    accelerator_inline GpuVectorCD operator()(double a, double b){
      GpuVectorCD ret;
      for(int i=0;i<GpuVectorCD::N;i++){
	ret.v[i] = typename GpuVectorCD::datum(a,b);
      }
      return ret;
    }
    //Real double
    accelerator_inline GpuVectorRD operator()(double a){
      GpuVectorRD ret; 
      for(int i=0;i<GpuVectorRD::N;i++){
	ret.v[i] = typename GpuVectorRD::datum(a);
      }
      return ret;
    }
    //Integer
    accelerator_inline GpuVectorI operator()(Integer a){
      GpuVectorI ret;
      for(int i=0;i<GpuVectorI::N;i++){
	ret.v[i] = typename GpuVectorI::datum(a);
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
  };

  struct Vstream{
    template<int N,class datum, class P>
    accelerator_inline void operator()(P* F,GpuVector<N,datum> a){
      GpuVector<N,datum> *vF = (GpuVector<N,datum> *)F;
      *vF = a;
    }
  };

  struct Vset{
    // Complex float 
    accelerator_inline GpuVectorCF operator()(Grid::ComplexF *a){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i].real(),a[i].imag());
      }
      return ret;
    }
    // Complex double 
    accelerator_inline GpuVectorCD operator()(Grid::ComplexD *a){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i].real(),a[i].imag());
      }
      return ret;
    }
    // Real float 
    accelerator_inline GpuVectorRF operator()(float *a){
      typedef GpuVectorRF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i]);
      }
      return ret;
    }
    // Real double
    accelerator_inline GpuVectorRD operator()(double *a){
      typedef GpuVectorRD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i]);
      }
      return ret;
    }
    // Integer
    accelerator_inline GpuVectorI operator()(Integer *a){
      typedef GpuVectorI vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = vec::datum(a[i]);
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
	ret.v[i] = real_mult(a.v[i],b.v[i]);
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = real_mult(a.v[i],b.v[i]);
      }
      return ret;
    }
  };

  struct MaddRealPart{
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a,GpuVectorCF b,GpuVectorCF c){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = real_mult(a.v[i],b.v[i]) +c.v[i];
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a,GpuVectorCD b,GpuVectorCD c){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i] = real_mult(a.v[i],b.v[i]) +c.v[i];
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
    // Danger -- element wise divide fro complex, not complex div. 
    // See Grid_vector_types.h lines around 735, applied after "toReal"
    accelerator_inline GpuVectorCF operator()(GpuVectorCF a, GpuVectorCF b){
      GpuVectorCF ret;
      for(int i=0;i< GpuVectorCF::N;i++){
	ret.v[i].z.x = a.v[i].z.x / b.v[i].z.x;
	ret.v[i].z.y = a.v[i].z.y / b.v[i].z.y;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD a, GpuVectorCD b){
      GpuVectorCD ret;
      for(int i=0;i< GpuVectorCD::N;i++){
	ret.v[i].z.x = a.v[i].z.x / b.v[i].z.x;
	ret.v[i].z.y = a.v[i].z.y / b.v[i].z.y;
      }
      return ret;
    }
  };


  struct Conj{
    // Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x = in.v[i].z.x;
	ret.v[i].z.y =-in.v[i].z.y;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x = in.v[i].z.x;
	ret.v[i].z.y =-in.v[i].z.y;
      }
      return ret;
    }
  };

  struct TimesMinusI{
    //Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in,GpuVectorCF dummy){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x = in.v[i].z.y;
	ret.v[i].z.y =-in.v[i].z.x;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in,GpuVectorCD dummy){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x = in.v[i].z.y;
	ret.v[i].z.y =-in.v[i].z.x;
      }
      return ret;
    }
  };

  struct TimesI{
    //Complex single
    accelerator_inline GpuVectorCF operator()(GpuVectorCF in,GpuVectorCF dummy){
      typedef GpuVectorCF vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x =-in.v[i].z.y;
	ret.v[i].z.y = in.v[i].z.x;
      }
      return ret;
    }
    accelerator_inline GpuVectorCD operator()(GpuVectorCD in,GpuVectorCD dummy){
      typedef GpuVectorCD vec;
      vec ret;
      for(int i=0;i<vec::N;i++){
	ret.v[i].z.x =-in.v[i].z.y;
	ret.v[i].z.y = in.v[i].z.x;
      }
      return ret;
    }
  };

  struct Permute{

    template <int n,typename vec>				       
    static accelerator_inline vec PermuteN(vec in) {   
      vec out;					
      unsigned int _mask = vec::N >> (n + 1);	
      for(int i=0;i<vec::N;i++) {
	out.v[i] = in.v[i^_mask];
      }
      return out;	
    }
    
    template <typename vec>  static accelerator_inline vec Permute0(vec in) { return PermuteN<0,vec>(in);  }
    template <typename vec>  static accelerator_inline vec Permute1(vec in) { return PermuteN<1,vec>(in);  }
    template <typename vec>  static accelerator_inline vec Permute2(vec in) { return PermuteN<2,vec>(in);  }
    template <typename vec>  static accelerator_inline vec Permute3(vec in) { return PermuteN<3,vec>(in);  }
    
  };
  
  struct PrecisionChange {

    ////////////////////////////////////////////////////////////////////////////////////
    // Single / Half
    ////////////////////////////////////////////////////////////////////////////////////
    static accelerator_inline GpuVectorCH StoH (GpuVectorCF a,GpuVectorCF b) {
      int N = GpuVectorCF::N;
      GpuVectorCH h;
      for(int i=0;i<N;i++) {
        h.v[i  ].z.x = float2half(a.v[i].z.x);
        h.v[i  ].z.y = float2half(a.v[i].z.y);
	h.v[i+N].z.x = float2half(b.v[i].z.x);
	h.v[i+N].z.y = float2half(b.v[i].z.y);
      }
      return h;
    }
    static accelerator_inline void  HtoS (GpuVectorCH h,GpuVectorCF &sa,GpuVectorCF &sb) {
      int N = GpuVectorCF::N;
      for(int i=0;i<N;i++) {
	sa.v[i].z.x = half2float(h.v[i  ].z.x);
	sa.v[i].z.y = half2float(h.v[i  ].z.y);
	sb.v[i].z.x = half2float(h.v[i+N].z.x);
	sb.v[i].z.y = half2float(h.v[i+N].z.y);
      }
    }
    static accelerator_inline GpuVectorRH StoH (GpuVectorRF a,GpuVectorRF b) {
      int N = GpuVectorRF::N;
      GpuVectorRH h;
      for(int i=0;i<N;i++) {
        h.v[i  ] = float2half(a.v[i]);
	h.v[i+N] = float2half(b.v[i]);
      }
      return h;
    }
    static accelerator_inline void  HtoS (GpuVectorRH h,GpuVectorRF &sa,GpuVectorRF &sb) {
      int N = GpuVectorRF::N;
      for(int i=0;i<N;i++) {
	sa.v[i] = half2float(h.v[i  ]);
	sb.v[i] = half2float(h.v[i+N]);
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // Double Single
    ////////////////////////////////////////////////////////////////////////////////////
    static accelerator_inline GpuVectorCF DtoS (GpuVectorCD a,GpuVectorCD b) {
      int N = GpuVectorCD::N;
      GpuVectorCF h;
      for(int i=0;i<N;i++) {
        h.v[i  ].z.x = a.v[i].z.x;
        h.v[i  ].z.y = a.v[i].z.y;
	h.v[i+N].z.x = b.v[i].z.x;
	h.v[i+N].z.y = b.v[i].z.y;
      }
      return h;
    }

    static accelerator_inline void  StoD (GpuVectorCF h,GpuVectorCD &sa,GpuVectorCD &sb) {
      int N = GpuVectorCD::N;
      for(int i=0;i<N;i++) {
	sa.v[i].z.x = h.v[i  ].z.x;
	sa.v[i].z.y = h.v[i  ].z.y;
	sb.v[i].z.x = h.v[i+N].z.x;
	sb.v[i].z.y = h.v[i+N].z.y;
      }
    }

    static accelerator_inline GpuVectorRF DtoS (GpuVectorRD a,GpuVectorRD b) {
      int N = GpuVectorRD::N;
      GpuVectorRF h;
      for(int i=0;i<N;i++) {
        h.v[i  ] = a.v[i];
	h.v[i+N] = b.v[i];
      }
      return h;
    }

    static accelerator_inline void  StoD (GpuVectorRF h,GpuVectorRD &sa,GpuVectorRD &sb) {
      int N = GpuVectorRD::N;
      for(int i=0;i<N;i++) {
	sa.v[i] = h.v[i  ];
	sb.v[i] = h.v[i+N];
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

  template <typename vec,int n>
  static accelerator_inline void ExchangeN(vec &out1,vec &out2,vec &in1,vec &in2){
    unsigned int mask = vec::N >> (n + 1);
    for(int i=0;i<vec::N;i++) {
      int j1 = i&(~mask);
      if  ( (i&mask) == 0 ) { out1.v[i]=in1.v[j1];}
      else                  { out1.v[i]=in2.v[j1];}
      int j2 = i|mask;
      if  ( (i&mask) == 0 ) { out2.v[i]=in1.v[j2];}
      else                  { out2.v[i]=in2.v[j2];}
    }      
  }
  template <typename vec>
  static accelerator_inline void Exchange0(vec &out1,vec &out2,vec &in1,vec &in2){
    ExchangeN<vec,0>(out1,out2,in1,in2);
  };
  template <typename vec>
  static accelerator_inline void Exchange1(vec &out1,vec &out2,vec &in1,vec &in2){
    ExchangeN<vec,1>(out1,out2,in1,in2);
  };
  template <typename vec>
  static accelerator_inline void Exchange2(vec &out1,vec &out2,vec &in1,vec &in2){
    ExchangeN<vec,2>(out1,out2,in1,in2);
  };
  template <typename vec>
  static accelerator_inline void Exchange3(vec &out1,vec &out2,vec &in1,vec &in2){
    ExchangeN<vec,3>(out1,out2,in1,in2);
  };

};

struct Rotate{

  template <int n, typename vec> static accelerator_inline vec tRotate(vec in){
    return rotate(in, n);
  }
    
  template <typename vec>
  static accelerator_inline vec rotate_template(vec in, int n){
    vec out;
    for(int i=0;i<vec::N;i++){
      out.v[i] = in.v[(i + n)%vec::N];
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
    GpuComplexF greduce = in.v[0];
    for(int i=1;i<GpuVectorCF::N;i++) {
      greduce = greduce+in.v[i];
    }
    Grid::ComplexF ret(greduce.z.x,greduce.z.y);
    return ret;
  }

  template<>
  accelerator_inline Grid::ComplexD
  Reduce<Grid::ComplexD, GpuVectorCD>::operator()(GpuVectorCD in)
  {
    GpuComplexD greduce = in.v[0];
    for(int i=1;i<GpuVectorCD::N;i++) {
      greduce = greduce+in.v[i];
    }
    Grid::ComplexD ret(greduce.z.x,greduce.z.y);
    return ret;
  }

  // Real
  template<>
  accelerator_inline Grid::RealF 
  Reduce<RealF, GpuVectorRF>::operator()(GpuVectorRF in)
  {
    RealF ret = in.v[0];
    for(int i=1;i<GpuVectorRF::N;i++) {
      ret = ret+in.v[i];
    }
    return ret;
  }

  template<>
  accelerator_inline Grid::RealD 
  Reduce<RealD, GpuVectorRD>::operator()(GpuVectorRD in)
  {
    RealD ret = in.v[0];
    for(int i=1;i<GpuVectorRD::N;i++) {
      ret = ret+in.v[i];
    }
    return ret;
  }

  template<>
  accelerator_inline Integer
  Reduce<Integer, GpuVectorI>::operator()(GpuVectorI in)
  {
    Integer ret = in.v[0];
    for(int i=1;i<GpuVectorI::N;i++) {
      ret = ret+in.v[i];
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
