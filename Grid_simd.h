#ifndef GRID_SIMD_H
#define GRID_SIMD_H

////////////////////////////////////////////////////////////////////////
// Define scalar and vector floating point types
//
// Scalar:   RealF, RealD, ComplexF, ComplexD
//
// Vector:  vRealF, vRealD, vComplexF, vComplexD
//
// Vector types are arch dependent
////////////////////////////////////////////////////////////////////////
    

#ifdef SSE4
#include <pmmintrin.h>
#endif
#if defined(AVX1) || defined (AVX2)
#include <immintrin.h>
#endif
#ifdef AVX512
#include <immintrin.h>
#endif

namespace Grid {

  typedef  float  RealF;
  typedef  double RealD;
  
  typedef std::complex<RealF> ComplexF;
  typedef std::complex<RealD> ComplexD;


  inline RealF adj(const RealF  & r){ return r; }
  inline RealF conj(const RealF  & r){ return r; }
  inline RealF real(const RealF  & r){ return r; }

  inline RealD adj(const RealD  & r){ return r; }
  inline RealD conj(const RealD  & r){ return r; }
  inline RealD real(const RealD  & r){ return r; }

  inline ComplexD innerProduct(const ComplexD & l, const ComplexD & r) { return conj(l)*r; }
  inline ComplexF innerProduct(const ComplexF & l, const ComplexF & r) { return conj(l)*r; }
  inline RealD innerProduct(const RealD & l, const RealD & r) { return l*r; }
  inline RealF innerProduct(const RealF & l, const RealF & r) { return l*r; }

    ////////////////////////////////////////////////////////////////////////////////
    //Provide support functions for basic real and complex data types required by Grid
    //Single and double precision versions. Should be able to template this once only.
    ////////////////////////////////////////////////////////////////////////////////
    inline void mac (ComplexD * __restrict__ y,const ComplexD * __restrict__ a,const ComplexD *__restrict__ x){ *y = (*a) * (*x)+(*y); };
    inline void mult(ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) * (*r);}
    inline void sub (ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) - (*r);}
    inline void add (ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) + (*r);}
    inline ComplexD adj(const ComplexD& r){ return(conj(r)); }
    // conj already supported for complex
    
    inline void mac (ComplexF * __restrict__ y,const ComplexF * __restrict__ a,const ComplexF *__restrict__ x){ *y = (*a) * (*x)+(*y); }
    inline void mult(ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) * (*r); }
    inline void sub (ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) - (*r); }
    inline void add (ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) + (*r); }
    inline ComplexF  adj(const ComplexF& r ){ return(conj(r)); }
    //conj already supported for complex
    
    inline void mac (RealD * __restrict__ y,const RealD * __restrict__ a,const RealD *__restrict__ x){  *y = (*a) * (*x)+(*y);}
    inline void mult(RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) * (*r);}
    inline void sub (RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) - (*r);}
    inline void add (RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) + (*r);}
    
    inline void mac (RealF * __restrict__ y,const RealF * __restrict__ a,const RealF *__restrict__ x){  *y = (*a) * (*x)+(*y); }
    inline void mult(RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) * (*r); }
    inline void sub (RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) - (*r); }
    inline void add (RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) + (*r); }
    


  class Zero{};
  static Zero zero;
  template<class itype> inline void zeroit(itype &arg){ arg=zero;};
  template<>            inline void zeroit(ComplexF &arg){ arg=0; };
  template<>            inline void zeroit(ComplexD &arg){ arg=0; };
  template<>            inline void zeroit(RealF &arg){ arg=0; };
  template<>            inline void zeroit(RealD &arg){ arg=0; };


#if defined (SSE4)
    typedef __m128 fvec;
    typedef __m128d dvec;
    typedef __m128 cvec;
    typedef __m128d zvec;
    typedef __m128i ivec;
#endif
#if defined (AVX1) || defined (AVX2)
    typedef __m256 fvec;
    typedef __m256d dvec;
    typedef __m256  cvec;
    typedef __m256d zvec;
    typedef __m256i ivec;
#endif
#if defined (AVX512)
    typedef __m512  fvec;
    typedef __m512d dvec;
    typedef __m512  cvec;
    typedef __m512d zvec;
    typedef __m512i ivec;
#endif
#if defined (QPX)
    typedef float  fvec __attribute__ ((vector_size (16))); // QPX has same SIMD width irrespective of precision
    typedef float  cvec __attribute__ ((vector_size (16)));
    
    typedef vector4double dvec;
    typedef vector4double zvec;
#endif
#if defined (AVX1) || defined (AVX2) || defined (AVX512)
    inline void v_prefetch0(int size, const char *ptr){
          for(int i=0;i<size;i+=64){ //  Define L1 linesize above// What about SSE?
            _mm_prefetch(ptr+i+4096,_MM_HINT_T1);
            _mm_prefetch(ptr+i+512,_MM_HINT_T0);
          }
    }
#else 
    inline void v_prefetch0(int size, const char *ptr){};
#endif


/////////////////////////////////////////////////////////////////
// Generic extract/merge/permute
/////////////////////////////////////////////////////////////////
template<class vsimd,class scalar>
inline void Gextract(const vsimd &y,std::vector<scalar *> &extracted){
  // FIXME: bounce off stack is painful
  // temporary hack while I figure out better way.
  // There are intrinsics to do this work without the storage.
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  std::vector<scalar,alignedAllocator<scalar> > buf(Nsimd); 
  vstore(y,&buf[0]);
  for(int i=0;i<Nextr;i++){
    *extracted[i] = buf[i*s];
    extracted[i]++;
  }
};
template<class vsimd,class scalar>
inline void Gmerge(vsimd &y,std::vector<scalar *> &extracted){
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  std::vector<scalar> buf(Nsimd); 
  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      buf[i*s+ii]=*extracted[i];
    }
    extracted[i]++;
  }
  vset(y,&buf[0]); 
};
template<class vsimd,class scalar>
inline void Gextract(const vsimd &y,std::vector<scalar> &extracted){
  // FIXME: bounce off stack is painful
  // temporary hack while I figure out better way.
  // There are intrinsics to do this work without the storage.
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  std::vector<scalar,alignedAllocator<scalar> > buf(Nsimd); 

  vstore(y,&buf[0]);

  for(int i=0;i<Nextr;i++){
    extracted[i] = buf[i*s];
  }
};
template<class vsimd,class scalar>
inline void Gmerge(vsimd &y,std::vector<scalar> &extracted){
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  std::vector<scalar> buf(Nsimd); 
  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      buf[i*s+ii]=extracted[i];
    }
  }
  vset(y,&buf[0]); 
};

//////////////////////////////////////////////////////////
// Permute
// Permute 0 every ABCDEFGH -> BA DC FE HG
// Permute 1 every ABCDEFGH -> CD AB GH EF
// Permute 2 every ABCDEFGH -> EFGH ABCD
// Permute 3 possible on longer iVector lengths (512bit = 8 double = 16 single)
// Permute 4 possible on half precision @512bit vectors.
//////////////////////////////////////////////////////////
template<class vsimd>
inline void Gpermute(vsimd &y,const vsimd &b,int perm){
      switch (perm){
#if defined(AVX1)||defined(AVX2)
      // 8x32 bits=>3 permutes
      case 2: y.v = _mm256_shuffle_ps(b.v,b.v,_MM_SHUFFLE(2,3,0,1)); break;
      case 1: y.v = _mm256_shuffle_ps(b.v,b.v,_MM_SHUFFLE(1,0,3,2)); break;
      case 0: y.v = _mm256_permute2f128_ps(b.v,b.v,0x01); break;
#endif
#ifdef SSE4
      case 1: y.v = _mm_shuffle_ps(b.v,b.v,_MM_SHUFFLE(2,3,0,1)); break;
      case 0: y.v = _mm_shuffle_ps(b.v,b.v,_MM_SHUFFLE(1,0,3,2));break;
#endif
#ifdef AVX512
	// 16 floats=> permutes
        // Permute 0 every abcd efgh ijkl mnop -> badc fehg jilk nmpo 
        // Permute 1 every abcd efgh ijkl mnop -> cdab ghef jkij opmn 
        // Permute 2 every abcd efgh ijkl mnop -> efgh abcd mnop ijkl
        // Permute 3 every abcd efgh ijkl mnop -> ijkl mnop abcd efgh
      case 3: y.v = _mm512_swizzle_ps(b.v,_MM_SWIZ_REG_CDAB); break;
      case 2: y.v = _mm512_swizzle_ps(b.v,_MM_SWIZ_REG_BADC); break;
      case 1: y.v = _mm512_permute4f128_ps(b.v,(_MM_PERM_ENUM)_MM_SHUFFLE(2,3,0,1)); break;
      case 0: y.v = _mm512_permute4f128_ps(b.v,(_MM_PERM_ENUM)_MM_SHUFFLE(1,0,3,2)); break;
#endif
#ifdef QPX
#error not implemented
#endif
      default: assert(0); break;
      }
    };
};

#include <Grid_vInteger.h>
#include <Grid_vRealF.h>
#include <Grid_vRealD.h>
#include <Grid_vComplexF.h>
#include <Grid_vComplexD.h>

namespace Grid {

  // NB: Template the following on "type Complex" and then implement *,+,- for 
  // ComplexF, ComplexD, RealF, RealD above to
  // get full generality of binops with scalars.
   inline void mac (vComplexF *__restrict__ y,const ComplexF *__restrict__ a,const vComplexF *__restrict__ x){ *y = (*a)*(*x)+(*y); };
   inline void mult(vComplexF *__restrict__ y,const ComplexF *__restrict__ l,const vComplexF *__restrict__ r){ *y = (*l) * (*r); }
   inline void sub (vComplexF *__restrict__ y,const ComplexF *__restrict__ l,const vComplexF *__restrict__ r){ *y = (*l) - (*r); }
   inline void add (vComplexF *__restrict__ y,const ComplexF *__restrict__ l,const vComplexF *__restrict__ r){ *y = (*l) + (*r); }
   inline void mac (vComplexF *__restrict__ y,const vComplexF *__restrict__ a,const ComplexF *__restrict__ x){ *y = (*a)*(*x)+(*y); };
   inline void mult(vComplexF *__restrict__ y,const vComplexF *__restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) * (*r); }
   inline void sub (vComplexF *__restrict__ y,const vComplexF *__restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) - (*r); }
   inline void add (vComplexF *__restrict__ y,const vComplexF *__restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) + (*r); }

   inline void mac (vComplexD *__restrict__ y,const ComplexD *__restrict__ a,const vComplexD *__restrict__ x){ *y = (*a)*(*x)+(*y); };
   inline void mult(vComplexD *__restrict__ y,const ComplexD *__restrict__ l,const vComplexD *__restrict__ r){ *y = (*l) * (*r); }
   inline void sub (vComplexD *__restrict__ y,const ComplexD *__restrict__ l,const vComplexD *__restrict__ r){ *y = (*l) - (*r); }
   inline void add (vComplexD *__restrict__ y,const ComplexD *__restrict__ l,const vComplexD *__restrict__ r){ *y = (*l) + (*r); }
   inline void mac (vComplexD *__restrict__ y,const vComplexD *__restrict__ a,const ComplexD *__restrict__ x){ *y = (*a)*(*x)+(*y); };
   inline void mult(vComplexD *__restrict__ y,const vComplexD *__restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) * (*r); }
   inline void sub (vComplexD *__restrict__ y,const vComplexD *__restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) - (*r); }
   inline void add (vComplexD *__restrict__ y,const vComplexD *__restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) + (*r); }

   inline void mac (vRealF *__restrict__ y,const RealF *__restrict__ a,const vRealF *__restrict__ x){ *y = (*a)*(*x)+(*y); };
   inline void mult(vRealF *__restrict__ y,const RealF *__restrict__ l,const vRealF *__restrict__ r){ *y = (*l) * (*r); }
   inline void sub (vRealF *__restrict__ y,const RealF *__restrict__ l,const vRealF *__restrict__ r){ *y = (*l) - (*r); }
   inline void add (vRealF *__restrict__ y,const RealF *__restrict__ l,const vRealF *__restrict__ r){ *y = (*l) + (*r); }
   inline void mac (vRealF *__restrict__ y,const vRealF *__restrict__ a,const RealF *__restrict__ x){ *y = (*a)*(*x)+(*y); };
   inline void mult(vRealF *__restrict__ y,const vRealF *__restrict__ l,const RealF *__restrict__ r){ *y = (*l) * (*r); }
   inline void sub (vRealF *__restrict__ y,const vRealF *__restrict__ l,const RealF *__restrict__ r){ *y = (*l) - (*r); }
   inline void add (vRealF *__restrict__ y,const vRealF *__restrict__ l,const RealF *__restrict__ r){ *y = (*l) + (*r); }

   inline void mac (vRealD *__restrict__ y,const RealD *__restrict__ a,const vRealD *__restrict__ x){ *y = (*a)*(*x)+(*y); };
   inline void mult(vRealD *__restrict__ y,const RealD *__restrict__ l,const vRealD *__restrict__ r){ *y = (*l) * (*r); }
   inline void sub (vRealD *__restrict__ y,const RealD *__restrict__ l,const vRealD *__restrict__ r){ *y = (*l) - (*r); }
   inline void add (vRealD *__restrict__ y,const RealD *__restrict__ l,const vRealD *__restrict__ r){ *y = (*l) + (*r); }
   inline void mac (vRealD *__restrict__ y,const vRealD *__restrict__ a,const RealD *__restrict__ x){ *y = (*a)*(*x)+(*y); };
   inline void mult(vRealD *__restrict__ y,const vRealD *__restrict__ l,const RealD *__restrict__ r){ *y = (*l) * (*r); }
   inline void sub (vRealD *__restrict__ y,const vRealD *__restrict__ l,const RealD *__restrict__ r){ *y = (*l) - (*r); }
   inline void add (vRealD *__restrict__ y,const vRealD *__restrict__ l,const RealD *__restrict__ r){ *y = (*l) + (*r); }

  // Default precision
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
  typedef RealD   Real;
  typedef vRealD vReal;
  typedef vComplexD vComplex;
  typedef std::complex<Real>  Complex;
#else
  typedef RealF  Real;
  typedef vRealF vReal;
  typedef vComplexF vComplex;
  typedef std::complex<Real>  Complex;
#endif
}
#endif
