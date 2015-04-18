#ifndef VCOMPLEXD_H
#define  VCOMPLEXD_H
#include "Grid.h"
#include "Grid_vComplexF.h"

namespace Grid {
    class vComplexD {
    public:
        zvec v;
    public:
	typedef zvec     vector_type;
	typedef ComplexD scalar_type;

        vComplexD & operator = ( Zero & z){
            vzero(*this);
            return (*this);
        }
        vComplexD(){};
        vComplexD(ComplexD a){
	  vsplat(*this,a);
	};
        vComplexD(double a){
	  vsplat(*this,ComplexD(a));
	};
 
        ///////////////////////////////////////////////
        // mac, mult, sub, add, adj
        // Should do an AVX2 version with mac.
       ///////////////////////////////////////////////
        friend inline void mac (vComplexD * __restrict__ y,const vComplexD * __restrict__ a,const vComplexD *__restrict__ x) {*y = (*a)*(*x)+(*y);};
        friend inline void mult(vComplexD * __restrict__ y,const vComplexD * __restrict__ l,const vComplexD *__restrict__ r) {*y = (*l) * (*r);}
        friend inline void sub (vComplexD * __restrict__ y,const vComplexD * __restrict__ l,const vComplexD *__restrict__ r) {*y = (*l) - (*r);}
        friend inline void add (vComplexD * __restrict__ y,const vComplexD * __restrict__ l,const vComplexD *__restrict__ r) {*y = (*l) + (*r);}
        friend inline vComplexD adj(const vComplexD &in){ return conj(in); }

        //////////////////////////////////
        // Initialise to 1,0,i
        //////////////////////////////////
        friend inline void vone      (vComplexD &ret){ vsplat(ret,1.0,0.0);}
        friend inline void vzero     (vComplexD &ret){ vsplat(ret,0.0,0.0);}
        friend inline void vcomplex_i(vComplexD &ret){ vsplat(ret,0.0,1.0);}
        
        ////////////////////////////////////
        // Arithmetic operator overloads +,-,*
        ////////////////////////////////////
        friend inline vComplexD operator + (vComplexD a, vComplexD b)
        {
            vComplexD ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_add_pd(a.v,b.v);
#endif
#ifdef SSE4
            ret.v = _mm_add_pd(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_add_pd(a.v,b.v);
#endif
#ifdef QPX
            ret.v = vec_add(a.v,b.v);
#endif
            return ret;
        };
        
        friend inline vComplexD operator - (vComplexD a, vComplexD b)
        {
            vComplexD ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_sub_pd(a.v,b.v);
#endif
#ifdef SSE4
            ret.v = _mm_sub_pd(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_sub_pd(a.v,b.v);
#endif
#ifdef QPX
            ret.v = vec_sub(a.v,b.v);
#endif
            return ret;
        };
        
        friend inline vComplexD operator * (vComplexD a, vComplexD b)
        {
            vComplexD ret;
            
            //Multiplicationof (ak+ibk)*(ck+idk)
            // a + i b can be stored as a data structure
            //From intel optimisation reference guide
            /*
             movsldup xmm0, Src1; load real parts into the destination,
             ; a1, a1, a0, a0
             movaps xmm1, src2; load the 2nd pair of complex values, ; i.e. d1, c1, d0, c0
             mulps xmm0, xmm1; temporary results, a1d1, a1c1, a0d0, ; a0c0
             shufps xmm1, xmm1, b1; reorder the real and imaginary ; parts, c1, d1, c0, d0
             movshdup xmm2, Src1; load the imaginary parts into the ; destination, b1, b1, b0, b0
             mulps xmm2, xmm1; temporary results, b1c1, b1d1, b0c0, ; b0d0
             addsubps xmm0, xmm2; b1c1+a1d1, a1c1 -b1d1, b0c0+a0d
             VSHUFPD (VEX.256 encoded version)
             IF IMM0[0] = 0
             THEN DEST[63:0]=SRC1[63:0] ELSE DEST[63:0]=SRC1[127:64] FI;
             IF IMM0[1] = 0
             THEN DEST[127:64]=SRC2[63:0] ELSE DEST[127:64]=SRC2[127:64] FI;
             IF IMM0[2] = 0
             THEN DEST[191:128]=SRC1[191:128] ELSE DEST[191:128]=SRC1[255:192] FI;
             IF IMM0[3] = 0
             THEN DEST[255:192]=SRC2[191:128] ELSE DEST[255:192]=SRC2[255:192] FI;
             */
#if defined (AVX1)|| defined (AVX2)
            zvec ymm0,ymm1,ymm2;
            ymm0 = _mm256_shuffle_pd(a.v,a.v,0x0); // ymm0 <- ar ar, ar,ar b'00,00
            ymm0 = _mm256_mul_pd(ymm0,b.v);        // ymm0 <- ar bi, ar br
            ymm1 = _mm256_shuffle_pd(b.v,b.v,0x5); // ymm1 <- br,bi  b'01,01
            ymm2 = _mm256_shuffle_pd(a.v,a.v,0xF); // ymm2 <- ai,ai  b'11,11
            ymm1 = _mm256_mul_pd(ymm1,ymm2);       // ymm1 <- br ai, ai bi
            ret.v= _mm256_addsub_pd(ymm0,ymm1);
#endif
#ifdef SSE4
            zvec ymm0,ymm1,ymm2;
            ymm0 = _mm_shuffle_pd(a.v,a.v,0x0); // ymm0 <- ar ar,
            ymm0 = _mm_mul_pd(ymm0,b.v);        // ymm0 <- ar bi, ar br
            ymm1 = _mm_shuffle_pd(b.v,b.v,0x1); // ymm1 <- br,bi   b01
            ymm2 = _mm_shuffle_pd(a.v,a.v,0x3); // ymm2 <- ai,ai   b11
            ymm1 = _mm_mul_pd(ymm1,ymm2);       // ymm1 <- br ai, ai bi
            ret.v= _mm_addsub_pd(ymm0,ymm1);
#endif
#ifdef AVX512
            /* This is from
             * Automatic SIMD Vectorization of Fast Fourier Transforms for the Larrabee and AVX Instruction Sets 
             * @inproceedings{McFarlin:2011:ASV:1995896.1995938,
             * author = {McFarlin, Daniel S. and Arbatov, Volodymyr and Franchetti, Franz and P\"{u}schel, Markus},
             * title = {Automatic SIMD Vectorization of Fast Fourier Transforms for the Larrabee and AVX Instruction Sets},
             * booktitle = {Proceedings of the International Conference on Supercomputing},
             * series = {ICS '11},
             * year = {2011},
             * isbn = {978-1-4503-0102-2},
             * location = {Tucson, Arizona, USA},
             * pages = {265--274},
             * numpages = {10},
             * url = {http://doi.acm.org/10.1145/1995896.1995938},
             * doi = {10.1145/1995896.1995938},
             * acmid = {1995938},
             * publisher = {ACM},
             * address = {New York, NY, USA},
             * keywords = {autovectorization, fourier transform, program generation, simd, super-optimization},
 *                } 
             */
            zvec vzero,ymm0,ymm1,real,imag;
            vzero =  _mm512_setzero();
            ymm0 =  _mm512_swizzle_pd(a.v, _MM_SWIZ_REG_CDAB); // 
            real =  _mm512_mask_or_epi64(a.v, 0xAAAA,vzero, ymm0);
            imag =  _mm512_mask_sub_pd(a.v, 0x5555,vzero, ymm0);
            ymm1 =  _mm512_mul_pd(real, b.v);
            ymm0 =  _mm512_swizzle_pd(b.v, _MM_SWIZ_REG_CDAB); // OK
            ret.v=  _mm512_fmadd_pd(ymm0,imag,ymm1);
             /* Imag OK */
#endif
#ifdef QPX
            ret.v = vec_mul(a.v,b.v);
#endif
            return ret;
        };

	////////////////////////////////////////////////////////////////////
	// General permute; assumes vector length is same across 
	// all subtypes; may not be a good assumption, but could
	// add the vector width as a template param for BG/Q for example
	////////////////////////////////////////////////////////////////////
	friend inline void permute(vComplexD &y,vComplexD b,int perm)
	{
	  Gpermute<vComplexD>(y,b,perm);
	}
	friend inline void merge(vComplexD &y,std::vector<ComplexD *> &extracted)
	{
	  Gmerge<vComplexD,ComplexD >(y,extracted);
	}
	friend inline void extract(const vComplexD &y,std::vector<ComplexD *> &extracted)
	{
	  Gextract<vComplexD,ComplexD>(y,extracted);
	}
	friend inline void merge(vComplexD &y,std::vector<ComplexD > &extracted)
	{
	  Gmerge<vComplexD,ComplexD >(y,extracted);
	}
	friend inline void extract(const vComplexD &y,std::vector<ComplexD > &extracted)
	{
	  Gextract<vComplexD,ComplexD>(y,extracted);
	}

        ///////////////////////
        // Splat
        ///////////////////////
        friend inline void vsplat(vComplexD &ret,ComplexD c){
            float a= real(c);
            float b= imag(c);
            vsplat(ret,a,b);
        }


        friend inline void vsplat(vComplexD &ret,double rl,double ig){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_pd(ig,rl,ig,rl);
#endif
#ifdef SSE4
            ret.v = _mm_set_pd(ig,rl);
#endif
#ifdef AVX512
            ret.v = _mm512_set_pd(ig,rl,ig,rl,ig,rl,ig,rl);
#endif
#ifdef QPX
            ret.v = {ig,rl,ig,rl};
#endif
        }

        friend inline void vset(vComplexD &ret,ComplexD *a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_pd(a[1].imag(),a[1].real(),a[0].imag(),a[0].real());
#endif
#ifdef SSE4
            ret.v = _mm_set_pd(a[0].imag(),a[0].real());
#endif
#ifdef AVX512
            ret.v = _mm512_set_pd(a[3].imag(),a[3].real(),a[2].imag(),a[2].real(),a[1].imag(),a[1].real(),a[0].imag(),a[0].real());
            // Note v has a0 a1 a2 a3 
#endif
#ifdef QPX
            ret.v = {a[0].real(),a[0].imag(),a[1].real(),a[3].imag()};
#endif
        }

friend inline void vstore(const vComplexD &ret, ComplexD *a){
#if defined (AVX1)|| defined (AVX2)
       _mm256_store_pd((double *)a,ret.v);
#endif
#ifdef SSE4
       _mm_store_pd((double *)a,ret.v);
#endif
#ifdef AVX512
       _mm512_store_pd((double *)a,ret.v);
   //Note v has a3 a2 a1 a0
#endif
#ifdef QPX
	assert(0);
#endif
        }
      friend inline void vprefetch(const vComplexD &v)
        {
            _mm_prefetch((const char*)&v.v,_MM_HINT_T0);
        }

        ////////////////////////
        // Conjugate
        ////////////////////////
        friend inline vComplexD conj(const vComplexD &in){
            vComplexD ret ; vzero(ret);
#if defined (AVX1)|| defined (AVX2)
            // addsubps 0, inv=>0+in.v[3] 0-in.v[2], 0+in.v[1], 0-in.v[0], ...
            __m256d tmp = _mm256_addsub_pd(ret.v,_mm256_shuffle_pd(in.v,in.v,0x5));
             ret.v=_mm256_shuffle_pd(tmp,tmp,0x5);
#endif
#ifdef SSE4
            ret.v = _mm_addsub_pd(ret.v,in.v);
#endif
#ifdef AVX512
             // Xeon does not have fmaddsub or addsub 
             // with mask 0xa (1010), v[0] -v[1] v[2] -v[3] ....
             ret.v = _mm512_mask_sub_pd(in.v, 0xaaaa,ret.v, in.v);
             
#endif
#ifdef QPX
	     assert(0);
#endif
            return ret;
        }
// REDUCE FIXME must be a cleaner implementation
       friend inline ComplexD Reduce(const vComplexD & in)
       { 
#if defined (AVX1) || defined(AVX2)
	 //            return std::complex<double>(_mm256_mask_reduce_add_pd(0x55, in.v),_mm256_mask_reduce_add_pd(0xAA, in.v));
	 __attribute__ ((aligned(32))) double c_[4];
         _mm256_store_pd(c_,in.v);
	 return ComplexD(c_[0]+c_[2],c_[1]+c_[3]);
#endif 
#ifdef AVX512
            return ComplexD(_mm512_mask_reduce_add_pd(0x5555, in.v),_mm512_mask_reduce_add_pd(0xAAAA, in.v));
#endif 
#ifdef QPX
#endif
        }
        
        // Unary negation
        friend inline vComplexD operator -(const vComplexD &r) {
            vComplexD ret;
            vzero(ret);
            ret = ret - r;
            return ret;
        }
        // *=,+=,-= operators
        inline vComplexD &operator *=(const vComplexD &r) {
            *this = (*this)*r;
            return *this;
        }
        inline vComplexD &operator +=(const vComplexD &r) {
            *this = *this+r;
            return *this;
        }
        inline vComplexD &operator -=(const vComplexD &r) {
            *this = *this-r;
            return *this;
        }

    public:
        static int Nsimd(void) { return sizeof(zvec)/sizeof(double)/2;}
    };


    inline vComplexD innerProduct(const vComplexD & l, const vComplexD & r) { return conj(l)*r; }


    typedef  vComplexD vDComplex;
    inline void zeroit(vComplexD &z){ vzero(z);}

    inline vComplexD outerProduct(const vComplexD &l, const vComplexD& r)
    {
        return l*r;
    }
    inline vComplexD trace(const vComplexD &arg){
        return arg;
    }
/////////////////////////////////////////////////////////////////////////
//// Generic routine to promote object<complex> -> object<vcomplex>
//// Supports the array reordering transformation that gives me SIMD utilisation
///////////////////////////////////////////////////////////////////////////
/*
template<template<class> class object>
inline object<vComplex> splat(object<Complex >s){
      object<vComplex> ret;
      vComplex * v_ptr = (vComplex *)& ret;
      Complex * s_ptr = (Complex *) &s;
      for(int i=0;i<sizeof(ret);i+=sizeof(vComplex)){
          vsplat(*(v_ptr++),*(s_ptr++));
      }
      return ret;
    }
*/
}
#endif
