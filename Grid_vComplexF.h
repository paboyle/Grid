#ifndef VCOMPLEXF
#define VCOMPLEXF
#include "Grid.h"

namespace Grid {

  /*
  inline void Print(const char *A,cvec c) { 
	  float *fp=(float *)&c; 
	  printf(A); 
	  printf(" %le %le %le %le %le %le %le %le\n",
		 fp[0],fp[1],fp[2],fp[3],fp[4],fp[5],fp[6],fp[7]);
	}
  */

    class vComplexF {
      //    protected:

    public:
        cvec v;
        
    public:
        static inline int Nsimd(void) { return sizeof(cvec)/sizeof(float)/2;}
    public:
	typedef cvec     vector_type;
	typedef ComplexF scalar_type;

        vComplexF & operator = ( Zero & z){
            vzero(*this);
            return (*this);
        }
        vComplexF(){};
        vComplexF(ComplexF a){
	  vsplat(*this,a);
	};
        vComplexF(double a){
	  vsplat(*this,ComplexF(a));
	};
       
        ///////////////////////////////////////////////
        // mac, mult, sub, add, adj
        // Should do an AVX2 version with mac.
        ///////////////////////////////////////////////
        friend inline void mac (vComplexF * __restrict__ y,const vComplexF * __restrict__ a,const vComplexF *__restrict__ x){ *y = (*a)*(*x)+(*y); };
        friend inline void mult(vComplexF * __restrict__ y,const vComplexF * __restrict__ l,const vComplexF *__restrict__ r){ *y = (*l) * (*r); }
        friend inline void sub (vComplexF * __restrict__ y,const vComplexF * __restrict__ l,const vComplexF *__restrict__ r){ *y = (*l) - (*r); }
        friend inline void add (vComplexF * __restrict__ y,const vComplexF * __restrict__ l,const vComplexF *__restrict__ r){ *y = (*l) + (*r); }
        friend inline vComplexF adj(const vComplexF &in){ return conj(in); }
        
        //////////////////////////////////
        // Initialise to 1,0,i
        //////////////////////////////////
        friend inline void vone(vComplexF &ret)      { vsplat(ret,1.0,0.0); }
        friend inline void vzero(vComplexF &ret)     { vsplat(ret,0.0,0.0); }
        friend inline void vcomplex_i(vComplexF &ret){ vsplat(ret,0.0,1.0);}
          
        ////////////////////////////////////
        // Arithmetic operator overloads +,-,*
        ////////////////////////////////////
        friend inline vComplexF operator + (vComplexF a, vComplexF b)
        {
            vComplexF ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_add_ps(a.v,b.v);
#endif
#ifdef SSE2
            ret.v = _mm_add_ps(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_add_ps(a.v,b.v);
#endif
#ifdef QPX
#error
#endif
            return ret;
        };
        
        friend inline vComplexF operator - (vComplexF a, vComplexF b)
        {
            vComplexF ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_sub_ps(a.v,b.v);
#endif
#ifdef SSE2
            ret.v = _mm_sub_ps(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_sub_ps(a.v,b.v);
#endif
#ifdef QPX
#error
#endif
            return ret;
        };
        
        friend inline vComplexF operator * (vComplexF a, vComplexF b)
        {
            vComplexF ret;
            
            //Multiplicationof (ak+ibk)*(ck+idk)
            // a + i b can be stored as a data structure
            //From intel optimisation reference
            /*
             movsldup xmm0, Src1; load real parts into the destination,
             ; a1, a1, a0, a0
             movaps xmm1, src2; load the 2nd pair of complex values, ; i.e. d1, c1, d0, c0
             mulps xmm0, xmm1; temporary results, a1d1, a1c1, a0d0, ; a0c0
             shufps xmm1, xmm1, b1; reorder the real and imaginary ; parts, c1, d1, c0, d0
             movshdup xmm2, Src1; load the imaginary parts into the ; destination, b1, b1, b0, b0
             mulps xmm2, xmm1; temporary results, b1c1, b1d1, b0c0, ; b0d0
             addsubps xmm0, xmm2; b1c1+a1d1, a1c1 -b1d1, b0c0+a0d
             */
#if defined (AVX1)|| defined (AVX2)
            cvec ymm0,ymm1,ymm2;
            ymm0 = _mm256_shuffle_ps(a.v,a.v,_MM_SHUFFLE(2,2,0,0)); // ymm0 <- ar ar,
            ymm0 = _mm256_mul_ps(ymm0,b.v);        // ymm0 <- ar bi, ar br
            // FIXME AVX2 could MAC
            ymm1 = _mm256_shuffle_ps(b.v,b.v,_MM_SHUFFLE(2,3,0,1)); // ymm1 <- br,bi
            ymm2 = _mm256_shuffle_ps(a.v,a.v,_MM_SHUFFLE(3,3,1,1)); // ymm2 <- ai,ai
            ymm1 = _mm256_mul_ps(ymm1,ymm2);       // ymm1 <- br ai, ai bi
            ret.v= _mm256_addsub_ps(ymm0,ymm1);    
#endif
#ifdef SSE2
            cvec ymm0,ymm1,ymm2;
            ymm0 = _mm_shuffle_ps(a.v,a.v,_MM_SHUFFLE(2,2,0,0)); // ymm0 <- ar ar,
            ymm0 = _mm_mul_ps(ymm0,b.v);        // ymm0 <- ar bi, ar br
            ymm1 = _mm_shuffle_ps(b.v,b.v,_MM_SHUFFLE(2,3,0,1)); // ymm1 <- br,bi
            ymm2 = _mm_shuffle_ps(a.v,a.v,_MM_SHUFFLE(3,3,1,1)); // ymm2 <- ai,ai
            ymm1 = _mm_mul_ps(ymm1,ymm2);       // ymm1 <- br ai, ai bi
            ret.v= _mm_addsub_ps(ymm0,ymm1);
#endif
#ifdef AVX512
//
            cvec vzero,ymm0,ymm1,real, imag;
            vzero = _mm512_setzero();
            ymm0  = _mm512_swizzle_ps(a.v, _MM_SWIZ_REG_CDAB); // 
            real  = _mm512_mask_or_epi32(a.v, 0xAAAA,vzero, ymm0);
            imag  = _mm512_mask_sub_ps(a.v, 0x5555,vzero, ymm0);
            ymm1  = _mm512_mul_ps(real, b.v);
            ymm0  = _mm512_swizzle_ps(b.v, _MM_SWIZ_REG_CDAB); // OK
            ret.v = _mm512_fmadd_ps(ymm0,imag,ymm1);


#endif
#ifdef QPX
            ret.v = vec_mul(a.v,b.v);
#endif
            return ret;
        };
      

	////////////////////////////////////////////////////////////////////////
	// FIXME:  gonna remove these load/store, get, set, prefetch
	////////////////////////////////////////////////////////////////////////
        friend inline void vset(vComplexF &ret, ComplexF *a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_ps(a[3].imag(),a[3].real(),a[2].imag(),a[2].real(),a[1].imag(),a[1].real(),a[0].imag(),a[0].real());
#endif
#ifdef SSE2
            ret.v = _mm_set_ps(a[1].imag, a[1].real(),a[0].imag(),a[0].real());
#endif
#ifdef AVX512
            ret.v = _mm512_set_ps(a[7].imag(),a[7].real(),a[6].imag(),a[6].real(),a[5].imag(),a[5].real(),a[4].imag(),a[4].real(),a[3].imag(),a[3].real(),a[2].imag(),a[2].real(),a[1].imag(),a[1].real(),a[0].imag(),a[0].real());
            // Note v has a0 a1 a2 a3 a4 a5 a6 a7
#endif
#ifdef QPX
	    ret.v = {a[0].real(),a[0].imag(),a[1].real(),a[1].imag(),a[2].real(),a[2].imag(),a[3].real(),a[3].imag()};
#endif
        }
        
        ///////////////////////
        // Splat
        ///////////////////////
        friend inline void vsplat(vComplexF &ret,ComplexF c){
            float a= real(c);
            float b= imag(c);
            vsplat(ret,a,b);
        }

friend inline void vstore(const vComplexF &ret, ComplexF *a){
#if defined (AVX1)|| defined (AVX2)
        _mm256_store_ps((float *)a,ret.v);
#endif
#ifdef SSE2
        _mm_store_ps((float *)a,ret.v);
#endif
#ifdef AVX512
	_mm512_store_ps((float *)a,ret.v);
//Note v has a3 a2 a1 a0
#endif
#ifdef QPX
	assert(0);
#endif
}
      friend inline void vprefetch(const vComplexF &v)
        {
            _mm_prefetch((const char*)&v.v,_MM_HINT_T0);
        }

        friend inline void vsplat(vComplexF &ret,float a,float b){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_ps(b,a,b,a,b,a,b,a);
#endif
#ifdef SSE2
            ret.v = _mm_set_ps(a,b,a,b);
#endif
#ifdef AVX512
            ret.v = _mm512_set_ps(b,a,b,a,b,a,b,a,b,a,b,a,b,a,b,a);
#endif
#ifdef QPX
            ret.v = {a,b,a,b};
#endif
        }
       friend inline ComplexF Reduce(const vComplexF & in)
       {
#if defined (AVX1) || defined(AVX2)
         __attribute__ ((aligned(32))) float c_[8];
         _mm256_store_ps(c_,in.v);
         return ComplexF(c_[0]+c_[2]+c_[4]+c_[6],c_[1]+c_[3]+c_[5]+c_[7]);

#endif
#ifdef AVX512
            return ComplexF(_mm512_mask_reduce_add_ps(0x5555, in.v),_mm512_mask_reduce_add_ps(0xAAAA, in.v));
#endif
#ifdef QPX
#endif
        }

        friend inline vComplexF operator * (const ComplexF &a, vComplexF b){
            vComplexF va;
            vsplat(va,a);
            return va*b;
        }
        friend inline vComplexF operator * (vComplexF b,const ComplexF &a){
	  return a*b;
        }

       /*
	template<class real>
        friend inline vComplexF operator * (vComplexF b,const real &a){
            vComplexF va;
	    Complex ca(a,0);
            vsplat(va,ca);
            return va*b;
        }
	template<class real>
	friend inline vComplexF operator * (const real &a,vComplexF b){
	  return a*b;
	}

        friend inline vComplexF operator + (const Complex &a, vComplexF b){
            vComplexF va;
            vsplat(va,a);
            return va+b;
        }
        friend inline vComplexF operator + (vComplexF b,const Complex &a){
            return a+b;
        }
	template<class real>
        friend inline vComplexF operator + (vComplexF b,const real &a){
            vComplexF va;
	    Complex ca(a,0);
            vsplat(va,ca);
            return va+b;
        }
	template<class real>
	friend inline vComplexF operator + (const real &a,vComplexF b){
	  return a+b;
	}
        friend inline vComplexF operator - (const Complex &a, vComplexF b){
            vComplexF va;
            vsplat(va,a);
            return va-b;
        }
        friend inline vComplexF operator - (vComplexF b,const Complex &a){
            vComplexF va;
            vsplat(va,a);
            return b-va;
        }
	template<class real>
        friend inline vComplexF operator - (vComplexF b,const real &a){
            vComplexF va;
	    Complex ca(a,0);
            vsplat(va,ca);
            return b-va;
        }
	template<class real>
	friend inline vComplexF operator - (const real &a,vComplexF b){
            vComplexF va;
	    Complex ca(a,0);
            vsplat(va,ca);
            return va-b;
	}
       */
       

        ///////////////////////
        // Conjugate
        ///////////////////////
								     
        friend inline vComplexF conj(const vComplexF &in){
            vComplexF ret ; vzero(ret);
#if defined (AVX1)|| defined (AVX2)
             cvec tmp;
             tmp = _mm256_addsub_ps(ret.v,_mm256_shuffle_ps(in.v,in.v,_MM_SHUFFLE(2,3,0,1))); // ymm1 <- br,bi
             ret.v=_mm256_shuffle_ps(tmp,tmp,_MM_SHUFFLE(2,3,0,1));
#endif
#ifdef SSE2
            ret.v = _mm_addsub_ps(ret.v,in.v);
#endif
#ifdef AVX512
            ret.v = _mm512_mask_sub_ps(in.v,0xaaaa,ret.v,in.v); // Zero out 0+real 0-imag 
#endif
#ifdef QPX
            assert(0);
#endif
            return ret;
        }
        
        // Unary negation
        friend inline vComplexF operator -(const vComplexF &r) {
            vComplexF ret;
            vzero(ret);
            ret = ret - r;
            return ret;
        }
        // *=,+=,-= operators
        inline vComplexF &operator *=(const vComplexF &r) {
            *this = (*this)*r;
            return *this;
        }
        inline vComplexF &operator +=(const vComplexF &r) {
            *this = *this+r;
            return *this;
        }
        inline vComplexF &operator -=(const vComplexF &r) {
            *this = *this-r;
            return *this;
        }

      friend inline void permute(vComplexF &y,vComplexF b,int perm)
      {
	Gpermute<vComplexF>(y,b,perm);
      }
      friend inline void merge(vComplexF &y,std::vector<ComplexF *> &extracted)
      {
	Gmerge<vComplexF,ComplexF >(y,extracted);
      }
      friend inline void extract(const vComplexF &y,std::vector<ComplexF *> &extracted)
      {
	Gextract<vComplexF,ComplexF>(y,extracted);
      }
      friend inline void merge(vComplexF &y,std::vector<ComplexF > &extracted)
      {
	Gmerge<vComplexF,ComplexF >(y,extracted);
      }
      friend inline void extract(const vComplexF &y,std::vector<ComplexF > &extracted)
      {
	Gextract<vComplexF,ComplexF>(y,extracted);
      }


    };

    inline vComplexF innerProduct(const vComplexF & l, const vComplexF & r) 
    {
      return conj(l)*r; 
    }

    inline void zeroit(vComplexF &z){ vzero(z);}

    inline vComplexF outerProduct(const vComplexF &l, const vComplexF& r)
    {
        return l*r;
    }
    inline vComplexF trace(const vComplexF &arg){
        return arg;
    }
}
#endif
