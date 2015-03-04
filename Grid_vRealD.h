#ifndef VREALD_H
#define VREALD_H

#include "Grid.h"

namespace dpo{
    class vRealD  {
    protected:
        dvec v; // dvec is double precision vector
    public:
        vRealD(){};

        friend inline void mult(vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) * (*r);}
        friend inline void sub (vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) - (*r);}
        friend inline void add (vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) + (*r);}
        friend inline vRealD adj(const vRealD &in) { return in; }
        friend inline vRealD conj(const vRealD &in){ return in; }
        
        friend inline void mac (vRealD &y,const vRealD a,const vRealD x){
#if defined (AVX1) || defined (SSE2)
            y = a*x+y;
#endif
#ifdef AVX2     // AVX 2 introduced FMA support. FMA4 eliminates a copy, but AVX only has FMA3
            // accelerates multiply accumulate, but not general multiply add
            y.v = _mm256_fmadd_pd(a.v,x.v,y.v);
#endif
#ifdef AVX512
            // here precision of vector are still single
            y.v = _mm512_fmadd_pd(a.v,x.v,y.v);
#endif
#ifdef QPX
            y.v = vec_madd(a.v,x.v,y.v);
#endif
        }
        //////////////////////////////////
        // Initialise to 1,0
        //////////////////////////////////
        friend inline void vone (vRealD &ret){ vsplat(ret,1.0);}
        friend inline void vzero(vRealD &ret){ vsplat(ret,0.0);}
        
        
        ////////////////////////////////////
        // Arithmetic operator overloads +,-,*
        ////////////////////////////////////
        friend inline vRealD operator + (vRealD a, vRealD b)
        {
            vRealD ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_add_pd(a.v,b.v);
#endif
#ifdef SSE2
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
        friend inline vRealD operator - (vRealD a, vRealD b)
        {
            vRealD ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_sub_pd(a.v,b.v);
#endif
#ifdef SSE2
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
        
        friend inline vRealD operator * (vRealD a, vRealD b)
        {
            vRealD ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_mul_pd(a.v,b.v);
#endif
#ifdef SSE2
            ret.v = _mm_mul_pd(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_mul_pd(a.v,b.v);
#endif
#ifdef QPX
            ret.v = vec_mul(a.v,b.v);
#endif
            return ret;
        };
        
        // Permute plans
        // Permute 0 every ABCDEFGH -> BA DC FE HG
        // Permute 1 every ABCDEFGH -> CD AB GH EF
        // Permute 2 every ABCDEFGH -> EFGH ABCD
        // Permute 3 possible on longer iVector lengths (512bit = 8 double = 16 single)
        // Permute 4 possible on half precision @512bit vectors.
        friend inline void permute(vRealD &y,vRealD b,int perm){
            switch (perm){
                    // 4 doubles=>2 permutes
#if defined(AVX1)||defined(AVX2)
                case 0: y.v = _mm256_shuffle_pd(b.v,b.v,0x5); break;
                case 1: y.v = _mm256_permute2f128_pd(b.v,b.v,0x01); break;
#endif
#ifdef SSE2
                case 0: y.v = _mm_shuffle_pd(b.v,b.v,0x1); break;
#endif
#ifdef AVX512
                    // 8 double => 3 permutes
        // Permute 0 every abcd efgh -> badc fehg 
        // Permute 1 every abcd efgh -> cdab ghef 
        // Permute 2 every abcd efgh -> efgh abcd 
        // NOTE: mm_512_permutex_pd not implemented
        // NOTE: ignore warning
                case 0: y.v = _mm512_swizzle_pd(b.v,_MM_SWIZ_REG_CDAB); break;
                case 1: y.v = _mm512_swizzle_pd(b.v,_MM_SWIZ_REG_BADC); break;
                case 2: y.v = _mm512_permute4f128_ps(b.v,(_MM_PERM_ENUM)_MM_SHUFFLE(1,0,3,2)); break;
      
#endif
#ifdef QPX
#error
#endif
                default: exit(EXIT_FAILURE); break;
            }
        };
// gona be bye bye
        void vload(dvec& a){
          this->v = a;
        }
        dvec vget(){
          return this->v ;
        }
        
        friend inline void vsplat(vRealD &ret,double a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_pd(a,a,a,a);
#endif
#ifdef SSE2
            ret.v = _mm_set_pd(a,a);
#endif
#ifdef AVX512
            ret.v = _mm512_set1_pd(a);
#endif
#ifdef QPX
            ret.v = {a,a,a,a};
#endif
        }
	friend inline void vset(vRealD &ret, double *a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_pd(a[3],a[2],a[1],a[0]);
#endif
#ifdef SSE2
            ret.v = _mm_set_pd(a[0],a[1]);
#endif
#ifdef AVX512
            ret.v = _mm512_set_pd(a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
            // Note v has a0 a1 a2 a3 a4 a5 a6 a7
#endif
#ifdef QPX
            ret.v = {a[0],a[1],a[2],a[3]};
#endif
	}

	friend inline void vstore(vRealD &ret, double *a){
#if defined (AVX1)|| defined (AVX2)
            _mm256_store_pd(a,ret.v);
#endif
#ifdef SSE2
            _mm_store_pd(a,ret.v);
#endif
#ifdef AVX512
            _mm512_store_pd(a,ret.v);
            // Note v has a7 a6 a5ba4 a3 a2 a1 a0
#endif
#ifdef QPX
          printf("%s Not implemented\n",__func__); 
          exit(-1);
#endif
	}
        friend inline void vprefetch(const vRealD &v)
        {
            _mm_prefetch((const char*)&v.v,_MM_HINT_T0);
        }
        // Unary negation
        friend inline vRealD operator -(const vRealD &r) {
            vRealD ret;
            vzero(ret);
            ret = ret - r;
            return ret;
        }

       friend inline RealD Reduce(const vRealD & in)
       {
#if defined (AVX1) || defined(AVX2)
	 typedef union  {
	   uint64_t l;
	   double   d;
	 } my_conv_t;
	 my_conv_t converter;
// more reduce_add
/*
            __attribute__ ((aligned(32))) double c_[16];
	    __m256d tmp  = _mm256_permute2f128_pd(in.v,in.v,0x01); // tmp 1032; in= 3210
            __m256d hadd = _mm256_hadd_pd(in.v,tmp);              // hadd = 1+0,3+2,3+2,1+0
  	             tmp = _mm256_permute2f128_pd(hadd,hadd,0x01);// tmp  = 3+2,1+0,1+0,3+2
                    hadd = _mm256_hadd_pd(tmp,tmp);               // tmp  = 3+2+1+0,3+2+1+0,1+0+3+2,1+0+3+2
                    _mm256_store_pd(c_,hadd);ô
             return c[0]
*/
	    __m256d tmp  = _mm256_permute2f128_pd(in.v,in.v,0x01); // tmp 1032; in= 3210
            __m256d hadd = _mm256_hadd_pd(in.v,tmp);              // hadd = 1+0,3+2,3+2,1+0
                    hadd = _mm256_hadd_pd(hadd,hadd);             // hadd = 1+0+3+2...
	    converter.l = _mm256_extract_epi64(hadd,0);
            return converter.d;
#endif
#ifdef AVX512
            return _mm512_reduce_add_pd(in.v);
/*
            __attribute__ ((aligned(32))) double c_[8];
           _mm512_store_pd(c_,in.v);
            return c_[0]+c_[1]+c_[2]+c_[3]+c_[4]+c_[5]+c_[6]+c_[7];
*/
#endif
#ifdef QPX
#endif
        }

        // *=,+=,-= operators
        inline vRealD &operator *=(const vRealD &r) {
            *this = (*this)*r;
            return *this;
        }
        inline vRealD &operator +=(const vRealD &r) {
            *this = *this+r;
            return *this;
        }
        inline vRealD &operator -=(const vRealD &r) {
            *this = *this-r;
            return *this;
        }

    public:
        static int Nsimd(void) { return sizeof(dvec)/sizeof(double);}
    };

   inline vRealD localInnerProduct(const vRealD & l, const vRealD & r) { return conj(l)*r; }
    inline void zeroit(vRealD &z){ vzero(z);}

    inline vRealD outerProduct(const vRealD &l, const vRealD& r)
    {
        return l*r;
    }


}
#endif
