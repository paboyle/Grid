#ifndef VREALF_H
#define VREALF_H

#include "Grid.h"

namespace Grid {
    class vRealF  {
    protected:
        fvec v;

    public:

	typedef fvec  vector_type;
	typedef RealF scalar_type;

        vRealF(){};
        ////////////////////////////////////
        // Arithmetic operator overloads +,-,*
        ////////////////////////////////////
        friend inline vRealF operator + ( vRealF a,  vRealF b)
        {
            vRealF ret;
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
            vector4double aa,bb,cc;
            aa = vec_lda(0,(float *)&a);
            bb = vec_lda(0,(float *)&b);
            cc = vec_add(aa,bb);
            vec_sta(cc,0,(float *)&ret.v);
#endif
            return ret;
        };
        
        friend inline vRealF operator - ( vRealF a, vRealF b)
        {
            vRealF ret;
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
            vector4double aa,bb,cc;
            aa = vec_lda(0,(float *)&a);
            bb = vec_lda(0,(float *)&b);
            cc = vec_sub(aa,bb);
            vec_sta(cc,0,(float *)&ret.v);
#endif
            return ret;
        };

        friend inline vRealF operator * ( vRealF a, vRealF b)
        {
            vRealF ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_mul_ps(a.v,b.v);
#endif
#ifdef SSE2
            ret.v = _mm_mul_ps(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_mul_ps(a.v,b.v);
#endif
#ifdef QPX
            vector4double aa,bb,cc; // QPX single we are forced to load as this promotes single mem->double regs.
            aa = vec_lda(0,(float *)&a);
            bb = vec_lda(0,(float *)&b);
            cc = vec_mul(aa,bb);
            vec_sta(cc,0,(float *)&ret.v);
#endif
            return ret;
        };
        
        ///////////////////////////////////////////////
        // mult, sub, add, adj,conj, mac functions
        ///////////////////////////////////////////////
        friend inline void mult(vRealF * __restrict__ y,const vRealF * __restrict__ l,const vRealF *__restrict__ r) {*y = (*l) * (*r);}
        friend inline void sub (vRealF * __restrict__ y,const vRealF * __restrict__ l,const vRealF *__restrict__ r) {*y = (*l) - (*r);}
        friend inline void add (vRealF * __restrict__ y,const vRealF * __restrict__ l,const vRealF *__restrict__ r) {*y = (*l) + (*r);}
        friend inline vRealF adj(const vRealF &in) { return in; }
        friend inline vRealF conj(const vRealF &in){ return in; }

        friend inline void mac (vRealF &y,const vRealF a,const vRealF x){
#if defined (AVX1) || defined (SSE2)
            y = a*x+y;
#endif
#ifdef AVX2     // AVX 2 introduced FMA support. FMA4 eliminates a copy, but AVX only has FMA3
            // accelerates multiply accumulate, but not general multiply add
            y.v = _mm256_fmadd_ps(a.v,x.v,y.v);
#endif
#ifdef AVX512
            y.v = _mm512_fmadd_ps(a.v,x.v,y.v);
#endif
#ifdef QPX
            vector4double aa,xx,yy; // QPX single we are forced to load as this promotes single mem->double regs.
            aa = vec_lda(0,(float *)&a.v);
            xx = vec_lda(0,(float *)&x.v);
            yy = vec_lda(0,(float *)&y.v);
            yy = vec_madd(aa,xx,yy);
            vec_sta(yy,0,(float *)&y.v);
#endif
        }
        
        //////////////////////////////////
        // Initialise to 1,0,i
        //////////////////////////////////
        friend inline void vone (vRealF &ret){vsplat(ret,1.0);}
        friend inline void vzero(vRealF &ret){vsplat(ret,0.0);}


        /////////////////////////////////////////////////////////////////
        // Extract
        /////////////////////////////////////////////////////////////////
        friend inline void extract(vRealF &y,std::vector<RealF *> &extracted){
	  // Bounce off stack is painful
	  // temporary hack while I figure out the right interface
	  const int Nsimd = vRealF::Nsimd();
	  RealF buf[Nsimd]; 

	  vstore(y,buf);

	  for(int i=0;i<Nsimd;i++){
	    *extracted[i] = buf[i];
	    extracted[i]++;
	  }
        };

        friend inline void merge(vRealF &y,std::vector<RealF *> &extracted){
	  // Bounce off stack is painful
	  // temporary hack while I figure out the right interface
	  const int Nsimd = vRealF::Nsimd();
	  RealF buf[Nsimd]; 

	  for(int i=0;i<Nsimd;i++){
	    buf[i]=*extracted[i];
	    extracted[i]++;
	  }
	  vset(y,buf); 
        };
        
        //////////////////////////////////////////////////////////
        // Permute
        // Permute 0 every ABCDEFGH -> BA DC FE HG
        // Permute 1 every ABCDEFGH -> CD AB GH EF
        // Permute 2 every ABCDEFGH -> EFGH ABCD
        // Permute 3 possible on longer iVector lengths (512bit = 8 double = 16 single)
        // Permute 4 possible on half precision @512bit vectors.
        //////////////////////////////////////////////////////////
        friend inline void permute(vRealF &y,vRealF b,int perm){
            switch (perm){
                    // 8 floats=>3 permutes
#if defined(AVX1)||defined(AVX2)
                case 0: y.v = _mm256_shuffle_ps(b.v,b.v,_MM_SHUFFLE(2,3,0,1)); break;
                case 1: y.v = _mm256_shuffle_ps(b.v,b.v,_MM_SHUFFLE(1,0,3,2)); break;
                case 2: y.v = _mm256_permute2f128_ps(b.v,b.v,0x01); break;
#endif
#ifdef SSE2
                case 0: y.v = _mm_shuffle_ps(b.v,b.v,_MM_SHUFFLE(2,3,0,1)); break;
                case 1: y.v = _mm_shuffle_ps(b.v,b.v,_MM_SHUFFLE(1,0,3,2));break;
#endif
#ifdef AVX512
                    // 16 floats=> permutes
        // Permute 0 every abcd efgh ijkl mnop -> badc fehg jilk nmpo 
        // Permute 1 every abcd efgh ijkl mnop -> cdab ghef jkij opmn 
        // Permute 2 every abcd efgh ijkl mnop -> efgh abcd mnop ijkl
        // Permute 3 every abcd efgh ijkl mnop -> ijkl mnop abcd efgh
//#error not implemented should do something
                case 0: y.v = _mm512_swizzle_ps(b.v,_MM_SWIZ_REG_CDAB); break;
                case 1: y.v = _mm512_swizzle_ps(b.v,_MM_SWIZ_REG_BADC); break;
                case 2: y.v = _mm512_permute4f128_ps(b.v,(_MM_PERM_ENUM)_MM_SHUFFLE(2,3,0,1)); break;
                case 3: y.v = _mm512_permute4f128_ps(b.v,(_MM_PERM_ENUM)_MM_SHUFFLE(1,0,3,2)); break;
#endif
#ifdef QPX
#error not implemented
#endif
	    default: assert(0); break;
            }
        };
        
        /////////////////////////////////////////////////////
        // Broadcast a value across Nsimd copies.
        /////////////////////////////////////////////////////
        friend inline void vsplat(vRealF &ret,float a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_ps(a,a,a,a,a,a,a,a);
#endif
#ifdef SSE2
            ret.v = _mm_set_ps(a,a,a,a);
#endif
#ifdef AVX512
            //ret.v = _mm512_set_ps(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a);
            ret.v = _mm512_set1_ps(a);
#endif
#ifdef QPX
            ret.v = {a,a,a,a};
#endif
        }
        friend inline void vset(vRealF &ret, float *a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_ps(a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
#endif
#ifdef SSE2
            ret.v = _mm_set_ps(a[0],a[1],a[2],a[3]);
#endif
#ifdef AVX512
            ret.v = _mm512_set_ps( a[15],a[14],a[13],a[12],a[11],a[10],a[9],a[8],
                                   a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
            // Note v has a0 a1 a2 a3 a4 a5 a6 a7
#endif
#ifdef QPX
            ret.v = {a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]};
#endif
	}

friend inline void vstore(vRealF &ret, float *a){
#if defined (AVX1)|| defined (AVX2)
	_mm256_store_ps(a,ret.v);
#endif
#ifdef SSE2
	_mm_store_ps(a,ret.v);
#endif
#ifdef AVX512
	_mm512_store_ps(a,ret.v);
	// Note v has a7 a6 a5ba4 a3 a2 a1 a0
#endif
#ifdef QPX
	assert(0);
#endif
        }


        friend inline void vprefetch(const vRealF &v)
        {
            _mm_prefetch((const char*)&v.v,_MM_HINT_T0);
        }
        // Unary negation
        friend inline vRealF operator -(const vRealF &r) {
            vRealF ret;
            vzero(ret);
            ret = ret - r;
            return ret;
        }
       friend inline RealF Reduce(const vRealF & in)
       {
#if defined (AVX1) || defined(AVX2)
            __attribute__ ((aligned(32))) float c_[16];
            __m256 tmp = _mm256_permute2f128_ps(in.v,in.v,0x01);
            __m256 hadd = _mm256_hadd_ps(in.v,tmp);
                   tmp = _mm256_permute2f128_ps(hadd,hadd,0x01);
                   hadd = _mm256_hadd_ps(tmp,tmp);
                  _mm256_store_ps(c_,hadd);
         return (float)c_[0];

#endif
#ifdef AVX512
            return _mm512_reduce_add_ps(in.v);
/*
             __attribute__ ((aligned(64))) float c_[16];
             _mm512_store_ps(c_,in.v);
             return c_[0]+c_[1]+c_[2]+c_[3]+c_[4]+c_[5]+c_[6]+c_[7]
                    +c_[8]+c_[9]+c_[10]+c_[11]+c_[12]+c_[13]+c_[14]+c_[15];
*/
#endif
#ifdef QPX
#endif
        }

        // *=,+=,-= operators
        inline vRealF &operator *=(const vRealF &r) {
            *this = (*this)*r;
            return *this;
        }
        inline vRealF &operator +=(const vRealF &r) {
            *this = *this+r;
            return *this;
        }
        inline vRealF &operator -=(const vRealF &r) {
            *this = *this-r;
            return *this;
        }
    public:
        static inline int Nsimd(void) { return sizeof(fvec)/sizeof(float);}
    };
    inline vRealF localInnerProduct(const vRealF & l, const vRealF & r) { return conj(l)*r; }
    inline void  zeroit(vRealF &z){ vzero(z);}

    inline vRealF outerProduct(const vRealF &l, const vRealF& r)
    {
        return l*r;
    }


    
}
#endif
