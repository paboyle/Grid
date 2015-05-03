#ifndef GRID_VREALF_H
#define GRID_VREALF_H


namespace Grid {
    class vRealF  {
    public:
        fvec v;

    public:
	typedef fvec  vector_type;
	typedef RealF scalar_type;

        vRealF()=default;
        vRealF(RealF a){
	  vsplat(*this,a);
	};
        vRealF(Zero &zero){
	  zeroit(*this);
	}
        ////////////////////////////////////
        // Arithmetic operator overloads +,-,*
        ////////////////////////////////////
        friend inline vRealF operator + ( vRealF a,  vRealF b)
        {
            vRealF ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_add_ps(a.v,b.v);
#endif
#ifdef SSE4
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
#ifdef SSE4
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
#ifdef SSE4
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
#if defined (AVX1) || defined (SSE4)
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


	////////////////////////////////////////////////////////////////////
	// General permute; assumes vector length is same across 
	// all subtypes; may not be a good assumption, but could
	// add the vector width as a template param for BG/Q for example
	////////////////////////////////////////////////////////////////////
	/*
	friend inline void permute(vRealF &y,vRealF b,int perm)
	{
	  Gpermute<vRealF>(y,b,perm);
	}
	friend inline void merge(vRealF &y,std::vector<RealF *> &extracted)
	{
	  Gmerge<vRealF,RealF >(y,extracted);
	}
	friend inline void extract(const vRealF &y,std::vector<RealF *> &extracted)
	{
	  Gextract<vRealF,RealF>(y,extracted);
	}
	friend inline void merge(vRealF &y,std::vector<RealF> &extracted)
	{
	  Gmerge<vRealF,RealF >(y,extracted);
	}
	friend inline void extract(const vRealF &y,std::vector<RealF> &extracted)
	{
	  Gextract<vRealF,RealF>(y,extracted);
	}
	*/

        
        /////////////////////////////////////////////////////
        // Broadcast a value across Nsimd copies.
        /////////////////////////////////////////////////////
        friend inline void vsplat(vRealF &ret,float a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_ps(a,a,a,a,a,a,a,a);
#endif
#ifdef SSE4
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
#ifdef SSE4
            ret.v = _mm_set_ps(a[3],a[2],a[1],a[0]);
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

	////////////////////////////////////////////////////////////////////////
	// FIXME:  gonna remove these load/store, get, set, prefetch
	////////////////////////////////////////////////////////////////////////
friend inline void vstore(const vRealF &ret, float *a){
#if defined (AVX1)|| defined (AVX2)
	_mm256_store_ps(a,ret.v);
#endif
#ifdef SSE4
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
#if defined (SSE4)
	 // FIXME Hack
	 const RealF * ptr = (const RealF *) &in;
	 RealF ret = 0; 
	 for(int i=0;i<vRealF::Nsimd();i++){
	   ret = ret+ptr[i];
	 }
	 return ret;
#endif
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
    inline vRealF innerProduct(const vRealF & l, const vRealF & r) { return conj(l)*r; }
    inline void  zeroit(vRealF &z){ vzero(z);}

    inline vRealF outerProduct(const vRealF &l, const vRealF& r)
    {
        return l*r;
    }
    inline vRealF trace(const vRealF &arg){
        return arg;
    }
    inline vRealF real(const vRealF &arg){
        return arg;
    }

    
}
#endif
