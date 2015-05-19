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
        vRealF & operator = ( Zero & z){
	  vzero(*this);
	  return (*this);
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
        // mult, sub, add, adj,conjugate, mac functions
        ///////////////////////////////////////////////
        friend inline void mult(vRealF * __restrict__ y,const vRealF * __restrict__ l,const vRealF *__restrict__ r) {*y = (*l) * (*r);}
        friend inline void sub (vRealF * __restrict__ y,const vRealF * __restrict__ l,const vRealF *__restrict__ r) {*y = (*l) - (*r);}
        friend inline void add (vRealF * __restrict__ y,const vRealF * __restrict__ l,const vRealF *__restrict__ r) {*y = (*l) + (*r);}
        friend inline vRealF adj(const vRealF &in) { return in; }
        friend inline vRealF conjugate(const vRealF &in){ return in; }

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

	friend inline void permute(vRealF &y,vRealF b,int perm)
	{
	  Gpermute<vRealF>(y,b,perm);
	}
	/*
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
        friend inline void vstream(vRealF &out,const vRealF &in){
#if defined (AVX1)|| defined (AVX2)
	  _mm256_stream_ps((float *)&out.v,in.v);
#endif
#ifdef SSE4
	  _mm_stream_ps((float *)&out.v,in.v);
#endif
#ifdef AVX512
	  _mm512_storenrngo_ps((float *)&out.v,in.v);
	  //	  _mm512_stream_ps((float *)&out.v,in.v);
	  //Note v has a3 a2 a1 a0
#endif
#ifdef QPX
	  assert(0);
#endif
	}


        friend inline void prefetch(const vRealF &v)
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
	 vRealF v1,v2;
	 union { 
	   fvec v;
	   float f[sizeof(fvec)/sizeof(double)];
	 } conv;
#ifdef SSE4
	 permute(v1,in,0); // sse 128; quad single
	 v1=v1+in;
	 permute(v2,v1,1); 
	 v1=v1+v2;
#endif
#if defined(AVX1) || defined (AVX2)
	 permute(v1,in,0); // avx 256; octo-double
	 v1=v1+in;
	 permute(v2,v1,1); 
	 v1=v1+v2;
	 permute(v2,v1,2); 
	 v1=v1+v2;
#endif
#ifdef AVX512
	 permute(v1,in,0); // avx 256; octo-double
	 v1=v1+in;
	 permute(v2,v1,1); 
	 v1=v1+v2;
	 permute(v2,v1,2); 
	 v1=v1+v2;
	 permute(v2,v1,3); 
	 v1=v1+v2;
#endif
#ifdef QPX
#endif
	 conv.v=v1.v;
	 return conv.f[0];
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
    inline vRealF innerProduct(const vRealF & l, const vRealF & r) { return conjugate(l)*r; }
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
