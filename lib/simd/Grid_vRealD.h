#ifndef GRID_VREALD_H
#define GRID_VREALD_H

namespace Grid {
    class vRealD  {
    public:
        dvec v; // dvec is double precision vector

    public:
	typedef dvec  vector_type;
	typedef RealD scalar_type;

        vRealD()=default;
        vRealD(RealD a){
	  vsplat(*this,a);
	};
        vRealD(Zero &zero){
	  zeroit(*this);
	}
        vRealD & operator = ( Zero & z){
	  vzero(*this);
	  return (*this);
        }

        friend inline void mult(vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) * (*r);}
        friend inline void sub (vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) - (*r);}
        friend inline void add (vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) + (*r);}
        friend inline vRealD adj(const vRealD &in) { return in; }
        friend inline vRealD conjugate(const vRealD &in){ return in; }
        
        friend inline void mac (vRealD &y,const vRealD a,const vRealD x){
#if defined (AVX1) || defined (SSE4)
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
        friend inline vRealD operator - (vRealD a, vRealD b)
        {
            vRealD ret;
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
        
        friend inline vRealD operator * (vRealD a, vRealD b)
        {
            vRealD ret;
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_mul_pd(a.v,b.v);
#endif
#ifdef SSE4
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

	////////////////////////////////////////////////////////////////////
	// General permute; assumes vector length is same across 
	// all subtypes; may not be a good assumption, but could
	// add the vector width as a template param for BG/Q for example
	////////////////////////////////////////////////////////////////////

	friend inline void permute(vRealD &y,vRealD b,int perm)
	{
	  Gpermute<vRealD>(y,b,perm);
	}
	/*
	friend inline void merge(vRealD &y,std::vector<RealD *> &extracted)
	{
	  Gmerge<vRealD,RealD >(y,extracted);
	}
	friend inline void extract(const vRealD &y,std::vector<RealD *> &extracted)
	{
	  Gextract<vRealD,RealD>(y,extracted);
	}
	friend inline void merge(vRealD &y,std::vector<RealD > &extracted)
	{
	  Gmerge<vRealD,RealD >(y,extracted);
	}
	friend inline void extract(const vRealD &y,std::vector<RealD > &extracted)
	{
	  Gextract<vRealD,RealD>(y,extracted);
	}
	*/
        
        friend inline void vsplat(vRealD &ret,double a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set_pd(a,a,a,a);
#endif
#ifdef SSE4
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
#ifdef SSE4
            ret.v = _mm_set_pd(a[1],a[0]);
#endif
#ifdef AVX512
            ret.v = _mm512_set_pd(a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
            // Note v has a0 a1 a2 a3 a4 a5 a6 a7
#endif
#ifdef QPX
            ret.v = {a[0],a[1],a[2],a[3]};
#endif
	}

	friend inline void vstore(const vRealD &ret, double *a){
#if defined (AVX1)|| defined (AVX2)
            _mm256_store_pd(a,ret.v);
#endif
#ifdef SSE4
            _mm_store_pd(a,ret.v);
#endif
#ifdef AVX512
            _mm512_store_pd(a,ret.v);
            // Note v has a7 a6 a5ba4 a3 a2 a1 a0
#endif
#ifdef QPX
	    assert(0);
#endif
	}
        friend inline void vstream(vRealD &out,const vRealD &in){
#if defined (AVX1)|| defined (AVX2)
	  _mm256_stream_pd((double *)&out.v,in.v);
#endif
#ifdef SSE4
	  _mm_stream_pd((double *)&out.v,in.v);
#endif
#ifdef AVX512
	  _mm512_storenrngo_pd((double *)&out.v,in.v);
	  //Note v has a3 a2 a1 a0
#endif
#ifdef QPX
	  assert(0);
#endif
	}
        friend inline void prefetch(const vRealD &v)
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
	 vRealD v1,v2;
	 union { 
	   dvec v;
	   double f[sizeof(dvec)/sizeof(double)];
	 } conv;
#ifdef SSE4
	 permute(v1,in,0); // sse 128; paired real double
	 v1=v1+in;
#endif
#if defined(AVX1) || defined (AVX2)
	 permute(v1,in,0); // avx 256; quad double
	 v1=v1+in;
	 permute(v2,v1,1); 
	 v1=v1+v2;
#endif
#ifdef AVX512
	 permute(v1,in,0); // avx 512; octo-double
	 v1=v1+in;
	 permute(v2,v1,1); 
	 v1=v1+v2;
	 permute(v2,v1,2); 
	 v1=v1+v2;
#endif
#ifdef QPX
#endif
	 conv.v=v1.v;
	 return conv.f[0];
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

    inline vRealD innerProduct(const vRealD & l, const vRealD & r) { return conjugate(l)*r; }
    inline void zeroit(vRealD &z){ vzero(z);}

    inline vRealD outerProduct(const vRealD &l, const vRealD& r)
    {
        return l*r;
    }
    inline vRealD trace(const vRealD &arg){
        return arg;
    }
    inline vRealD real(const vRealD &arg){
        return arg;
    }


}
#endif
