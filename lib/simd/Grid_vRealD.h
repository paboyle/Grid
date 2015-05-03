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

        friend inline void mult(vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) * (*r);}
        friend inline void sub (vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) - (*r);}
        friend inline void add (vRealD * __restrict__ y,const vRealD * __restrict__ l,const vRealD *__restrict__ r) {*y = (*l) + (*r);}
        friend inline vRealD adj(const vRealD &in) { return in; }
        friend inline vRealD conj(const vRealD &in){ return in; }
        
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
	/*
	friend inline void permute(vRealD &y,vRealD b,int perm)
	{
	  Gpermute<vRealD>(y,b,perm);
	}
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
#if defined (SSE4)
	 // FIXME Hack
	 const RealD * ptr =(const RealD *)  &in;
	 RealD ret = 0; 
	 for(int i=0;i<vRealD::Nsimd();i++){
	   ret = ret+ptr[i];
	 }
	 return ret;
#endif
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

    inline vRealD innerProduct(const vRealD & l, const vRealD & r) { return conj(l)*r; }
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
