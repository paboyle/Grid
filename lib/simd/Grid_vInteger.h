#ifndef GRID_VINTEGER_H
#define GRID_VINTEGER_H

namespace Grid {

#define _mm256_set_m128i(hi,lo) _mm256_insertf128_si256(_mm256_castsi128_si256(lo),(hi),1)
// _mm256_set_m128i(hi,lo); // not defined in all versions of immintrin.h

  typedef uint32_t Integer;

  class vInteger {
    protected:

    public:

      ivec v;

	typedef ivec     vector_type;
	typedef Integer scalar_type;

        vInteger(){};
        vInteger & operator = (const Zero & z){
            vzero(*this);
            return (*this);
        }
        vInteger(Integer a){
	  vsplat(*this,a);
	};
        ////////////////////////////////////
        // Arithmetic operator overloads +,-,*
        ////////////////////////////////////
        friend inline vInteger operator + ( vInteger a,  vInteger b)
        {
	  vInteger ret;
#if defined (AVX1) 
	  __m128i a0,a1;
	  __m128i b0,b1;
	  a0 = _mm256_extractf128_si256(a.v,0);
	  b0 = _mm256_extractf128_si256(b.v,0);
	  a1 = _mm256_extractf128_si256(a.v,1);
	  b1 = _mm256_extractf128_si256(b.v,1);
	  a0 = _mm_add_epi32(a0,b0);
	  a1 = _mm_add_epi32(a1,b1);
	  ret.v = _mm256_set_m128i(a1,a0);
#endif
#if defined (AVX2)
            ret.v = _mm256_add_epi32(a.v,b.v);
#endif
#ifdef SSE4
            ret.v = _mm_add_epi32(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_add_epi32(a.v,b.v);
#endif
#ifdef QPX
            // Implement as array of ints is only option
#error
#endif
            return ret;
        };
        
        friend inline vInteger operator - ( vInteger a, vInteger b)
        {
            vInteger ret;
#if defined (AVX1) 
	  __m128i a0,a1;
	  __m128i b0,b1;
	  a0 = _mm256_extractf128_si256(a.v,0);
	  b0 = _mm256_extractf128_si256(b.v,0);
	  a1 = _mm256_extractf128_si256(a.v,1);
	  b1 = _mm256_extractf128_si256(b.v,1);
	  a0 = _mm_sub_epi32(a0,b0);
	  a1 = _mm_sub_epi32(a1,b1);
	  ret.v = _mm256_set_m128i(a1,a0);
#endif
#if defined (AVX2)
            ret.v = _mm256_sub_epi32(a.v,b.v);
#endif
#ifdef SSE4
            ret.v = _mm_sub_epi32(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_sub_epi32(a.v,b.v);
#endif
#ifdef QPX
            // Implement as array of ints is only option
#error
#endif
            return ret;
        };

        friend inline vInteger operator * ( vInteger a, vInteger b)
        {
            vInteger ret;
#if defined (AVX1) 
	  __m128i a0,a1;
	  __m128i b0,b1;
	  a0 = _mm256_extractf128_si256(a.v,0);
	  b0 = _mm256_extractf128_si256(b.v,0);
	  a1 = _mm256_extractf128_si256(a.v,1);
	  b1 = _mm256_extractf128_si256(b.v,1);
	  a0 = _mm_mul_epi32(a0,b0);
	  a1 = _mm_mul_epi32(a1,b1);
	  ret.v = _mm256_set_m128i(a1,a0);
#endif
#if defined (AVX2)
            ret.v = _mm256_mul_epi32(a.v,b.v);
#endif
#ifdef SSE4
            ret.v = _mm_mul_epi32(a.v,b.v);
#endif
#ifdef AVX512
            ret.v = _mm512_mul_epi32(a.v,b.v);
#endif
#ifdef QPX
            // Implement as array of ints is only option
#error
#endif
            return ret;
        };
        
        ///////////////////////////////////////////////
        // mult, sub, add, adj,conj, mac functions
        ///////////////////////////////////////////////
        friend inline void mult(vInteger * __restrict__ y,const vInteger * __restrict__ l,const vInteger *__restrict__ r) {*y = (*l) * (*r);}
        friend inline void sub (vInteger * __restrict__ y,const vInteger * __restrict__ l,const vInteger *__restrict__ r) {*y = (*l) - (*r);}
        friend inline void add (vInteger * __restrict__ y,const vInteger * __restrict__ l,const vInteger *__restrict__ r) {*y = (*l) + (*r);}
        friend inline void mac (vInteger &y,const vInteger a,const vInteger x){
            y = a*x+y;
	}
        
        //////////////////////////////////
        // Initialise to 1,0,i
        //////////////////////////////////
        friend inline void vone (vInteger &ret){vsplat(ret,1);}
        friend inline void vzero(vInteger &ret){vsplat(ret,0);}
        friend inline void vtrue (vInteger &ret){vsplat(ret,0xFFFFFFFF);}
        friend inline void vfalse(vInteger &ret){vsplat(ret,0);}

        
        /////////////////////////////////////////////////////
        // Broadcast a value across Nsimd copies.
        /////////////////////////////////////////////////////
        friend inline void vsplat(vInteger &ret,scalar_type a){
#if defined (AVX1)|| defined (AVX2)
            ret.v = _mm256_set1_epi32(a);
#endif
#ifdef SSE4
            ret.v = _mm_set1_epi32(a);
#endif
#ifdef AVX512
            ret.v = _mm512_set1_epi32(a);
#endif
#ifdef QPX
#error
#endif
        }
        friend inline void vset(vInteger &ret,scalar_type *a){
#if defined (AVX1)|| defined (AVX2)
	  ret.v = _mm256_set_epi32(a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
#endif
#ifdef SSE4
	  ret.v = _mm_set_epi32(a[0],a[1],a[2],a[3]);
#endif
#ifdef AVX512
	  ret.v = _mm512_set_epi32( a[15],a[14],a[13],a[12],a[11],a[10],a[9],a[8],
                                   a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
#endif
#ifdef QPX
#error
#endif
	}

	friend inline void vstore(const vInteger &ret, Integer *a){
#if defined (AVX1)|| defined (AVX2)
	  _mm256_store_si256((__m256i*)a,ret.v);
#endif
#ifdef SSE4
	  _mm_store_si128((__m128i *)a,ret.v);
#endif
#ifdef AVX512
	  _mm512_store_si512(a,ret.v);
#endif
#ifdef QPX
	  assert(0);
#endif
        }

	friend inline void vstream(vInteger & out,const vInteger &in){
	  out=in;
        }

        friend inline void prefetch(const vInteger &v)
        {
            _mm_prefetch((const char*)&v.v,_MM_HINT_T0);
        }

        // Unary negation
        friend inline vInteger operator -(const vInteger &r) {
            vInteger ret;
            vzero(ret);
            ret = ret - r;
            return ret;
        }
       friend inline Integer Reduce(const vInteger & in)
       {
	 // unimplemented
	 assert(0);
       }
        // *=,+=,-= operators
        inline vInteger &operator *=(const vInteger &r) {
            *this = (*this)*r;
            return *this;
        }
        inline vInteger &operator +=(const vInteger &r) {
            *this = *this+r;
            return *this;
        }
        inline vInteger &operator -=(const vInteger &r) {
            *this = *this-r;
            return *this;
        }

      friend inline void permute(vInteger &y,const vInteger b,int perm)
      {
	Gpermute<vInteger>(y,b,perm);
      }
      /*
      friend inline void merge(vInteger &y,std::vector<Integer *> &extracted)
      {
	Gmerge<vInteger,Integer>(y,extracted);
      }
      friend inline void extract(const vInteger &y,std::vector<Integer *> &extracted)
      {
	Gextract<vInteger,Integer>(y,extracted);
      }
      friend inline void merge(vInteger &y,std::vector<Integer> &extracted)
      {
	Gmerge<vInteger,Integer>(y,extracted);
      }
      friend inline void extract(const vInteger &y,std::vector<Integer> &extracted)
      {
	Gextract<vInteger,Integer>(y,extracted);
      }
      */

    public:
        static inline int Nsimd(void) { return sizeof(ivec)/sizeof(Integer);}
  };

    inline vInteger localInnerProduct(const vInteger & l, const vInteger & r) { return l*r; }

    inline void  zeroit(vInteger &z){ vzero(z);}

    inline vInteger outerProduct(const vInteger &l, const vInteger& r)
    {
        return l*r;
    }
    
}

#endif
