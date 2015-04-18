#ifndef GRID_COMPARISON_H
#define GRID_COMPARISON_H
namespace Grid {

    // Generic list of functors
    template<class lobj,class robj> class veq {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs == rhs;
	}
    };
    template<class lobj,class robj> class vne {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs != rhs;
	}
    };
    template<class lobj,class robj> class vlt {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs < rhs;
	}
    };
    template<class lobj,class robj> class vle {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs <= rhs;
	}
    };
    template<class lobj,class robj> class vgt {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs > rhs;
	}
    };
    template<class lobj,class robj> class vge {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs >= rhs;
	}
    };

    // Generic list of functors
    template<class lobj,class robj> class seq {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs == rhs;
	}
    };
    template<class lobj,class robj> class sne {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs != rhs;
	}
    };
    template<class lobj,class robj> class slt {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs < rhs;
	}
    };
    template<class lobj,class robj> class sle {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs <= rhs;
	}
    };
    template<class lobj,class robj> class sgt {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs > rhs;
	}
    };
    template<class lobj,class robj> class sge {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return lhs >= rhs;
	}
    };


    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Integer gets extra relational functions. Could also implement these for RealF, RealD etc..
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class sfunctor> 
    inline vInteger Comparison(sfunctor sop,const vInteger & lhs, const vInteger & rhs)
    {
      std::vector<Integer> vlhs(vInteger::Nsimd());   // Use functors to reduce this to single implementation
      std::vector<Integer> vrhs(vInteger::Nsimd());
      vInteger ret;
      extract(lhs,vlhs);
      extract(rhs,vrhs);
      for(int s=0;s<vInteger::Nsimd();s++){
	vlhs[s] = sop(vlhs[s],vrhs[s]);
      }
      merge(ret,vlhs);
      return ret;
    }
    inline vInteger operator < (const vInteger & lhs, const vInteger & rhs)
    {
      return Comparison(slt<Integer,Integer>(),lhs,rhs);
    }
    inline vInteger operator <= (const vInteger & lhs, const vInteger & rhs)
    {
      return Comparison(sle<Integer,Integer>(),lhs,rhs);
    }
    inline vInteger operator > (const vInteger & lhs, const vInteger & rhs)
    {
      return Comparison(sgt<Integer,Integer>(),lhs,rhs);
    }
    inline vInteger operator >= (const vInteger & lhs, const vInteger & rhs)
    {
      return Comparison(sge<Integer,Integer>(),lhs,rhs);
    }
    inline vInteger operator == (const vInteger & lhs, const vInteger & rhs)
    {
      return Comparison(seq<Integer,Integer>(),lhs,rhs);
    }
    inline vInteger operator != (const vInteger & lhs, const vInteger & rhs)
    {
      return Comparison(sne<Integer,Integer>(),lhs,rhs);
    }

    //////////////////////////////////////////////////////////////////////////
    // relational operators
    // 
    // Support <,>,<=,>=,==,!=
    //
    //Query supporting bitwise &, |, ^, !
    //Query supporting logical &&, ||, 
    //////////////////////////////////////////////////////////////////////////
    template<class vfunctor,class lobj,class robj> 
    inline Lattice<vInteger> LLComparison(vfunctor op,const Lattice<lobj> &lhs,const Lattice<robj> &rhs)
    {
      Lattice<vInteger> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  ret._odata[ss]=op(lhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
    }
    template<class vfunctor,class lobj,class robj> 
    inline Lattice<vInteger> LSComparison(vfunctor op,const Lattice<lobj> &lhs,const robj &rhs)
    {
      Lattice<vInteger> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites(); ss++){
	  ret._odata[ss]=op(lhs._odata[ss],rhs);
        }
        return ret;
    }
    template<class vfunctor,class lobj,class robj> 
    inline Lattice<vInteger> SLComparison(vfunctor op,const lobj &lhs,const Lattice<robj> &rhs)
    {
      Lattice<vInteger> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  ret._odata[ss]=op(lhs._odata[ss],rhs);
        }
        return ret;
    }

    // Less than
   template<class lobj,class robj>
   inline Lattice<vInteger> operator < (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
     return LLComparison(vlt<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator < (const Lattice<lobj> & lhs, const robj & rhs) {
     return LSComparison(vlt<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator < (const lobj & lhs, const Lattice<robj> & rhs) {
     return SLComparison(vlt<lobj,robj>(),lhs,rhs);
   }

   // Less than equal
   template<class lobj,class robj>
   inline Lattice<vInteger> operator <= (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
     return LLComparison(vle<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator <= (const Lattice<lobj> & lhs, const robj & rhs) {
     return LSComparison(vle<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator <= (const lobj & lhs, const Lattice<robj> & rhs) {
     return SLComparison(vle<lobj,robj>(),lhs,rhs);
   }

   // Greater than 
   template<class lobj,class robj>
   inline Lattice<vInteger> operator > (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
     return LLComparison(vgt<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator > (const Lattice<lobj> & lhs, const robj & rhs) {
     return LSComparison(vgt<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator > (const lobj & lhs, const Lattice<robj> & rhs) {
     return SLComparison(vgt<lobj,robj>(),lhs,rhs);
   }


   // Greater than equal
   template<class lobj,class robj>
   inline Lattice<vInteger> operator >= (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
     return LLComparison(vge<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator >= (const Lattice<lobj> & lhs, const robj & rhs) {
     return LSComparison(vge<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator >= (const lobj & lhs, const Lattice<robj> & rhs) {
     return SLComparison(vge<lobj,robj>(),lhs,rhs);
   }


   // equal
   template<class lobj,class robj>
   inline Lattice<vInteger> operator == (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
     return LLComparison(veq<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator == (const Lattice<lobj> & lhs, const robj & rhs) {
     return LSComparison(veq<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator == (const lobj & lhs, const Lattice<robj> & rhs) {
     return SLComparison(veq<lobj,robj>(),lhs,rhs);
   }


   // not equal
   template<class lobj,class robj>
   inline Lattice<vInteger> operator != (const Lattice<lobj> & lhs, const Lattice<robj> & rhs) {
     return LLComparison(vne<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator != (const Lattice<lobj> & lhs, const robj & rhs) {
     return LSComparison(vne<lobj,robj>(),lhs,rhs);
   }
   template<class lobj,class robj>
   inline Lattice<vInteger> operator != (const lobj & lhs, const Lattice<robj> & rhs) {
     return SLComparison(vne<lobj,robj>(),lhs,rhs);
   }


}
#endif
