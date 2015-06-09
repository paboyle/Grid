#ifndef GRID_COMPARISON_H
#define GRID_COMPARISON_H

namespace Grid {

  /////////////////////////////////////////
  // This implementation is a bit poor.
  // Only support logical operations (== etc)
  // on scalar objects. Strip any tensor structures.
  // Should guard this with isGridTensor<> enable if?
  /////////////////////////////////////////
    // Generic list of functors
    template<class lobj,class robj> class veq {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) == TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class vne {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) != TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class vlt {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) < TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class vle {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) <= TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class vgt {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) > TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class vge {
    public:
      vInteger operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) >= TensorRemove(rhs);
	}
    };

    // Generic list of functors
    template<class lobj,class robj> class seq {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) == TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class sne {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) != TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class slt {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) < TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class sle {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) <= TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class sgt {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) > TensorRemove(rhs);
	}
    };
    template<class lobj,class robj> class sge {
    public:
      Integer operator()(const lobj &lhs, const robj &rhs)
	{ 
	  return TensorRemove(lhs) >= TensorRemove(rhs);
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
      extract<vInteger,Integer>(lhs,vlhs);
      extract<vInteger,Integer>(rhs,vrhs);
      for(int s=0;s<vInteger::Nsimd();s++){
	vlhs[s] = sop(vlhs[s],vrhs[s]);
      }
      merge<vInteger,Integer>(ret,vlhs);
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
}


#endif
