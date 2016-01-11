    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_comparison_utils.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_COMPARISON_H
#define GRID_COMPARISON_H

namespace Grid {

  /////////////////////////////////////////
  // This implementation is a bit poor.
  //
  // Only support relational logical operations (<, >  etc)
  // on scalar objects. Therefore can strip any tensor structures.
  //
  // Should guard this with isGridTensor<> enable if?
  /////////////////////////////////////////
  //
  // Generic list of functors
  //
  template<class lobj,class robj> class veq {
  public:
    vInteger operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) == (rhs);
    }
  };
  template<class lobj,class robj> class vne {
  public:
    vInteger operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) != (rhs);
    }
  };
  template<class lobj,class robj> class vlt {
  public:
    vInteger operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) < (rhs);
    }
  };
  template<class lobj,class robj> class vle {
  public:
    vInteger operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) <= (rhs);
    }
  };
  template<class lobj,class robj> class vgt {
  public:
    vInteger operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) > (rhs);
    }
  };
  template<class lobj,class robj> class vge {
    public:
    vInteger operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) >= (rhs);
    }
  };
  
  // Generic list of functors
  template<class lobj,class robj> class seq {
  public:
    Integer operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) == (rhs);
    }
  };
  template<class lobj,class robj> class sne {
  public:
    Integer operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) != (rhs);
    }
  };
  template<class lobj,class robj> class slt {
  public:
    Integer operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) < (rhs);
    }
  };
  template<class lobj,class robj> class sle {
  public:
    Integer operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) <= (rhs);
    }
  };
  template<class lobj,class robj> class sgt {
  public:
    Integer operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) > (rhs);
    }
  };
  template<class lobj,class robj> class sge {
  public:
    Integer operator()(const lobj &lhs, const robj &rhs)
    { 
      return (lhs) >= (rhs);
    }
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Integer and real get extra relational functions.
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class sfunctor, class vsimd,IfNotComplex<vsimd> = 0> 
    inline vInteger Comparison(sfunctor sop,const vsimd & lhs, const vsimd & rhs)
    {
      typedef typename vsimd::scalar_type scalar;
      std::vector<scalar> vlhs(vsimd::Nsimd());   // Use functors to reduce this to single implementation
      std::vector<scalar> vrhs(vsimd::Nsimd());
      std::vector<Integer> vpred(vsimd::Nsimd());
      vInteger ret;
      extract<vsimd,scalar>(lhs,vlhs);
      extract<vsimd,scalar>(rhs,vrhs);
      for(int s=0;s<vsimd::Nsimd();s++){
	vpred[s] = sop(vlhs[s],vrhs[s]);
      }
      merge<vInteger,Integer>(ret,vpred);
      return ret;
    }

  template<class sfunctor, class vsimd,IfNotComplex<vsimd> = 0> 
    inline vInteger Comparison(sfunctor sop,const vsimd & lhs, const typename vsimd::scalar_type & rhs)
    {
      typedef typename vsimd::scalar_type scalar;
      std::vector<scalar> vlhs(vsimd::Nsimd());   // Use functors to reduce this to single implementation
      std::vector<Integer> vpred(vsimd::Nsimd());
      vInteger ret;
      extract<vsimd,scalar>(lhs,vlhs);
      for(int s=0;s<vsimd::Nsimd();s++){
	vpred[s] = sop(vlhs[s],rhs);
      }
      merge<vInteger,Integer>(ret,vpred);
      return ret;
    }

  template<class sfunctor, class vsimd,IfNotComplex<vsimd> = 0> 
    inline vInteger Comparison(sfunctor sop,const typename vsimd::scalar_type & lhs, const vsimd & rhs)
    {
      typedef typename vsimd::scalar_type scalar;
      std::vector<scalar> vrhs(vsimd::Nsimd());   // Use functors to reduce this to single implementation
      std::vector<Integer> vpred(vsimd::Nsimd());
      vInteger ret;
      extract<vsimd,scalar>(rhs,vrhs);
      for(int s=0;s<vsimd::Nsimd();s++){
	vpred[s] = sop(lhs,vrhs[s]);
      }
      merge<vInteger,Integer>(ret,vpred);
      return ret;
    }

#define DECLARE_RELATIONAL(op,functor) \
  template<class vsimd,IfSimd<vsimd> = 0>\
    inline vInteger operator op (const vsimd & lhs, const vsimd & rhs)\
    {\
      typedef typename vsimd::scalar_type scalar;\
      return Comparison(functor<scalar,scalar>(),lhs,rhs);\
    }\
  template<class vsimd,IfSimd<vsimd> = 0>\
    inline vInteger operator op (const vsimd & lhs, const typename vsimd::scalar_type & rhs) \
    {\
      typedef typename vsimd::scalar_type scalar;\
      return Comparison(functor<scalar,scalar>(),lhs,rhs);\
    }\
  template<class vsimd,IfSimd<vsimd> = 0>\
    inline vInteger operator op (const typename vsimd::scalar_type & lhs, const vsimd & rhs) \
    {\
      typedef typename vsimd::scalar_type scalar;\
      return Comparison(functor<scalar,scalar>(),lhs,rhs);\
    }\
  template<class vsimd>\
    inline vInteger operator op(const iScalar<vsimd> &lhs,const iScalar<vsimd> &rhs)\
    {									\
      return lhs._internal op rhs._internal;				\
    }									\
  template<class vsimd>\
    inline vInteger operator op(const iScalar<vsimd> &lhs,const typename vsimd::scalar_type &rhs) \
    {									\
      return lhs._internal op rhs;					\
    }									\
  template<class vsimd>\
    inline vInteger operator op(const typename vsimd::scalar_type &lhs,const iScalar<vsimd> &rhs) \
    {									\
      return lhs op rhs._internal;					\
    }									


DECLARE_RELATIONAL(<,slt);
DECLARE_RELATIONAL(<=,sle);
DECLARE_RELATIONAL(>,sgt);
DECLARE_RELATIONAL(>=,sge);
DECLARE_RELATIONAL(==,seq);
DECLARE_RELATIONAL(!=,sne);

#undef DECLARE_RELATIONAL

}


#endif
