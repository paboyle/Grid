#ifndef GRID_LATTICE_COMPARISON_H
#define GRID_LATTICE_COMPARISON_H

namespace Grid {

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
