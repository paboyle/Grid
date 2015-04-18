#ifndef GRID_LATTICE_LOCALREDUCTION_H
#define GRID_LATTICE_LOCALREDUCTION_H

///////////////////////////////////////////////
// localInner, localNorm, outerProduct
///////////////////////////////////////////////

namespace Grid {

    /////////////////////////////////////////////////////
    // Non site, reduced locally reduced routines
    /////////////////////////////////////////////////////

    // localNorm2,
    template<class vobj>
    inline auto localNorm2 (const Lattice<vobj> &rhs)-> Lattice<typename vobj::tensor_reduced>
    {
      Lattice<typename vobj::tensor_reduced> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  ret._odata[ss]=innerProduct(rhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
    }
    
    // localInnerProduct
    template<class vobj>
    inline auto localInnerProduct (const Lattice<vobj> &lhs,const Lattice<vobj> &rhs)
      -> Lattice<typename vobj::tensor_reduced>
    {
      Lattice<typename vobj::tensor_reduced> ret(rhs._grid);
#pragma omp parallel for
      for(int ss=0;ss<rhs._grid->oSites(); ss++){
	ret._odata[ss]=innerProduct(lhs._odata[ss],rhs._odata[ss]);
      }
      return ret;
    }
    
    // outerProduct Scalar x Scalar -> Scalar
    //              Vector x Vector -> Matrix
    template<class ll,class rr>
    inline auto outerProduct (const Lattice<ll> &lhs,const Lattice<rr> &rhs) -> Lattice<decltype(outerProduct(lhs._odata[0],rhs._odata[0]))>
    {
        Lattice<decltype(outerProduct(lhs._odata[0],rhs._odata[0]))> ret(rhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
            ret._odata[ss]=outerProduct(lhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
     }

}

#endif
