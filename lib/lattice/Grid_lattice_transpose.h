#ifndef GRID_LATTICE_TRANSPOSE_H
#define GRID_LATTICE_TRANSPOSE_H

///////////////////////////////////////////////
// Transpose
///////////////////////////////////////////////

namespace Grid {

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Transpose
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class vobj>
    inline Lattice<vobj> transpose(const Lattice<vobj> &lhs){
        Lattice<vobj> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = transpose(lhs._odata[ss]);
        }
        return ret;
    };
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Index level dependent transpose
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<int Index,class vobj>
    inline auto transposeIndex(const Lattice<vobj> &lhs)
      -> Lattice<decltype(transposeIndex<Index>(lhs._odata[0]))>
    {
      Lattice<decltype(transposeIndex<Index>(lhs._odata[0]))> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = transposeIndex<Index>(lhs._odata[ss]);
        }
        return ret;
    };

}
#endif
