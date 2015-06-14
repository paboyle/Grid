#ifndef GRID_LATTICE_REALITY_H
#define GRID_LATTICE_REALITY_H


// FIXME .. this is the sector of the code 
// I am most worried about the directions
// The choice of burying complex in the SIMD
// is making the use of "real" and "imag" very cumbersome

namespace Grid {

    template<class vobj> inline Lattice<vobj> adj(const Lattice<vobj> &lhs){
        Lattice<vobj> ret(lhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = adj(lhs._odata[ss]);
        }
        return ret;
    };

    template<class vobj> inline Lattice<vobj> conjugate(const Lattice<vobj> &lhs){
        Lattice<vobj> ret(lhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = conjugate(lhs._odata[ss]);
        }
        return ret;
    };


}
#endif
