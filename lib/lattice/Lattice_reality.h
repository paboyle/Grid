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

    template<class vobj> inline auto real(const Lattice<vobj> &z) -> Lattice<decltype(real(z._odata[0]))>
    {
      Lattice<decltype(real(z._odata[0]))> ret(z._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<z._grid->oSites();ss++){
            ret._odata[ss] = real(z._odata[ss]);
        }
      return ret;
    }

    template<class vobj> inline auto imag(const Lattice<vobj> &z) -> Lattice<decltype(imag(z._odata[0]))>
    {
      Lattice<decltype(imag(z._odata[0]))> ret(z._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<z._grid->oSites();ss++){
            ret._odata[ss] = imag(z._odata[ss]);
        }
      return ret;
    }


}
#endif
