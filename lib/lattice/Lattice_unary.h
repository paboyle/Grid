#ifndef GRID_LATTICE_UNARY_H
#define GRID_LATTICE_UNARY_H

namespace Grid {

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  avoid copy back routines for mult, mac, sub, add
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class obj> Lattice<obj> sqrt(const Lattice<obj> &rhs){
    Lattice<obj> ret(rhs._grid);
    ret.checkerboard = rhs.checkerboard;
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      ret._odata[ss]=sqrt(rhs._odata[ss]);
    }
    return ret;
  }
  template<class obj> Lattice<obj> rsqrt(const Lattice<obj> &rhs){
    Lattice<obj> ret(rhs._grid);
    ret.checkerboard = rhs.checkerboard;
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      ret._odata[ss]=rsqrt(rhs._odata[ss]);
    }
    return ret;
  }


}
#endif
