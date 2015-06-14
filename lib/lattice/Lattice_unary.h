#ifndef GRID_LATTICE_UNARY_H
#define GRID_LATTICE_UNARY_H

namespace Grid {

  template<class obj> Lattice<obj> pow(const Lattice<obj> &rhs,RealD y){
    Lattice<obj> ret(rhs._grid);
    ret.checkerboard = rhs.checkerboard;
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      ret._odata[ss]=pow(rhs._odata[ss],y);
    }
    return ret;
  }
  template<class obj> Lattice<obj> mod(const Lattice<obj> &rhs,Integer y){
    Lattice<obj> ret(rhs._grid);
    ret.checkerboard = rhs.checkerboard;
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      ret._odata[ss]=mod(rhs._odata[ss],y);
    }
    return ret;
  }

}
#endif
