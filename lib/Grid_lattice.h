#ifndef GRID_LATTICE_H
#define GRID_LATTICE_H

namespace Grid {

// TODO: 
//       mac,real,imag

// Functionality:
//     -=,+=,*=,()
//     add,+,sub,-,mult,mac,*
//     adj,conj
//     real,imag
//     transpose,transposeIndex  
//     trace,traceIndex
//     peekIndex
//     innerProduct,outerProduct,
//     localNorm2
//     localInnerProduct

extern int GridCshiftPermuteMap[4][16];

template<class vobj>
class Lattice
{
public:
    GridBase *_grid;
    int checkerboard;
    std::vector<vobj,alignedAllocator<vobj> > _odata;
public:

    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;

    //////////////////////////////////////////////////////////////////
    // Constructor requires "grid" passed.
    // what about a default grid?
    //////////////////////////////////////////////////////////////////
    Lattice(GridBase *grid) : _grid(grid) {
        _odata.reserve(_grid->oSites());
        assert((((uint64_t)&_odata[0])&0xF) ==0);
        checkerboard=0;
    }

    template<class sobj> inline Lattice<vobj> & operator = (const sobj & r){
#pragma omp parallel for
        for(int ss=0;ss<_grid->oSites();ss++){
            this->_odata[ss]=r;
        }
        return *this;
    }

    // *=,+=,-= operators inherit behvour from correspond */+/- operation
    template<class T> inline Lattice<vobj> &operator *=(const T &r) {
        *this = (*this)*r;
        return *this;
    }

    template<class T> inline Lattice<vobj> &operator -=(const T &r) {
        *this = (*this)-r;
        return *this;
    }
    template<class T> inline Lattice<vobj> &operator +=(const T &r) {
        *this = (*this)+r;
        return *this;
    }
    
    inline friend Lattice<vobj> operator / (const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
        conformable(lhs,rhs);
        Lattice<vobj> ret(lhs._grid);
#pragma omp parallel for
        for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = lhs._odata[ss]/rhs._odata[ss];
        }
        return ret;
    };
 }; // class Lattice
}

#include <lattice/Grid_lattice_conformable.h>
#include <lattice/Grid_lattice_arith.h>
#include <lattice/Grid_lattice_trace.h>
#include <lattice/Grid_lattice_transpose.h>
#include <lattice/Grid_lattice_local.h>
#include <lattice/Grid_lattice_reduction.h>
#include <lattice/Grid_lattice_peekpoke.h>
#include <lattice/Grid_lattice_reality.h>
#include <lattice/Grid_lattice_coordinate.h>
#include <lattice/Grid_lattice_rng.h>
#include <lattice/Grid_lattice_transfer.h>



#endif
