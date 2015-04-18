#ifndef GRID_LATTICE_CONFORMABLE_H
#define GRID_LATTICE_CONFORMABLE_H

namespace Grid {

    template<class obj1,class obj2>
    void conformable(const Lattice<obj1> &lhs,const Lattice<obj2> &rhs)
    {
        assert(lhs._grid == rhs._grid);
        assert(lhs.checkerboard == rhs.checkerboard);
    }

}
#endif
