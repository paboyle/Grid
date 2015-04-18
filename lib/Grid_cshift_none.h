#ifndef _GRID_NONE_CSHIFT_H_
#define _GRID_NONE_CSHIFT_H_

friend Lattice<vobj> Cshift(Lattice<vobj> &rhs,int dimension,int shift)
{
  Lattice<vobj> ret(rhs._grid);
  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift);
  Cshift_local(ret,rhs,dimension,shift);
  return ret;
}
        
#endif
