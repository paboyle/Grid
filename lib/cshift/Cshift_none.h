#ifndef _GRID_CSHIFT_NONE_H_
#define _GRID_CSHIFT_NONE_H_
namespace Grid {
template<class vobj> Lattice<vobj> Cshift(const Lattice<vobj> &rhs,int dimension,int shift)
{
  Lattice<vobj> ret(rhs._grid);
  ret.checkerboard = rhs._grid->CheckerBoardDestination(rhs.checkerboard,shift,dimension);
  Cshift_local(ret,rhs,dimension,shift);
  return ret;
}
}
#endif
