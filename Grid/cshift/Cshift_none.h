    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cshift/Cshift_none.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
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
