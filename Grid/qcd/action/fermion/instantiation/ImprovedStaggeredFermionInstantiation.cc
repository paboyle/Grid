/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/ImprovedStaggeredFermion.cc

Copyright (C) 2015

Author: Azusa Yamaguchi, Peter Boyle

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
#include <Grid/Grid.h>

NAMESPACE_BEGIN(Grid);

const std::vector<int> ImprovedStaggeredFermionStatic::directions({0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3});
const std::vector<int> ImprovedStaggeredFermionStatic::displacements({1, 1, 1, 1, -1, -1, -1, -1, 3, 3, 3, 3, -3, -3, -3, -3});
std::vector<int> ImprovedStaggeredFermionStatic::MakeDirections(void)
{
  std::vector<int> directions(4*Nd);
  for(int d=0;d<Nd;d++){
    directions[d+Nd*0] = d;
    directions[d+Nd*1] = d;
    directions[d+Nd*2] = d;
    directions[d+Nd*3] = d;
  }
  return directions;
}
std::vector<int> ImprovedStaggeredFermionStatic::MakeDisplacements(void)
{
  std::vector<int> displacements(4*Nd);
  for(int d=0;d<Nd;d++){
    displacements[d+Nd*0] =+1;
    displacements[d+Nd*1] =-1;
    displacements[d+Nd*2] =+3;
    displacements[d+Nd*3] =-3;
  }
  return displacements;
}
NAMESPACE_END(Grid);
