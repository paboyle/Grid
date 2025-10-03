/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/WilsonKernels.cc

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/WilsonFermion.h>

NAMESPACE_BEGIN(Grid);

const std::vector<int> WilsonFermionStatic::directions(WilsonFermionStatic::MakeDirections());
const std::vector<int> WilsonFermionStatic::displacements(WilsonFermionStatic::MakeDisplacements());
int WilsonFermionStatic::HandOptDslash;
std::vector<int> WilsonFermionStatic::MakeDirections (void)
{
  std::vector<int> directions(2*Nd);
  for(int d=0;d<Nd;d++){
    directions[d]    = d;
    directions[d+Nd] = d;
  }
  return directions;
}
std::vector<int> WilsonFermionStatic::MakeDisplacements(void)
{
  std::vector<int> displacements(2*Nd);
  for(int d=0;d<Nd;d++){
    displacements[d]    = +1;
    displacements[d+Nd] = -1;
  }
  return displacements;
}

NAMESPACE_END(Grid);

