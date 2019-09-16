/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/ImprovedStaggeredFermion5D.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/ImprovedStaggeredFermion5D.h>

NAMESPACE_BEGIN(Grid);
  
// S-direction is INNERMOST and takes no part in the parity.
const std::vector<int> ImprovedStaggeredFermion5DStatic::directions({1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4});
const std::vector<int> ImprovedStaggeredFermion5DStatic::displacements({1, 1, 1, 1, -1, -1, -1, -1, 3, 3, 3, 3, -3, -3, -3, -3});
  
NAMESPACE_END(Grid);



