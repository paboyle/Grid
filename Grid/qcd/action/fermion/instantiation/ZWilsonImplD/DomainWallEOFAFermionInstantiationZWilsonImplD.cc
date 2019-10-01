/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DomainWallEOFAFermion.cc

Copyright (C) 2017

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

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

#include <Grid/Grid_Eigen_Dense.h>
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/DomainWallEOFAFermion.h>
#include <Grid/qcd/action/fermion/implementation/DomainWallEOFAFermionImplementation.h>
#include <Grid/qcd/action/fermion/implementation/DomainWallEOFAFermionCache.h>

NAMESPACE_BEGIN(Grid);

#include "impl.h"
template class DomainWallEOFAFermion<IMPLEMENTATION>; 

NAMESPACE_END(Grid);
