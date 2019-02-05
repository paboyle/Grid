/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Archive/Modules/WeakHamiltonianNonEye.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MContraction_WeakHamiltonianNonEye_hpp_
#define Hadrons_MContraction_WeakHamiltonianNonEye_hpp_

#include <Hadrons/Modules/MContraction/WeakHamiltonian.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         WeakHamiltonianNonEye                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

enum
{
    W_diag = 0,
    C_diag = 1,
    n_noneye_diag = 2
};

// Wing and Connected subdiagram contractions
#define MAKE_CW_SUBDIAG(Q_1, Q_2, gamma) (Q_1*adj(Q_2)*g5*gamma)

MAKE_WEAK_MODULE(WeakHamiltonianNonEye)

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakHamiltonianNonEye_hpp_
