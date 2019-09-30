    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/Fermion_base_aggregate.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>

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
#ifndef  GRID_QCD_FERMION_CORE_H
#define  GRID_QCD_FERMION_CORE_H

#include <Grid/GridCore.h>
#include <Grid/GridQCDcore.h>
#include <Grid/qcd/action/ActionCore.h>

////////////////////////////////////////////
// Fermion prereqs
////////////////////////////////////////////
#include <Grid/qcd/action/fermion/WilsonCompressor.h>     //used by all wilson type fermions
NAMESPACE_CHECK(Compressor);
#include <Grid/qcd/action/fermion/FermionOperatorImpl.h>
NAMESPACE_CHECK(FermionOperatorImpl);
#include <Grid/qcd/action/fermion/FermionOperator.h>
NAMESPACE_CHECK(FermionOperator);
#include <Grid/qcd/action/fermion/WilsonKernels.h>        //used by all wilson type fermions
#include <Grid/qcd/action/fermion/StaggeredKernels.h>        //used by all wilson type fermions
NAMESPACE_CHECK(Kernels);

#endif
