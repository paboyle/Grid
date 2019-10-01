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
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsAsmImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsHandImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsHandGparityImplementation.h>

NAMESPACE_BEGIN(Grid);

// Move these

#include "impl.h"

// G-parity requires more specialised implementation.
template <>
void WilsonKernels<IMPLEMENTATION>::ContractConservedCurrentSiteFwd(const SitePropagator &q_in_1, 
							  const SitePropagator &q_in_2, 
							  SitePropagator &q_out, 
							  DoubledGaugeFieldView &U,
							  unsigned int sU, 
							  unsigned int mu, 
							  bool switch_sign)
{			
  assert(0);		
}			
template <>		
void WilsonKernels<IMPLEMENTATION>::ContractConservedCurrentSiteBwd( const SitePropagator &q_in_1,
							   const SitePropagator &q_in_2,
							   SitePropagator &q_out,
							   DoubledGaugeFieldView &U,
							   unsigned int mu,
							   unsigned int sU,
							   bool switch_sign)
{			
  assert(0);		
}

HAND_SPECIALISE_GPARITY(IMPLEMENTATION);


template class WilsonKernels<IMPLEMENTATION>; 


NAMESPACE_END(Grid);

