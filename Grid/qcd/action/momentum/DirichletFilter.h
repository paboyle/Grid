/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/momentum/DirichletFilter.h

Copyright (C) 2021

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */

////////////////////////////////////////////////////
// Dirichlet filter with sub-block size B[mu]
////////////////////////////////////////////////////
#pragma once 

NAMESPACE_BEGIN(Grid);

template<typename MomentaField>
struct DirichletFilter: public MomentumFilterBase<MomentaField>
{
  Coordinate Block;
  
  DirichletFilter(const Coordinate &_Block): Block(_Block){}

  void applyFilter(MomentaField &P) const override
  {
    GridBase *grid = P.Grid();

    ////////////////////////////////////////////////////
    // Zero strictly links crossing between domains
    ////////////////////////////////////////////////////
    LatticeInteger coor(grid);
    typedef decltype(PeekIndex<LorentzIndex>(P,0)) MatrixType;
    MatrixType zz(grid); zz = Zero();
    Coordinate Global=grid->GlobalDimensions();
    for(int mu=0;mu<Nd;mu++) {
      if ( (Block[mu] <= Global[mu]) && (Block[mu]>1) ) {
	// If costly could provide Grid earlier and precompute masks
	LatticeCoordinate(coor,mu);
      
	auto P_mu = PeekIndex<LorentzIndex>(P, mu);
	P_mu = where(mod(coor,Block[mu])==Integer(Block[mu]-1),zz,P_mu);
	PokeIndex<LorentzIndex>(P, P_mu, mu);
      }
    }
  }
};

NAMESPACE_END(Grid);
