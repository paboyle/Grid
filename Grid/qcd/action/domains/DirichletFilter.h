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

#include <Grid/qcd/action/domains/DomainDecomposition.h>

NAMESPACE_BEGIN(Grid);


template<typename MomentaField>
struct DirichletFilter: public MomentumFilterBase<MomentaField>
{
  Coordinate Block;
  
  DirichletFilter(const Coordinate &_Block): Block(_Block) {}

  // Edge detect using domain projectors
  void applyFilter (MomentaField &U) const override
  {
    DomainDecomposition Domains(Block);
    GridBase *grid = U.Grid();
    LatticeInteger  coor(grid);
    LatticeInteger  face(grid);
    LatticeInteger  one(grid);   one = 1;
    LatticeInteger  zero(grid); zero = 0;
    LatticeInteger  omega(grid);
    LatticeInteger  omegabar(grid);
    LatticeInteger  tmp(grid);

    omega=one;    Domains.ProjectDomain(omega,0);
    omegabar=one; Domains.ProjectDomain(omegabar,1);
    
    LatticeInteger nface(grid); nface=Zero();
    
    MomentaField projected(grid); projected=Zero();
    typedef decltype(PeekIndex<LorentzIndex>(U,0)) MomentaLinkField;
    MomentaLinkField  Umu(grid);
    MomentaLinkField   zz(grid); zz=Zero();

    int dims = grid->Nd();
    Coordinate Global=grid->GlobalDimensions();
    assert(dims==Nd);

    for(int mu=0;mu<Nd;mu++){

      if ( Block[mu]!=0 ) {

	Umu = PeekIndex<LorentzIndex>(U,mu);

	// Upper face 
 	tmp = Cshift(omegabar,mu,1);
	tmp = tmp + omega;
	face = where(tmp == Integer(2),one,zero );

 	tmp = Cshift(omega,mu,1);
	tmp = tmp + omegabar;
	face = where(tmp == Integer(2),one,face );

	Umu = where(face,zz,Umu);

	PokeIndex<LorentzIndex>(U, Umu, mu);
      }
    }
  }
  
};

NAMESPACE_END(Grid);
