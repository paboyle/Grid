/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/momentum/Domains.h

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


struct DomainDecomposition
{
  Coordinate Block;
  static constexpr RealD factor = 0.0;

  DomainDecomposition(const Coordinate &_Block): Block(_Block){};
  
  template<class Field>
  void ProjectDomain(Field &f,Integer domain)
  {
    GridBase *grid = f.Grid();
    int dims = grid->Nd();
    int isDWF= (dims==Nd+1);
    assert((dims==Nd)||(dims==Nd+1));

    Field   zz(grid);  zz = Zero();
    LatticeInteger coor(grid);
    LatticeInteger domaincoor(grid);
    LatticeInteger mask(grid); mask = Integer(1);
    LatticeInteger zi(grid);     zi = Integer(0);
    for(int d=0;d<Nd;d++){
      Integer B= Block[d];
      if ( B ) {
	LatticeCoordinate(coor,d+isDWF);
	domaincoor = mod(coor,B);
	mask = where(domaincoor==Integer(0),zi,mask);
	mask = where(domaincoor==Integer(B-1),zi,mask);
      }
    }
    if ( domain )
      f = where(mask==Integer(1),f,zz);
    else 
      f = where(mask==Integer(0),f,zz);
  };
  template<class GaugeField>
  void ProjectDDHMC(GaugeField &U)
  {
    GridBase *grid = U.Grid();
    GaugeField Uslow(grid);
    Coordinate Global=grid->GlobalDimensions();
    GaugeField zzz(grid); zzz = Zero();
    LatticeInteger coor(grid); 

    ////////////////////////////////////////////////////
    // All links except BDY layers
    ////////////////////////////////////////////////////
    for(int mu=0;mu<Nd;mu++) {
      Integer B1 = Block[mu];
      if ( B1 && (B1 <= Global[mu]) ) {
	LatticeCoordinate(coor,mu);

	////////////////////////////////
	// OmegaBar - zero all links contained
	////////////////////////////////
	U = where(mod(coor,B1)==Integer(B1-1),zzz,U); 

	auto zzz_mu = PeekIndex<LorentzIndex>(zzz,mu);
	auto U_mu = PeekIndex<LorentzIndex>(U,mu);
	U_mu = where(mod(coor,B1)==Integer(0),zzz_mu,U_mu); 
	U = where(mod(coor,B1)==Integer(0),zzz,U); 
	PokeIndex<LorentzIndex>(U, U_mu, mu);

	Uslow = U * factor;
	////////////////////////////////////////////
	// Omega interior slow the evolution
	////////////////////////////////////////////
	U = where(mod(coor,B1)==Integer(B1-4),Uslow,U); // Could scale down???
	U = where(mod(coor,B1)==Integer(B1-3),Uslow,U); // Could scale down???
	U = where(mod(coor,B1)==Integer(B1-2),Uslow,U); 

	U_mu = PeekIndex<LorentzIndex>(U,mu);
	U_mu = where(mod(coor,B1)==Integer(1),U_mu*factor,U_mu);
	U = where(mod(coor,B1)==Integer(1),zzz,U); 
	PokeIndex<LorentzIndex>(U, U_mu, mu);


	U = where(mod(coor,B1)==Integer(2),Uslow,U); 

	U_mu = PeekIndex<LorentzIndex>(U,mu);                  
	U = where(mod(coor,B1)==Integer(3),Uslow,U);           
	PokeIndex<LorentzIndex>(U, U_mu, mu);
      }
    }
  }
};

NAMESPACE_END(Grid);
