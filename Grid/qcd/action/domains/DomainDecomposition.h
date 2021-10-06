/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/domains/DomainDecomposition.h

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
  static constexpr RealD factor = 0.6;

  DomainDecomposition(const Coordinate &_Block): Block(_Block){ assert(Block.size()==Nd);};
  
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
    if ( !domain )
      f = where(mask==Integer(1),f,zz);
    else 
      f = where(mask==Integer(0),f,zz);
  };
  template<class GaugeField>
  void ProjectDDHMC(GaugeField &U)
  {
    GridBase *grid = U.Grid();
    Coordinate Global=grid->GlobalDimensions();
    GaugeField zzz(grid); zzz = Zero();
    LatticeInteger coor(grid); 

    GaugeField Uorg(grid); Uorg = U;
    
    auto zzz_mu = PeekIndex<LorentzIndex>(zzz,0);
    ////////////////////////////////////////////////////
    // Zero BDY layers
    ////////////////////////////////////////////////////
    for(int mu=0;mu<Nd;mu++) {
      Integer B1 = Block[mu];
      if ( B1 && (B1 <= Global[mu]) ) {
	LatticeCoordinate(coor,mu);


	////////////////////////////////
	// OmegaBar - zero all links contained in slice B-1,0 and
	// mu links connecting to Omega
	////////////////////////////////

	U    = where(mod(coor,B1)==Integer(B1-1),zzz,U);
	U    = where(mod(coor,B1)==Integer(0)   ,zzz,U); 

	auto U_mu   = PeekIndex<LorentzIndex>(U,mu);
	U_mu = where(mod(coor,B1)==Integer(B1-2),zzz_mu,U_mu); 
	PokeIndex<LorentzIndex>(U, U_mu, mu);

      }
    }
   
    ////////////////////////////////////////////
    // Omega interior slow the evolution
    // Tricky as we need to take the smallest of values imposed by each cut
    // Do them in order or largest to smallest and smallest writes last
    ////////////////////////////////////////////
    RealD f= factor;
#if 0    
    for(int mu=0;mu<Nd;mu++) {
      Integer B1 = Block[mu];
      if ( B1 && (B1 <= Global[mu]) ) {

	auto U_mu   = PeekIndex<LorentzIndex>(U,mu);
	auto Uorg_mu= PeekIndex<LorentzIndex>(Uorg,mu);
	// In the plane
	U = where(mod(coor,B1)==Integer(B1-5),Uorg*f,U); 
	U = where(mod(coor,B1)==Integer(4)   ,Uorg*f,U); 

	// Perp links
       	U_mu = where(mod(coor,B1)==Integer(B1-6),Uorg_mu*f,U_mu);
	U_mu = where(mod(coor,B1)==Integer(4)   ,Uorg_mu*f,U_mu);

	PokeIndex<LorentzIndex>(U, U_mu, mu);
      }
    }
#endif
    for(int mu=0;mu<Nd;mu++) {
      Integer B1 = Block[mu];
      if ( B1 && (B1 <= Global[mu]) ) {

	auto U_mu   = PeekIndex<LorentzIndex>(U,mu);
	auto Uorg_mu= PeekIndex<LorentzIndex>(Uorg,mu);
	// In the plane
	U = where(mod(coor,B1)==Integer(B1-4),Uorg*f*f,U); 
	U = where(mod(coor,B1)==Integer(3)   ,Uorg*f*f,U); 

	// Perp links
       	U_mu = where(mod(coor,B1)==Integer(B1-5),Uorg_mu*f*f,U_mu);
	U_mu = where(mod(coor,B1)==Integer(3)   ,Uorg_mu*f*f,U_mu);

	PokeIndex<LorentzIndex>(U, U_mu, mu);
      }
    }
    for(int mu=0;mu<Nd;mu++) {
      Integer B1 = Block[mu];
      if ( B1 && (B1 <= Global[mu]) ) {

	auto U_mu   = PeekIndex<LorentzIndex>(U,mu);
	auto Uorg_mu= PeekIndex<LorentzIndex>(Uorg,mu);
	// In the plane
	U = where(mod(coor,B1)==Integer(B1-3),Uorg*f*f*f,U); 
	U = where(mod(coor,B1)==Integer(2)   ,Uorg*f*f*f,U); 

	// Perp links
       	U_mu = where(mod(coor,B1)==Integer(B1-4),Uorg_mu*f*f*f,U_mu);
	U_mu = where(mod(coor,B1)==Integer(2)   ,Uorg_mu*f*f*f,U_mu);

	PokeIndex<LorentzIndex>(U, U_mu, mu);
      }
    }
    for(int mu=0;mu<Nd;mu++) {
      Integer B1 = Block[mu];
      if ( B1 && (B1 <= Global[mu]) ) {

	auto U_mu   = PeekIndex<LorentzIndex>(U,mu);
	auto Uorg_mu= PeekIndex<LorentzIndex>(Uorg,mu);
	// In the plane
	U = where(mod(coor,B1)==Integer(B1-2),zzz,U); 
	U = where(mod(coor,B1)==Integer(1)   ,zzz,U); 

	// Perp links
	U_mu = where(mod(coor,B1)==Integer(B1-3),Uorg_mu*f*f*f*f,U_mu);
	U_mu = where(mod(coor,B1)==Integer(1)   ,Uorg_mu*f*f*f*f,U_mu);

	PokeIndex<LorentzIndex>(U, U_mu, mu);
      }
    }
  }
};

NAMESPACE_END(Grid);
