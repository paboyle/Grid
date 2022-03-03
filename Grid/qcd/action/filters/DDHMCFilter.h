/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/integrators/DirichletFilter.h

Copyright (C) 2015

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
//--------------------------------------------------------------------
#pragma once

NAMESPACE_BEGIN(Grid);
////////////////////////////////////////////////////
// DDHMC filter with sub-block size B[mu]
////////////////////////////////////////////////////

template<typename GaugeField>
struct DDHMCFilter: public MomentumFilterBase<GaugeField>
{
  Coordinate Block;
  int Width;
  
  DDHMCFilter(const Coordinate &_Block,int _Width=2): Block(_Block) { Width=_Width; }

  void applyFilter(GaugeField &U) const override
  {
    GridBase *grid = U.Grid();
    Coordinate Global=grid->GlobalDimensions();
    GaugeField zzz(grid); zzz = Zero();
    LatticeInteger coor(grid); 
    
    auto zzz_mu = PeekIndex<LorentzIndex>(zzz,0);
    ////////////////////////////////////////////////////
    // Zero BDY layers
    ////////////////////////////////////////////////////
    std::cout<<GridLogMessage<<" DDHMC Force Filter Block "<<Block<<" width " <<Width<<std::endl;
    for(int mu=0;mu<Nd;mu++) {

      Integer B1 = Block[mu];
      if ( B1 && (B1 <= Global[mu]) ) {
	LatticeCoordinate(coor,mu);

	////////////////////////////////
	// OmegaBar - zero all links contained in slice B-1,0 and
	// mu links connecting to Omega
	////////////////////////////////
	if ( Width==1) { 
	  U    = where(mod(coor,B1)==Integer(B1-1),zzz,U);
	  U    = where(mod(coor,B1)==Integer(0)   ,zzz,U); 
	  auto U_mu   = PeekIndex<LorentzIndex>(U,mu);
	  U_mu = where(mod(coor,B1)==Integer(B1-2),zzz_mu,U_mu); 
	  PokeIndex<LorentzIndex>(U, U_mu, mu);
	}
	if ( Width==2) { 
	  U    = where(mod(coor,B1)==Integer(B1-2),zzz,U);
	  U    = where(mod(coor,B1)==Integer(B1-1),zzz,U);
	  U    = where(mod(coor,B1)==Integer(0)   ,zzz,U); 
	  U    = where(mod(coor,B1)==Integer(1)   ,zzz,U); 
	  auto U_mu   = PeekIndex<LorentzIndex>(U,mu);
	  U_mu = where(mod(coor,B1)==Integer(B1-3),zzz_mu,U_mu); 
	  PokeIndex<LorentzIndex>(U, U_mu, mu);
	}
	if ( Width==3) { 
	  U    = where(mod(coor,B1)==Integer(B1-3),zzz,U);
	  U    = where(mod(coor,B1)==Integer(B1-2),zzz,U);
	  U    = where(mod(coor,B1)==Integer(B1-1),zzz,U);
	  U    = where(mod(coor,B1)==Integer(0)   ,zzz,U); 
	  U    = where(mod(coor,B1)==Integer(1)   ,zzz,U); 
	  U    = where(mod(coor,B1)==Integer(2)   ,zzz,U); 
	  auto U_mu   = PeekIndex<LorentzIndex>(U,mu);
	  U_mu = where(mod(coor,B1)==Integer(B1-4),zzz_mu,U_mu); 
	  PokeIndex<LorentzIndex>(U, U_mu, mu);
	}
      }

    }
   
  }
};

NAMESPACE_END(Grid);

