/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/DDHMC.h

Copyright (C) 2021

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Christopher Kelly

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

NAMESPACE_BEGIN(Grid);
////////////////////////////////////////////////////
// DDHMC filter with sub-block size B[mu]
////////////////////////////////////////////////////

template<typename MomentaField>
struct DDHMCFilter: public MomentumFilterBase<MomentaField>
{
  //  typedef typename MomentaField::vector_type vector_type; //SIMD-vectorized complex type
  //  typedef typename MomentaField::scalar_type scalar_type; //scalar complex type
  //  typedef iVector<iScalar<iScalar<vector_type> >, Nd > LorentzScalarType; //complex phase for each site/direction
  //  typedef iScalar<iScalar<iScalar<vector_type> > >            ScalarType; //complex phase for each site
  //  typedef Lattice<LorentzScalarType> LatticeLorentzScalarType;
  //  typedef Lattice<ScalarType> LatticeScalarType;
  
  Coordinate Block;
  int Width;

  // Could also zero links in plane like Luscher advocates.
  
  DDHMCFilter(const Coordinate &_Block,int _Width): Block(_Block), Width(_Width)
  {
    assert( (Width==0) || (Width==1) || (Width==2) );
  }

  void applyFilter(MomentaField &P) const override
  {
    GridBase *grid = P.Grid();

    LatticeColourMatrix zz(grid); zz = Zero();
    MomentaField zzz(grid); zzz = Zero();

    ////////////////////////////////////////////////////
    // Zero strictly links crossing between domains
    // Luscher also zeroes links in plane of domain boundaries
    // Keeping interior only. This prevents force from plaquettes
    // crossing domains and keeps whole MD trajectory local.
    // This is the case Width=1.
    // Width = 2 should also decouple rectangle force.
    ////////////////////////////////////////////////////
    Coordinate Global=grid->GlobalDimensions();
    LatticeInteger coor(grid); 
    for(int mu=0;mu<Nd;mu++) {

#define PROJECT_DOMAIN      
      if ( (Block[mu] <= Global[mu])&&(Block[mu]>0) ) {
	LatticeCoordinate(coor,mu);
#if 1
	P = where(mod(coor,Block[mu])==Integer(Block[mu]-3),zzz,P); //width 4
	P = where(mod(coor,Block[mu])==Integer(Block[mu]-2),zzz,P); //width 4
	P = where(mod(coor,Block[mu])==Integer(Block[mu]-1),zzz,P); //width 2
	P = where(mod(coor,Block[mu])==Integer(0),zzz,P);           //width 2
	P = where(mod(coor,Block[mu])==Integer(1),zzz,P);           //width 4
	auto P_mu = PeekIndex<LorentzIndex>(P,mu);                  
	P = where(mod(coor,Block[mu])==Integer(2),zzz,P);           //width 6
	PokeIndex<LorentzIndex>(P, P_mu, mu);
#else
	P = where(mod(coor,Block[mu])==Integer(Block[mu]-1),zzz,P); //width 2
	auto P_mu = PeekIndex<LorentzIndex>(P,mu);                  
	P = where(mod(coor,Block[mu])==Integer(0),zzz,P);           //width 6
	PokeIndex<LorentzIndex>(P, P_mu, mu);
#endif	
      }      
    }
#ifdef PROJECT_DOMAIN    
    LatticeInteger domaincb(grid); domaincb=Zero();
    for(int d=0;d<Nd;d++){
      LatticeCoordinate(coor,d);
      domaincb = domaincb + div(coor,Block[d]);
    }
    P = where(mod(domaincb,2)==Integer(1),P,zzz);
#endif    
  }
};
NAMESPACE_END(Grid);

