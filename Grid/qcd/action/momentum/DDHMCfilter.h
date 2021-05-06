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

////////////////////////////////////////////////////
// DDHMC filter with sub-block size B[mu]
////////////////////////////////////////////////////

template<typename MomentaField>
struct DDHMCFilter: public MomentumFilterBase<MomentaField>
{
  typedef typename MomentaField::vector_type vector_type; //SIMD-vectorized complex type
  typedef typename MomentaField::scalar_type scalar_type; //scalar complex type
  typedef iVector<iScalar<iScalar<vector_type> >, Nd > LorentzScalarType; //complex phase for each site/direction
  typedef iScalar<iScalar<iScalar<vector_type> > >            ScalarType; //complex phase for each site
  typedef Lattice<LorentzScalarType> LatticeLorentzScalarType;
  typedef Lattice<ScalarType> LatticeScalarType;
  
  Coordinate Block;

  // Could also zero links in plane like Luscher advocates.
  
  DDHMCFilter(const Coordinate &_Block): Block(_Block){}

  void applyFilter(MomentaField &P) const override
  {
    GridBase *grid = P.Grid();
    LatticeScalarType mask_mu(grid);
    LatticeLorentzScalarType mask(grid); 

    ////////////////////////////////////////////////////
    // Zero strictly links crossing between domains
    // Luscher also zeroes links in plane of domain boundaries
    // Keeping interior only. This prevents force from plaquettes
    // crossing domains and keeps whole MD trajectory local.
    // Step further than I want to go.
    ////////////////////////////////////////////////////
    ScalarType zz = scalar_type(0.0);
    ScalarType one= scalar_type(1.0);

    LatticeInteger coor(grid); 
    for(int mu=0;mu<Nd;mu++) {

      LatticeCoordinate(coor,mu);

      mask_mu = where(mod(coor,Block[mu])==Block[mu]-1,zz,one);

      PokeIndex<LorentzIndex>(mask, mask_mu, mu);
      
    }
  }
};

