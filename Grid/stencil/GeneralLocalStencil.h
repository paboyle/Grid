/*************************************************************************************

     Grid physics library, www.github.com/paboyle/Grid 

     Source file: ./lib/GeneralLocalStencil.h

     Copyright (C) 2019

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

     See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#pragma once
NAMESPACE_BEGIN(Grid);

// Share with Cartesian Stencil
struct GeneralStencilEntry { 
  uint64_t _offset;            // 4 bytes 
  uint8_t _permute;            // 1 bytes // Horrible alignment properties
};
// Could pack to 8 + 4 + 4 = 128 bit and use 

class GeneralLocalStencilView {
 public:
  ////////////////////////////////////////
  // Basic Grid and stencil info
  ////////////////////////////////////////
  int                               _npoints; // Move to template param?
  GeneralStencilEntry*  _entries_p;

  accelerator_inline GeneralStencilEntry * GetEntry(int point,int osite) { 
    return & this->_entries_p[point+this->_npoints*osite]; 
  }

};
////////////////////////////////////////
// The Stencil Class itself
////////////////////////////////////////
class GeneralLocalStencil : public GeneralLocalStencilView {
public:
  typedef GeneralLocalStencilView View_type;

protected:
  GridBase *                        _grid;

public: 
  GridBase *Grid(void) const { return _grid; }

  View_type View(void) const {
    View_type accessor(*( (View_type *) this));
    return accessor;
  }

  // Resident in managed memory
  Vector<GeneralStencilEntry>  _entries; 

  GeneralLocalStencil(GridBase *grid, const std::vector<Coordinate> &shifts)
  {
    int npoints = shifts.size();
    int osites  = grid->oSites();
    
    this->_grid    = grid;
    this->_npoints = npoints;
    this->_entries.resize(npoints* osites);
    this->_entries_p = &_entries[0];


    Coordinate Coor;
    Coordinate NbrCoor;
    for(Integer site=0;site<osites;site++){
      for(Integer ii=0;ii<npoints;ii++){
	Integer lex = site*npoints+ii;
	GeneralStencilEntry SE;
	////////////////////////////////////////////////
	// Outer index of neighbour Offset calculation
	////////////////////////////////////////////////
	grid->oCoorFromOindex(Coor,site);
	for(int d=0;d<Coor.size();d++){
	  int rd = grid->_rdimensions[d];
	  NbrCoor[d] = (Coor[d] + shifts[ii][d] + rd )%rd;
	}
	SE._offset      = grid->oIndexReduced(NbrCoor);

	////////////////////////////////////////////////
	// Inner index permute calculation
	// Simpler version using icoor calculation
	////////////////////////////////////////////////
	SE._permute =0;
	for(int d=0;d<Coor.size();d++){

	  int fd = grid->_fdimensions[d];
	  int rd = grid->_rdimensions[d];
	  int ly = grid->_simd_layout[d];

	  assert((ly==1)||(ly==2));

	  int shift = (shifts[ii][d]+fd)%fd;  // make it strictly positive 0.. L-1
	  int x = Coor[d];                // x in [0... rd-1] as an oSite 

	  int permute_dim  = grid->PermuteDim(d);
	  int permute_slice=0;
	  if(permute_dim){    
	    int  num = shift%rd; // Slice within dest osite cell of slice zero
	    int wrap = shift/rd; // Number of osite local volume cells crossed through
                                  // x+num < rd dictates whether we are in same permute state as slice 0
	    if ( x< rd-num ) permute_slice=wrap;
	    else             permute_slice=(wrap+1)%ly;
	  }
	  if ( permute_slice ) {
	    int ptype       =grid->PermuteType(d);
	    uint8_t mask    =grid->Nsimd() >> (ptype + 1);		
	    SE._permute    |= mask;
	  }
	}	
	////////////////////////////////////////////////
	// Store in look up table
	////////////////////////////////////////////////
	this->_entries[lex] = SE;
      }
    }      
  }
  
};

NAMESPACE_END(Grid);

