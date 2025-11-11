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
struct IcosahedralStencilEntry { 
  uint64_t _offset;            // 8 bytes 
  uint8_t _is_local;           // is this with our lattice array (else in a comms buf)
  uint8_t _adjoint;            // is this with our lattice array (else in a comms buf)
  uint8_t _polarisation;       // which lorentz index on the neighbours patch
  uint8_t _missing_link;       // 
  uint8_t _permute;            // did this wrap (in T-direction)
};


inline int periAdd(int A,int inc,int L) { return (A+inc+L)%L ; }
inline int periAdd(int A,int inc,int L,int & wrap) {
  int r = (A+inc+L)%L;
  if ( r != (A+inc) ) wrap = 1;
  else wrap =0;
  return r;
}

class IcosahedralStencilView {
 public:
  ////////////////////////////////////////
  // Basic Grid and stencil info
  ////////////////////////////////////////
  int                       _npoints; // Move to template param?
  IcosahedralStencilEntry*  _entries_p;

  accelerator_inline IcosahedralStencilEntry * GetEntry(int point,int osite) const { 
    return & this->_entries_p[point+this->_npoints*osite]; 
  }
  void ViewClose(void){};
};

////////////////////////////////////////
// The Stencil Class itself
////////////////////////////////////////
class IcosahedralStencil : public IcosahedralStencilView {
public:
  typedef IcosahedralStencilView View_type;

protected:
  GridBase *                        _grid;
  GridBase *                        _vertexgrid;

public: 
  GridBase *Grid(void) const { return _grid; }

  View_type View(int mode) const {
    View_type accessor(*( (View_type *) this));
    return accessor;
  }
  // NB x+, y+ is ALWAYS present, as are the forward 3 directions for links owned by each site
  //
  // These are VERTEX mesh neigbours being returned, with isPole indicating we need N/S pole according to
  // hemisphere.
  //
  // If needing edge mesh "neigbours" to assemble loops we must find the mapping of a forward link
  // to a corresponding "backward link" on the pole
  deviceVector<IcosahedralStencilEntry>  _entries;

  
  void GetNbrForPlusX(GridBase *grid,Coordinate &Coor,Coordinate &NbrCoor, int &isPole)
  {
    int nd = grid->Nd();
    int L  = grid->LocalDimensions()[0];
    assert (grid->LocalDimensions()[0] == grid->LocalDimensions()[1]);
    int patch = Coor[nd-1];
    int north = patch/HemiPatches;
    int south = 1-north;
    int HemiPatch = patch%HemiPatches;
    int NbrPatch;
    
    NbrCoor = Coor;
    isPole = 0;
    if ( Coor[0]<(L-1) ) { 
      NbrCoor[0]=Coor[0]+1;
      // Nbr is inside THIS patch
      return;
    }
    if ( north ) {
      assert(Coor[0]==(L-1));
      // Simple connect to southern neighbour tile
      NbrCoor[0]=0;
      NbrCoor[nd-1]=periAdd(HemiPatch,+1,HemiPatches) + SouthernHemisphere;
      return;
    }
    ////////////////////////////////////////////////////////////
    // FIXME:
    // Can store the rotation of polarisation here: get xdir instead of ydir and must take adjoint
    ////////////////////////////////////////////////////////////
    if ( south ) {
      assert(Coor[0]==(L-1));
      if ( Coor[1] == 0 ) {
	isPole = 1;
	NbrCoor[0] = (L-1); // Coordinate of the edge graph site holding the edge for other vertex in pole triangle
	NbrCoor[1] = 0;
	NbrCoor[nd-1]=periAdd(HemiPatch,+1,HemiPatches) + SouthernHemisphere;
	return;
      } else {
	NbrCoor[1] = 0;
	NbrCoor[0] = L-Coor[1];
	NbrCoor[nd-1]=periAdd(HemiPatch,+1,HemiPatches) + SouthernHemisphere;
	return;
      }
    }
    assert(0);
  }
  void GetNbrForPlusY(GridBase *grid,Coordinate &Coor,Coordinate &NbrCoor, int &isPole)
  {
    int nd = grid->Nd();
    int L  = grid->LocalDimensions()[0];
    assert (grid->LocalDimensions()[0] == grid->LocalDimensions()[1]);
    int patch = Coor[nd-1];
    int north = patch/HemiPatches;
    int south = 1-north;
    int HemiPatch = patch%HemiPatches;
    int NbrPatch;
    
    NbrCoor = Coor;
    isPole = 0;
    if ( Coor[1]<(L-1) ) { 
      NbrCoor[1]=Coor[1]+1;
      // Nbr is inside THIS patch
      return;
    }
    if ( south ) {
      // Simple connect to northern neighbour tile
      assert(Coor[1]==(L-1));
      NbrCoor[1]=0;
      NbrCoor[nd-1]=HemiPatch + NorthernHemisphere;
      return;
    }
    ////////////////////////////////////////////////////////////
    // FIXME:
    // Could store the rotation of polarisation here: get xdir instead of ydir and must take adjoint
    // But probaby just write "getLinkPropertiesToCloseTriangleA" "getLinkPropertiesToCloseTriangleB"
    // Or write a 'double store' method to move edge graph to a vertex graph with all fermion transports
    ////////////////////////////////////////////////////////////
    if ( north ) {
      assert(Coor[1]==(L-1));
      if ( Coor[0] == 0 ) {
	isPole = 1;
	NbrCoor[1] = (L-1); // Coordinate of the edge graph site holding the edge for other vertex in pole triangle
	NbrCoor[0] = 0;
	NbrCoor[nd-1]=periAdd(HemiPatch,+1,HemiPatches) + NorthernHemisphere;
	return;
      } else {
	NbrCoor[0] = 0;
	NbrCoor[1] = L-Coor[0]; // x=1 --> y=L-1 for y+
	NbrCoor[nd-1]=periAdd(HemiPatch,+1,HemiPatches) + NorthernHemisphere;
	return;
      }
    }
    assert(0);
  }
  // Missing links are at (0,0) on local patch coordinates in -diagonal direction 
  // We are here returning VERTEX grid coordinates.
  void GetNbrForMinusX(GridBase *grid,Coordinate &Coor,Coordinate &NbrCoor)
  {
    int nd = grid->Nd();
    int L  = grid->LocalDimensions()[0];
    assert (grid->LocalDimensions()[0] == grid->LocalDimensions()[1]);
    int patch = Coor[nd-1];
    int north = patch/HemiPatches;
    int south = 1-north;
    int HemiPatch = patch%HemiPatches;
    int NbrPatch;
    
    NbrCoor = Coor;
    if ( Coor[0]>0 ) { 
      NbrCoor[0]=Coor[0]-1;
      return;
    }
    if ( south ) {
      assert(Coor[0]==0);
      // Simple connect to northern neighbour tile
      NbrCoor[0]=L-1;
      NbrCoor[nd-1]=periAdd(HemiPatch,-1,HemiPatches) + NorthernHemisphere;
      return;
    }
    if ( north ) {
	NbrCoor[0] = L-1-Coor[1];
	NbrCoor[1] = L-1;
	NbrCoor[nd-1]=periAdd(HemiPatch,-1,HemiPatches) + NorthernHemisphere;
	return;
    }      
    assert(0);
  }
  
  void GetNbrForMinusY(GridBase *grid,Coordinate &Coor,Coordinate &NbrCoor)
  {
    int nd = grid->Nd();
    int L  = grid->LocalDimensions()[0];
    assert (grid->LocalDimensions()[0] == grid->LocalDimensions()[1]);
    int patch = Coor[nd-1];
    int north = patch/HemiPatches;
    int south = 1-north;
    int HemiPatch = patch%HemiPatches;
    int NbrPatch;
    
    NbrCoor = Coor;
    if ( Coor[1]>0 ) { 
      NbrCoor[1]=Coor[1]-1;
      return;
    }
    if ( north ) {
      assert(Coor[1]==0);
      // Simple connect to northern neighbour tile
      NbrCoor[1]=L-1;
      NbrCoor[nd-1]=HemiPatch + SouthernHemisphere;
      return;
    }
    if ( south ) {
	NbrCoor[1] = L-1-Coor[0];
	NbrCoor[0] = L-1;
	NbrCoor[nd-1]=periAdd(HemiPatch,-1,HemiPatches) + SouthernHemisphere;
	return;
    }      
    assert(0);
  }
  void GetNbrForMinusDiagonal(GridBase *grid,Coordinate &Coor,Coordinate &NbrCoor,int &missingLink)
  {
    int nd = grid->Nd();
    int L  = grid->LocalDimensions()[0];
    assert (grid->LocalDimensions()[0] == grid->LocalDimensions()[1]);
    int patch = Coor[nd-1];
    int north = patch/HemiPatches;
    int south = 1-north;
    int HemiPatch = patch%HemiPatches;
    int NbrPatch;
    
    NbrCoor = Coor;
    
    missingLink = 0;
    if( Coor[0]==0 && Coor[1]==0) {
      missingLink=1;
      return;
    }
    if ( (Coor[0]>0) && (Coor[1]>0) ) { 
      // Nbr is inside THIS patch
      NbrCoor[0]=Coor[0]-1;
      NbrCoor[1]=Coor[1]-1;
      return;
    }
    if ( north ) {
      if ( Coor[0]==0 ) {
	// We are on -x edge
	// Maps to +y edge of hemipatch to the left
	NbrCoor[nd-1] = periAdd(HemiPatch,-1,HemiPatches) + NorthernHemisphere;
	NbrCoor[0]=(L-Coor[1]);
	NbrCoor[1]=(L-1);
	return;
      } else {
        // We are on the -y edge and NOT bottom left corner; Nbr is in the patch LEFT
	assert( (Coor[0]>0) && (Coor[1]==0) );
	NbrCoor[nd-1] = HemiPatch + SouthernHemisphere; // Map from north to south
	NbrCoor[0]=Coor[0]-1;
	NbrCoor[1]=(L-1);
	return;
      }
      assert(0);
    }
    if ( south ) {
      // We are on the -y edge
      if ( Coor[1]==0 ) {
	NbrCoor[nd-1] = periAdd(HemiPatch,-1,HemiPatches) + SouthernHemisphere;
        NbrCoor[0]=(L-1);
        NbrCoor[1]=(L-Coor[0]);
	return;
      } else {
      // We are on the -x edge and NOT bottom left corner; Nbr is in the patch LEFT
	assert( (Coor[0]==0) && (Coor[1]>0) );
	NbrCoor[nd-1] = periAdd(HemiPatch,-1,HemiPatches) + NorthernHemisphere; // south to north
	NbrCoor[0]=(L-1);
	NbrCoor[1]= Coor[1]-1;
	return;
      }
      assert(0);
    }
    assert(0);
  }
  void GetNbrForPlusDiagonal(GridBase *grid,Coordinate &Coor,Coordinate &NbrCoor)
  {
    int nd = grid->Nd();
    int L  = grid->LocalDimensions()[0];
    assert (grid->LocalDimensions()[0] == grid->LocalDimensions()[1]);
    int patch = Coor[nd-1];
    int north = patch/HemiPatches;
    int south = 1-north;
    int HemiPatch = patch%HemiPatches;
    int NbrPatch;
    
    NbrCoor = Coor;
    
    if ( (Coor[0]<L-1) && (Coor[1]<L-1) ) { 
      // Nbr is inside THIS patch
      NbrCoor[0]=Coor[0]+1; 
      NbrCoor[1]=Coor[1]+1;
      return;
    }

    if ( north ) {
      // We are on +y edge
      if ( Coor[1]==(L-1) ) {
	NbrCoor[nd-1] = periAdd(HemiPatch,+1,HemiPatches) + NorthernHemisphere;
	NbrCoor[0]=0;
	NbrCoor[1]=(L-1)-Coor[0];
	return;
      } else {
      // Else we are on the +x edge, not top right corner
	assert( Coor[0]==(L-1) && Coor[1]<(L-1) );
	NbrCoor[nd-1] = periAdd(HemiPatch,+1,HemiPatches) + SouthernHemisphere;
	NbrCoor[0]=0;
	NbrCoor[1]=Coor[1]+1;
	return;
      }
      assert(0);
    }

    if ( south ) {
      // We are on the +x edge
      if ( Coor[0]==(L-1) ) {
	NbrCoor[nd-1] = periAdd(HemiPatch,+1,HemiPatches) + SouthernHemisphere;
	NbrCoor[0]=(L-1)-Coor[1]; //y=(L-1) <-> x=0 ; y=0<->x=(L-1)  [ Sanity check ]
	NbrCoor[1]=0;
	return;
      } else { 
      // We are on the +y edge and NOT top right corner; Nbr is in the patch UP
	assert( (Coor[1]==L-1) && (Coor[0]<(L-1)) );
	NbrCoor[nd-1] = HemiPatch + NorthernHemisphere;
	NbrCoor[0]=Coor[0]+1;
	NbrCoor[1]=0;
	return;
      }
      assert(0);
    }
    assert(0);
  }
  
  void  TestGeometry(void)
  {
    GridBase *grid = this->_grid;
    uint64_t cart_sites = grid->CartesianOsites();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //    Loop over cart sites.
    //    Find two triangles per site.
    //    Check going forward in X, Up and forward in Diag match
    //    Check going Up, forward in X and forward Diag match; subtleties at poles and rotation in cross patch
    //////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << GridLogMessage<< "*************************************"<<std::endl;
    std::cout << GridLogMessage<< " Icosahedral Stencil Geometry Test !"<<std::endl;
    std::cout << GridLogMessage<< "*************************************"<<std::endl;

    const int triangle_ref = cart_sites;
    std::cout << GridLogMessage<< " Base triangle count for each type " <<triangle_ref;
    
    std::cout << GridLogMessage<< "------------------------------------"<<std::endl;
    std::cout << GridLogMessage<< "testing +x+y vs +diag"<<std::endl;
    std::cout << GridLogMessage<< "testing +y+x vs +diag"<<std::endl;
    std::cout << GridLogMessage<< "------------------------------------"<<std::endl;
    int xyd_pole_count=0;
    int xyd_count=0;
    int yxd_pole_count=0;
    int yxd_count=0;
    for(uint64_t site=0;site<cart_sites; site ++) {

      Coordinate Coor;
      Coordinate DiagCoor;

      int nd = grid->Nd();
      int L  = grid->LocalDimensions()[0];

      grid->oCoorFromOindex(Coor,site);

      int patch = Coor[nd-1];
      int north = patch/HemiPatches;
      int south = 1-north;
      int isPole      = 0;
      int discard;
      int missingLink = 0;
      
      int HemiPatch = patch%HemiPatches;

      //////////////////////////////
      // First test of triangle
      //////////////////////////////
      // Compare +x, +y to +diag
      Coordinate XpCoor;
      Coordinate YpXpCoor;
      GetNbrForPlusDiagonal(grid,Coor,DiagCoor);
      GetNbrForPlusX(grid,Coor,XpCoor,isPole);

      int XpHemiPatch  = XpCoor[nd-1]%HemiPatches;
      int XpHemisphere = XpCoor[nd-1]/HemiPatches;
      
      if (isPole) {
	YpXpCoor = XpCoor;
      } else if ( XpHemiPatch != HemiPatch && south ) {
	GetNbrForMinusX(grid,XpCoor,YpXpCoor);
      } else {
	GetNbrForPlusY(grid,XpCoor,YpXpCoor,discard);
      }
      
      if(isPole) {
	std::cout << GridLogDebug<<"Forward xyd triangle "<<Coor<<"-Pole["<<XpCoor[2]<<"]-"<<YpXpCoor<<" should be " <<DiagCoor<<std::endl;
	xyd_pole_count++;
      } else {
	std::cout << GridLogDebug<<"Forward xyd triangle "<<Coor<<"-"<<XpCoor<<"-"<<YpXpCoor<<" should be " <<DiagCoor<<std::endl;
	xyd_count++;
      }
      for(int d=0;d<DiagCoor.size();d++) {
	assert(DiagCoor[d]==YpXpCoor[d]);
      }

      Coordinate YpCoor;
      Coordinate XpYpCoor;
      GetNbrForPlusDiagonal(grid,Coor,DiagCoor);
      GetNbrForPlusY(grid,Coor,YpCoor,isPole);

      int YpHemiPatch  = YpCoor[nd-1]%HemiPatches;
      int YpHemisphere = YpCoor[nd-1]/HemiPatches;
      
      if(isPole) {
	XpYpCoor = YpCoor;
      } else if ( YpHemiPatch != HemiPatch  && north ) {
	GetNbrForMinusY(grid,YpCoor,XpYpCoor); // we hopped - this rotates the directions
      } else {
	GetNbrForPlusX(grid,YpCoor,XpYpCoor,discard);
      }

      if(isPole) {
	yxd_pole_count++;
	std::cout << GridLogDebug<<"Forward yxd triangle "<<Coor<<"-Pole["<<YpCoor[2]<<"]-"<<XpYpCoor<<" should be " <<DiagCoor<<std::endl;
      } else {
	yxd_count++;
	std::cout <<GridLogDebug << "Forward yxd triangle "<<Coor<<"-"<<YpCoor<<"-"<<XpYpCoor<<" should be " <<DiagCoor<<std::endl;
      }
      for(int d=0;d<DiagCoor.size();d++) {
	assert(DiagCoor[d]==XpYpCoor[d]);
      }
    }
    std::cout << GridLogMessage<< " xyd_count "<<xyd_count<<" + poles_count "<<xyd_pole_count<<" expect "<<triangle_ref<<" triangles "<<std::endl;
    std::cout << GridLogMessage<<" yxd_count "<<yxd_count<<" + poles_count "<<yxd_pole_count<<" expect "<<triangle_ref<<" triangles "<<std::endl;
    assert(xyd_count+xyd_pole_count == triangle_ref);
    assert(yxd_count+yxd_pole_count == triangle_ref);
    std::cout << GridLogMessage<< "------------------------------------"<<std::endl;
    std::cout << GridLogMessage<< "testing -diag +x+y = identity"<<std::endl;
    std::cout << GridLogMessage<< "testing -diag +y+x = identity"<<std::endl;
    std::cout << GridLogMessage<< "------------------------------------"<<std::endl;

    int dmxy_count=0;
    int dmyx_count=0;
    int dmxy_count_special=0;
    int dmyx_count_special=0;
    int num_missing=0;

    for(uint64_t site=0;site<cart_sites; site ++) {

      Coordinate Coor;

      int nd = grid->Nd();
      int L  = grid->LocalDimensions()[0];

      grid->oCoorFromOindex(Coor,site);

      int patch = Coor[nd-1];
      int north = patch/HemiPatches;
      int south = 1-north;
      int isPole      = 0;
      int discard;
      int missingLink = 0;
      int HemiPatch = patch%HemiPatches;

      Coordinate DmCoor;
      GetNbrForMinusDiagonal(grid,Coor,DmCoor,missingLink);
      if ( missingLink ) {
	std::cout << GridLogDebug<< Coor << " has no backwards diagonal link "<<std::endl;
	num_missing++;
      } else {

	int DmPatch  = DmCoor[nd-1];
	int DmHemiPatch  = DmCoor[nd-1]%HemiPatches;
	int DmHemisphere = DmCoor[nd-1]/HemiPatches;

	Coordinate XpDmCoor;
	Coordinate YpXpDmCoor;

	Coordinate YpDmCoor;
	Coordinate XpYpDmCoor;
	
	if ( DmHemiPatch != HemiPatch && north ) {
	  GetNbrForPlusDiagonal(grid,DmCoor,XpDmCoor);
	  GetNbrForPlusY(grid,XpDmCoor,YpXpDmCoor,isPole); assert(!isPole);

	  GetNbrForMinusX(grid,DmCoor,YpDmCoor);
	  GetNbrForPlusDiagonal(grid,YpDmCoor,XpYpDmCoor);
	  dmxy_count_special++;
	  dmyx_count_special++;
	} else if ( DmPatch == periAdd(HemiPatch,-1,HemiPatches) && south ) {
	  GetNbrForPlusDiagonal(grid,DmCoor,YpDmCoor);
	  GetNbrForPlusX(grid,YpDmCoor,XpYpDmCoor,isPole); assert(!isPole);

	  GetNbrForMinusY(grid,DmCoor,XpDmCoor);
	  GetNbrForPlusDiagonal(grid,XpDmCoor,YpXpDmCoor);
	  dmxy_count_special++;
	  dmyx_count_special++;
	} else {
	  GetNbrForPlusX(grid,DmCoor,XpDmCoor,isPole); assert(!isPole);
	  GetNbrForPlusY(grid,XpDmCoor,YpXpDmCoor,isPole); assert(!isPole);

	  GetNbrForPlusY(grid,DmCoor,YpDmCoor,isPole); assert(!isPole);
	  GetNbrForPlusX(grid,YpDmCoor,XpYpDmCoor,isPole); assert(!isPole);
	  dmxy_count++;
	  dmyx_count++;
	}
	std::cout<< GridLogDebug << Coor<<" DmXpYp triangle YpXpDm"<<YpXpDmCoor<<"-XpDm"<<XpDmCoor<<"-Dm"<<DmCoor<<" should be " <<Coor<<std::endl;
	for(int d=0;d<Coor.size();d++) {
	  assert(Coor[d]==YpXpDmCoor[d]);
	}

	std::cout << GridLogDebug<< Coor<<"DmXpYp triangle XpYpDm"<<XpYpDmCoor<<"-YpDm"<<YpDmCoor<<"-Dm"<<DmCoor<<" should be " <<Coor<<std::endl;
	for(int d=0;d<Coor.size();d++) {
	  assert(Coor[d]==XpYpDmCoor[d]);
	}
      }
    }
    std::cout <<GridLogMessage<<" dmxy_count "<<dmxy_count<<" + special "<<dmxy_count_special<<" + missing "<<num_missing<<" expect "<<triangle_ref<<" triangles "<<std::endl;
    std::cout <<GridLogMessage<<" dmyx_count "<<dmyx_count<<" + special "<<dmyx_count_special<<" + missing "<<num_missing<<" expect "<<triangle_ref<<" triangles "<<std::endl;
    assert(dmxy_count + dmxy_count_special + num_missing == triangle_ref);
    assert(dmyx_count + dmyx_count_special + num_missing == triangle_ref);
    
    std::cout << GridLogMessage<< "*************************************"<<std::endl;
    std::cout << GridLogMessage<< " Icosahedral Stencil Geometry Test Complete"<<std::endl;
    std::cout << GridLogMessage<< "*************************************"<<std::endl;
  }
  IcosahedralStencil(GridBase *grid,GridBase *vertexgrid) 
  {
    this->_grid          = grid;
    this->_vertexgrid    = vertexgrid;
    assert(grid->ProcessorCount() ==1);
    // Loop over L^2 x T x npatch and the
    assert(grid->isIcosahedral());
    assert(grid->isIcosahedral());
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // VertexInputs = true  implies the neighbour has vertex support
  // VertexInputs = false implies the neighbout has edge support
  // 
  // isVertex implies must generate stencil entries to evaluate result on north/south pole
  //
  // These are independent:
  //     can apply a vertex support gauge transform to edge supported gauge field
  //     can apply a vertex supported link double store to edge supported gauge field
  //     can apply a vertex supported laplace or dirac operator vertex supported matter field
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void NearestNeighbourStencil(int vertexInputs,int vertexOutputs)
  {
    GridBase * grid = this->_grid; // the edge grid
    GridBase * vertexgrid = this->_vertexgrid;
    
    uint64_t cart_sites = grid->CartesianOsites();
    uint64_t Npole_sites = vertexgrid->NorthPoleOsites();
    uint64_t Spole_sites = vertexgrid->SouthPoleOsites();

    Coordinate pcoor = grid->ThisProcessorCoor();
    Coordinate pgrid = grid->ProcessorGrid();
    /*
     * resize the stencil entries array and set npoints
     */
    const int np=8;
    this->_npoints=np; // Move to template param?
    if ( vertexOutputs ) {
      this->_entries.resize(this->_npoints * (cart_sites+Npole_sites+Spole_sites));
    } else {
      this->_entries.resize(this->_npoints * cart_sites);
    }
    this->_entries_p = &_entries[0];

    int nd = grid->Nd();
    int L  = grid->LocalDimensions()[0];

    for(uint64_t site=0;site<cart_sites; site ++) {

      Coordinate Coor;
      Coordinate NbrCoor;

      Integer lexXp = site*np  ;
      Integer lexYp = site*np+1;
      Integer lexDp = site*np+2;
      Integer lexXm = site*np+3;
      Integer lexYm = site*np+4;
      Integer lexDm = site*np+5;
      Integer lexTp = site*np+6;
      Integer lexTm = site*np+7;

      IcosahedralStencilEntry SE;
	
      ////////////////////////////////////////////////
      // Outer index of neighbour Offset calculation
      ////////////////////////////////////////////////
      grid->oCoorFromOindex(Coor,site);
      NbrCoor = Coor;
      assert( grid->LocalDimensions()[1]==grid->LocalDimensions()[0]);
      assert( grid->_simd_layout[0]==1); // Cannot vectorise in these dims
      assert( grid->_simd_layout[1]==1);
      assert( grid->_processors[0]==1); // Cannot mpi distribute in these dims
      assert( grid->_processors[1]==1);

      int Patch = Coor[nd-1];
      int HemiPatch = Patch%HemiPatches;
      int Hemisphere= Patch/HemiPatches;
      int north = Patch/HemiPatches;
      int south = 1-north;
      int isPoleY;
      int isPoleX;
      int missingLink;
      assert(Patch<IcosahedralPatches);
      assert((north==1)||(south==1));

      /*
       * Just get all six neighbours (if present).
       */
      Coordinate XpCoor;
      Coordinate YpCoor;
      Coordinate DpCoor;
      Coordinate XmCoor;
      Coordinate YmCoor;
      Coordinate DmCoor;
      Coordinate TpCoor;
      Coordinate TmCoor;

      GetNbrForPlusDiagonal(grid,Coor,DpCoor);
      GetNbrForPlusX(grid,Coor,XpCoor,isPoleX);
      GetNbrForPlusY(grid,Coor,YpCoor,isPoleY);

      GetNbrForMinusDiagonal(grid,Coor,DmCoor,missingLink);
      GetNbrForMinusX(grid,Coor,XmCoor);
      GetNbrForMinusY(grid,Coor,YmCoor);

      int tdim  = 2; int delta = 1;
      int permuteTp;
      int permuteTm;
      GetOrthogNbr(grid,Coor,TpCoor,tdim, delta,permuteTp);
      GetOrthogNbr(grid,Coor,TmCoor,tdim,-delta,permuteTm);
	
      int DpPatch      = DpCoor[nd-1];
      int DpHemiPatch  = DpCoor[nd-1]%HemiPatches;
      int DpHemisphere = DpCoor[nd-1]/HemiPatches;

      int YpHemiPatch  = YpCoor[nd-1]%HemiPatches;
      int XpHemiPatch  = XpCoor[nd-1]%HemiPatches;
      // For negative direction cannot use the Diagonal link
      // as this may not be present on the 5-points
      // Makes for a hemisphere dependent behaviour
      int XmHemiPatch  = XmCoor[nd-1]%HemiPatches;
      int XmHemisphere = XmCoor[nd-1]/HemiPatches;
      int YmHemiPatch  = YmCoor[nd-1]%HemiPatches;
      int YmHemisphere = YmCoor[nd-1]/HemiPatches;

      int DmPatch      = DmCoor[nd-1];
      int DmHemiPatch  = DmCoor[nd-1]%HemiPatches;
      int DmHemisphere = DmCoor[nd-1]/HemiPatches;

      SE._permute=0;
      if ( vertexInputs ) {// Neighbour will live on poles and peer point
	////////////////////////////////////////////////
	// XpCoor stencil entry; consider isPole case
	////////////////////////////////////////////////
	SE._is_local     = true;
	SE._missing_link = false;
	SE._offset       = grid->oIndex(XpCoor);
	if ( isPoleX ) {
	  SE._offset     = vertexgrid->PoleSiteForOcoor(Coor);
	  //	  std::cout << site<<" setting X-Pole site "<<SE._offset<<" for coor "<<Coor<<std::endl;
	}
	SE._polarisation = IcosahedronPatchY;
	SE._adjoint      = false;
	acceleratorPut(this->_entries[lexXp],SE);
	////////////////////////////////////////////////
	// for YpCoor
	////////////////////////////////////////////////
	SE._is_local     = true;
	SE._missing_link = false;
	SE._offset       = grid->oIndex(YpCoor);
	if ( isPoleY ) {
	  SE._offset     = vertexgrid->PoleSiteForOcoor(Coor);
	  //	  std::cout << site<<" setting Y-Pole site "<<SE._offset<<" for coor "<<Coor<<std::endl;
	}
	SE._polarisation = IcosahedronPatchX;
	SE._adjoint      = false;
	acceleratorPut(this->_entries[lexYp],SE);
      } else { // Neighbour will be a forward edge and connection may be more complicated
	////////////////////////////////////////////////
	// XpCoor stencil entry
	// Store in look up table
	////////////////////////////////////////////////
	// Basis rotates dictates BOTH adjoint and polarisation
	// Could reduce the amount of information stored here
	SE._is_local     = true;
	SE._missing_link = false;
	SE._offset       = grid->oIndex(XpCoor);
	SE._polarisation = IcosahedronPatchY;
	SE._adjoint      = false;
	if ( DpHemiPatch != HemiPatch && south ) { // These are the sneaky redirect for edge / faces
	  SE._offset       = grid->oIndex(DpCoor);
	  SE._polarisation = IcosahedronPatchX;
	  SE._adjoint      = true;
	}
	acceleratorPut(this->_entries[lexXp],SE);
	////////////////////////////////////////////////
	// for YpCoor
	////////////////////////////////////////////////
	SE._is_local     = true;
	SE._missing_link = false;
	SE._offset       = grid->oIndex(YpCoor);
	SE._polarisation = IcosahedronPatchX;
	SE._adjoint      = false;
	if ( YpHemiPatch != HemiPatch && north ) {  // These are the sneaky redirect for edge / faces
	  SE._offset       = grid->oIndex(DpCoor);
	  SE._polarisation = IcosahedronPatchY;
	  SE._adjoint      = true;
	}
	acceleratorPut(this->_entries[lexYp],SE);
      }

      SE._adjoint      = false;
      SE._is_local     = true;
      SE._missing_link = false;
      ////////////////////////////////////////////////
      // XmCoor stencil entry
      // Store in look up table
      ////////////////////////////////////////////////
      SE._offset       = grid->oIndex(XmCoor);
      SE._polarisation = IcosahedronPatchDiagonal;
      if ( XmHemiPatch != HemiPatch && north ) {
	SE._polarisation = IcosahedronPatchY; // nbrs Y instead of diagonal in North hemisphere exceptional case
      }
      acceleratorPut(this->_entries[lexXm],SE);
      
      ////////////////////////////////////////////////
      // for YmCoor
      ////////////////////////////////////////////////
      SE._offset       = grid->oIndex(YmCoor);
      SE._polarisation = IcosahedronPatchDiagonal;
      if ( YmHemiPatch != HemiPatch && south ) {
	SE._polarisation = IcosahedronPatchX;   // Basis rotates
      }
      acceleratorPut(this->_entries[lexYm],SE);
      
      /////////////////////////////////////////////////////////////////////
      // for DpCoor ; never needed for staples, only for vertex diff ops
      // no polarisation rotation
      /////////////////////////////////////////////////////////////////////
      SE._offset       = grid->oIndex(DpCoor);
      SE._polarisation = IcosahedronPatchDiagonal; // should ignore
      acceleratorPut(this->_entries[lexDp],SE);
    
      /////////////////////////////////////////////////////////////////////
      // for DmCoor ; never needed for staples, only for vertex diff ops
      // No polarisation rotation.
      // But polarisation rotation is needed for double storing.
      /////////////////////////////////////////////////////////////////////
      SE._offset       = grid->oIndex(DmCoor);
      SE._polarisation = IcosahedronPatchDiagonal; // default
      if ( (DmHemiPatch != HemiPatch) && (DmHemisphere==Hemisphere)  && south ) {
        SE._polarisation = IcosahedronPatchX;   // Basis rotates
      }
      if ( DmHemiPatch != HemiPatch && (DmHemisphere==Hemisphere) && north ) {
        SE._polarisation = IcosahedronPatchY;   // Basis rotates
      }
      SE._missing_link = missingLink;
      acceleratorPut(this->_entries[lexDm],SE);
      /////////////////////////////////////////////////////////////////////
      // for Tp/mCoor 
      /////////////////////////////////////////////////////////////////////
      SE._polarisation = 0;
      SE._offset       = grid->oIndex(TpCoor);
      SE._permute      = permuteTp;
      acceleratorPut(this->_entries[lexTp],SE);

      SE._offset       = grid->oIndex(TmCoor);
      SE._permute      = permuteTm;
      acceleratorPut(this->_entries[lexTm],SE);
      SE._permute      = 0;
    }

    if ( vertexOutputs ) {
      int ndm1 = grid->Nd()-1;
      if ( vertexgrid->ownsSouthPole() ) {
	IcosahedralStencilEntry SE;
	for(uint64_t site=0;site<cart_sites; site ++) { // loops over volume
	  Coordinate Coor;
	  Coordinate tCoor;
	  vertexgrid->oCoorFromOindex(Coor,site);
	  int north = Coor[ndm1]/HemiPatches;
	  int south = 1-north;
	  if( (Coor[0]==(L-1))&&(Coor[1]==0) &&south ) {
	    int64_t pole_site = vertexgrid->PoleSiteForOcoor(Coor);
	    int64_t lex       = pole_site*np+(Coor[ndm1]%HemiPatches);

	    SE._offset       = site;
	    SE._is_local     = true;
	    SE._polarisation = IcosahedronPatchX;  // ignored
	    SE._adjoint      = false;              // ignored
	    SE._missing_link = false;
	    SE._permute=0;
	    acceleratorPut(this->_entries[lex],SE);
	    
	    int64_t lex5      = pole_site*np+5; // We miss the backwards link
	    SE._missing_link = true;
	    acceleratorPut(this->_entries[lex5],SE);

	    int tdim = 2;
	    int Rt = vertexgrid->_rdimensions[tdim];
	    int permute;
	    int64_t nbr_pole_site;

	    lex = pole_site*np+6;// tp
	    tCoor = Coor;
	    tCoor[tdim] = periAdd(tCoor[tdim],1,Rt,permute);
	    nbr_pole_site = vertexgrid->PoleSiteForOcoor(tCoor);
	    SE._offset = nbr_pole_site;
	    SE._permute= permute;
	    acceleratorPut(this->_entries[lex],SE);

	    lex = pole_site*np+7;// tm
	    tCoor = Coor;
	    tCoor[tdim] = periAdd(tCoor[tdim],-1,Rt,permute);
	    //	    std::cout << " pole_site "<<pole_site<<" t"<<Coor[tdim]<<" tm "<<tCoor[tdim]<<" perm "<<permute<<std::endl;
	    nbr_pole_site = vertexgrid->PoleSiteForOcoor(tCoor);
	    SE._offset = nbr_pole_site;
	    SE._permute= permute;
	    acceleratorPut(this->_entries[lex],SE);

	  }
	}
      }
      if ( vertexgrid->ownsNorthPole() ) {
	IcosahedralStencilEntry SE;
	for(uint64_t site=0;site<cart_sites; site ++) {
	  Coordinate Coor;
	  Coordinate tCoor;
	  vertexgrid->oCoorFromOindex(Coor,site);
	  int north = Coor[ndm1]/HemiPatches;
	  if( (Coor[0]==0)&&(Coor[1]==(L-1))&&north ) {
	    int64_t pole_site = vertexgrid->PoleSiteForOcoor(Coor);
	    int64_t lex       = pole_site*np+(Coor[ndm1]%HemiPatches);
	    //	    std::cout << "Coor "<<Coor<<" connects to north pole_site "<<pole_site<<std::endl;
	    SE._offset       = site;
	    SE._is_local     = true;
	    SE._polarisation = IcosahedronPatchY;  // ignored
	    SE._adjoint      = false;              // ignored
	    SE._missing_link = false;
	    SE._permute=0;
	    acceleratorPut(this->_entries[lex],SE);
	    
	    int64_t lex5     = pole_site*np+5; // We miss the backwards link
	    SE._missing_link = true;
	    acceleratorPut(this->_entries[lex5],SE);

	    int tdim = 2;
	    int Rt = vertexgrid->_rdimensions[tdim];
	    int permute;
	    int64_t nbr_pole_site;

	    
	    lex = pole_site*np+6;// tp
	    tCoor = Coor;
	    tCoor[tdim] = periAdd(tCoor[tdim],1,Rt,permute);
	    nbr_pole_site = vertexgrid->PoleSiteForOcoor(tCoor);
	    SE._offset = nbr_pole_site;
	    SE._permute= permute;
	    acceleratorPut(this->_entries[lex],SE);
	    //	    std::cout << " Put nbr "<<SE._offset<<" for north site "<<lex<<std::endl;
	    
	    lex = pole_site*np+7;// tm
	    tCoor = Coor;
	    tCoor[tdim] = periAdd(tCoor[tdim],-1,Rt,permute);
	    nbr_pole_site = vertexgrid->PoleSiteForOcoor(tCoor);
	    SE._offset = nbr_pole_site;
	    SE._permute= permute;
	    acceleratorPut(this->_entries[lex],SE);
	    //	    std::cout << " Put nbr "<<SE._offset<<" for north site "<<lex<<std::endl;
	    
	  }
	}
      }
    }
  }
    /*************************************************************
     * For gauge action implementation
     *************************************************************
     */
  void FaceStencil(void)
  {
    GridBase * grid = this->_grid;

    int osites  = grid->oSites();

    uint64_t cart_sites = grid->CartesianOsites();
    uint64_t Npole_sites = grid->NorthPoleOsites();
    uint64_t Spole_sites = grid->SouthPoleOsites();
    Coordinate pcoor = grid->ThisProcessorCoor();
    Coordinate pgrid = grid->ProcessorGrid();

    /*
     * resize the stencil entries array and set npoints
     */
    this->_npoints=2; // Move to template param?
    this->_entries.resize(this->_npoints * cart_sites);
    this->_entries_p = &_entries[0];

    for(uint64_t site=0;site<cart_sites; site ++) {

      Coordinate Coor;

      int nd = grid->Nd();
      int L  = grid->LocalDimensions()[0];

      Integer lexXY = site*2;
      Integer lexYX = site*2+1;

      IcosahedralStencilEntry SE;
	
      ////////////////////////////////////////////////
      // Outer index of neighbour Offset calculation
      ////////////////////////////////////////////////
      grid->oCoorFromOindex(Coor,site);
      assert( grid->LocalDimensions()[1]==grid->LocalDimensions()[0]);
      assert( grid->_simd_layout[0]==1); // Cannot vectorise in these dims
      assert( grid->_simd_layout[1]==1);
      assert( grid->_processors[0]==1); // Cannot mpi distribute in these dims
      assert( grid->_processors[1]==1);

      int Patch = Coor[nd-1];
      int HemiPatch = Patch%HemiPatches;
      int north = Patch/HemiPatches;
      int south = 1-north;
      int isPoleY;
      int isPoleX;
      
      assert(Patch<IcosahedralPatches);
      assert((north==1)||(south==1));
	/*
	 * Minimal stencil for edge -> face triangle evaluation
	 *
	 * On each edge grid site, hold +x,+y,+diag links
	 * Must locate the closing link to form the two forward triangles
	 *
	 * case: +x,+d 
	 *
	 *   north: +x neighbours ydir link always
	 *   south: +x neighbours ydir link except
	 *            cross into a different patch in south
	 *         then
	 *          +d neighbours xdir link
	 *
	 * case: +y,+d 
	 *   south: +y neighbours xdir link always
	 *   north: +y neighbours xdir link unless
	 *            cross into a different patch in north
	 *         then
	 *          +d neighbours ydir link
	 *
	 */
      Coordinate XpCoor;
      Coordinate YpCoor;
      Coordinate DpCoor;

      GetNbrForPlusDiagonal(grid,Coor,DpCoor);
      GetNbrForPlusX(grid,Coor,XpCoor,isPoleX);
      GetNbrForPlusY(grid,Coor,YpCoor,isPoleY);

      //      int XpHemiPatch  = XpCoor[nd-1]%HemiPatches;
      //      int XpHemisphere = XpCoor[nd-1]/HemiPatches;

      int DpPatch      = DpCoor[nd-1];
      int DpHemiPatch  = DpCoor[nd-1]%HemiPatches;
      int DpHemisphere = DpCoor[nd-1]/HemiPatches;

      ////////////////////////////////////////////////
      // for trace [ U_x(z) U_y(z+\hat x) adj(U_d(z)) ]
      ////////////////////////////////////////////////
      if ( DpHemiPatch != HemiPatch && south ) {
	SE._offset       = grid->oIndex(DpCoor);
	SE._is_local     = true;
	SE._polarisation = IcosahedronPatchX;
	SE._adjoint      = true;
	SE._missing_link = false;
	SE._permute=0;
      } else {
	SE._offset       = grid->oIndex(XpCoor);
	SE._is_local     = true;
	SE._polarisation = IcosahedronPatchY;
	SE._adjoint      = false;
	SE._missing_link = false;
	SE._permute=0;
      }

      ////////////////////////////////////////////////
      // Store in look up table
      ////////////////////////////////////////////////
      acceleratorPut(this->_entries[lexXY],SE);
      
      // failed in the if case here
      ////////////////////////////////////////////////
      // for trace [  U_y(z) U_x(z+\hat y) adj(U_d(z))   ]
      ////////////////////////////////////////////////
      int YpHemiPatch  = YpCoor[nd-1]%HemiPatches;
      if ( YpHemiPatch != HemiPatch && north ) {
	SE._offset       = grid->oIndex(DpCoor);
	SE._is_local     = true;
	SE._polarisation = IcosahedronPatchY;
	SE._adjoint      = true;
	SE._missing_link = false;
	SE._permute=0;
      } else {
	SE._offset       = grid->oIndex(YpCoor);
	SE._is_local     = true;
	SE._polarisation = IcosahedronPatchX;
	SE._adjoint      = false;
	SE._missing_link = false;
	SE._permute=0;
      }
      ////////////////////////////////////////////////
      // Store in look up table
      ////////////////////////////////////////////////
      acceleratorPut(this->_entries[lexYX],SE);
    };
  }
  //
  // Orthogonal direction support
  //
  // Plaquettes:
  // Must be able to get the sites +T, and +X, +Y, +D
  //
  // Staples
  // Must be able to get the sites (+T, and +X, +Y, +D) (-T, and -T+X, -T+Y, -T+D)
  //
  // Laplacian:
  // Must be able to get the sites +-T
  //
  void GetOrthogNbr(GridBase *grid,Coordinate &Coor,Coordinate &NbrCoor,int dim,int delta, int & permute)
  {
    assert(delta==1 || delta==-1);
    int L = grid->_rdimensions[dim];

    NbrCoor = Coor;
    NbrCoor[dim] = periAdd(Coor[dim],delta,L,permute);
  }
  void TemporalPlaquetteStencil(void)
  {
    GridBase * grid = this->_grid; // the edge grid
    GridBase * vertexgrid = this->_vertexgrid;

    int osites  = grid->oSites();

    uint64_t cart_sites = grid->CartesianOsites();
    uint64_t Npole_sites = grid->NorthPoleOsites();
    uint64_t Spole_sites = grid->SouthPoleOsites();
    Coordinate pcoor = grid->ThisProcessorCoor();
    Coordinate pgrid = grid->ProcessorGrid();

    /*
     * resize the stencil entries array and set npoints
     */
    this->_npoints=4; 
    this->_entries.resize(this->_npoints * cart_sites);
    this->_entries_p = &_entries[0];

    for(uint64_t site=0;site<cart_sites; site ++) {

      int nd = grid->Nd();
      int L  = grid->LocalDimensions()[0];

      Coordinate Coor;
      IcosahedralStencilEntry SE;

      Integer lexT = site*4+0;
      Integer lexX = site*4+1;
      Integer lexY = site*4+2;
      Integer lexD = site*4+3;

      ////////////////////////////////////////////////
      // Outer index of neighbour Offset calculation
      ////////////////////////////////////////////////
      grid->oCoorFromOindex(Coor,site);
      assert( grid->LocalDimensions()[1]==grid->LocalDimensions()[0]);
      assert( grid->_simd_layout[0]==1); // Cannot vectorise in these dims
      assert( grid->_simd_layout[1]==1);
      assert( grid->_processors[0]==1); // Cannot mpi distribute in these dims
      assert( grid->_processors[1]==1);

      int Patch = Coor[nd-1];
      int HemiPatch = Patch%HemiPatches;
      int Hemisphere= Patch/HemiPatches;
      int north = Patch/HemiPatches;
      int south = 1-north;
      int isPoleY;
      int isPoleX;
      
      assert(Patch<IcosahedralPatches);
      assert((north==1)||(south==1));

      Coordinate XpCoor;
      Coordinate YpCoor;
      Coordinate DpCoor;
      Coordinate TpCoor;

      GetNbrForPlusDiagonal(grid,Coor,DpCoor);
      GetNbrForPlusX(grid,Coor,XpCoor,isPoleX);
      GetNbrForPlusY(grid,Coor,YpCoor,isPoleY);

      int DpPatch      = DpCoor[nd-1];
      int DpHemiPatch  = DpCoor[nd-1]%HemiPatches;
      int DpHemisphere = DpCoor[nd-1]/HemiPatches;

      SE._is_local     = true;
      SE._polarisation = 0;
      SE._adjoint      = false;
      SE._missing_link = false;
      SE._permute=0;
      //////////////////////////////////////////////////
      // Forward one site in time direction
      //////////////////////////////////////////////////
      int tdim  = 2;
      int delta = 1;
      int permute;
      GetOrthogNbr(grid,Coor,TpCoor,tdim,delta,permute);
      SE._offset = vertexgrid->oIndex(TpCoor);
      SE._permute=permute;
      //      std::cout << " Plaq stencil "<<Coor<<"  Tp "<<TpCoor<<" perm "<<permute<<std::endl;
      acceleratorPut(this->_entries[lexT],SE);
      SE._permute=0;
      ////////////////////////////////////////////////
      // X+ direction
      ////////////////////////////////////////////////
      if ( isPoleX ) {
	SE._offset     = vertexgrid->PoleSiteForOcoor(Coor);
      }	else {
	SE._offset     = vertexgrid->oIndex(XpCoor);
      }
      acceleratorPut(this->_entries[lexX],SE);
      ////////////////////////////////////////////////
      // Y+ direction
      ////////////////////////////////////////////////
      if ( isPoleY ) {
	SE._offset     = vertexgrid->PoleSiteForOcoor(Coor);
      }	else {
	SE._offset     = vertexgrid->oIndex(YpCoor);
      }
      acceleratorPut(this->_entries[lexY],SE);
      ////////////////////////////////////////////////
      // D+ direction
      ////////////////////////////////////////////////
      SE._offset     = vertexgrid->oIndex(DpCoor);
      //      std::cout << "Coor "<<Coor<<" DpCoor "<<DpCoor<<" site "<<SE._offset<<std::endl;
      acceleratorPut(this->_entries[lexD],SE);
    };
  }

  /*
   * enough to build staples in ico-T dir
   */
  /*
   * For gauge action derivative implementation
   * Staple 
   *
   * Case1: I x T loops
   * There is no complex rotation of links on other site
   *
   * Case2: I x I loops
   * Just use a general 6 point stencil and cherry pick terms
   */
  void TemporalStapleStencil(void)
  {
    GridBase * grid = this->_grid; // the edge grid
    GridBase * vertexgrid = this->_vertexgrid;

    int osites  = grid->oSites();

    uint64_t cart_sites = grid->CartesianOsites();
    uint64_t Npole_sites = grid->NorthPoleOsites();
    uint64_t Spole_sites = grid->SouthPoleOsites();
    Coordinate pcoor = grid->ThisProcessorCoor();
    Coordinate pgrid = grid->ProcessorGrid();

    /*
     * resize the stencil entries array and set npoints
     */
    this->_npoints=14; 
    this->_entries.resize(this->_npoints * cart_sites);
    this->_entries_p = &_entries[0];
    int np= this->_npoints;

    for(uint64_t site=0;site<cart_sites; site ++) {

      int nd = grid->Nd();
      int L  = grid->LocalDimensions()[0];

      Coordinate Coor;
      IcosahedralStencilEntry SE;

      Integer lexTp   = site*np+0;
      Integer lexXp   = site*np+1;
      Integer lexYp   = site*np+2;
      Integer lexDp   = site*np+3;
      Integer lexTm   = site*np+4;
      Integer lexTmXp = site*np+5;
      Integer lexTmYp = site*np+6;
      Integer lexTmDp = site*np+7;

      Integer lexXm   = site*np+8;
      Integer lexYm   = site*np+9;
      Integer lexDm   = site*np+10; // If !missingLink
      Integer lexXmTp = site*np+11;
      Integer lexYmTp = site*np+12;
      Integer lexDmTp = site*np+13; // If !missingLink
      
      ////////////////////////////////////////////////
      // Outer index of neighbour Offset calculation
      ////////////////////////////////////////////////
      grid->oCoorFromOindex(Coor,site);
      assert( grid->LocalDimensions()[1]==grid->LocalDimensions()[0]);
      assert( grid->_simd_layout[0]==1); // Cannot vectorise in these dims
      assert( grid->_simd_layout[1]==1);
      assert( grid->_processors[0]==1); // Cannot mpi distribute in these dims
      assert( grid->_processors[1]==1);

      int Patch = Coor[nd-1];
      int HemiPatch  = Coor[nd-1]%HemiPatches;
      int Hemisphere = Coor[nd-1]/HemiPatches;
      int north = Patch/HemiPatches;
      int south = 1-north;
      int isPoleY;
      int isPoleX;
      
      assert(Patch<IcosahedralPatches);
      assert((north==1)||(south==1));

      Coordinate XpCoor;
      Coordinate YpCoor;
      Coordinate DpCoor;
      Coordinate TpCoor;
      Coordinate TmCoor;
      Coordinate TmXpCoor;
      Coordinate TmYpCoor;
      Coordinate TmDpCoor;

      int missingLink;
      Coordinate XmCoor;
      Coordinate YmCoor;
      Coordinate DmCoor;
      Coordinate XmTpCoor;
      Coordinate YmTpCoor;
      Coordinate DmTpCoor;


      GetNbrForPlusDiagonal(grid,Coor,DpCoor);
      GetNbrForPlusX(grid,Coor,XpCoor,isPoleX);
      GetNbrForPlusY(grid,Coor,YpCoor,isPoleY);

      GetNbrForMinusDiagonal(grid,Coor,DmCoor,missingLink);
      GetNbrForMinusX(grid,Coor,XmCoor);
      GetNbrForMinusY(grid,Coor,YmCoor);

      
      int DpPatch      = DpCoor[nd-1];
      int DpHemiPatch  = DpCoor[nd-1]%HemiPatches;
      int DpHemisphere = DpCoor[nd-1]/HemiPatches;

      int XmHemiPatch  = XmCoor[nd-1]%HemiPatches;
      int XmHemisphere = XmCoor[nd-1]/HemiPatches;
      int YmHemiPatch  = YmCoor[nd-1]%HemiPatches;
      int YmHemisphere = YmCoor[nd-1]/HemiPatches;

      int DmPatch      = DmCoor[nd-1];
      int DmHemiPatch  = DmCoor[nd-1]%HemiPatches;
      int DmHemisphere = DmCoor[nd-1]/HemiPatches;
      
      SE._is_local     = true;
      SE._polarisation = 0;
      SE._adjoint      = false;
      SE._missing_link = false;
      SE._permute=0;
      //////////////////////////////////////////////////
      // Forward one site in time direction
      //////////////////////////////////////////////////
      int tdim  = 2;
      int delta = 1;
      int permute;
      GetOrthogNbr(grid,Coor,TpCoor,tdim,delta,permute);
      SE._offset = vertexgrid->oIndex(TpCoor);
      SE._permute= permute;
      acceleratorPut(this->_entries[lexTp],SE);
      SE._permute=0;
      ////////////////////////////////////////////////
      // X+ direction
      ////////////////////////////////////////////////
      if ( isPoleX ) {
	SE._offset     = vertexgrid->PoleSiteForOcoor(Coor);
      }	else {
	SE._offset     = vertexgrid->oIndex(XpCoor);
      }
      acceleratorPut(this->_entries[lexXp],SE);
      ////////////////////////////////////////////////
      // Y+ direction
      ////////////////////////////////////////////////
      if ( isPoleY ) {
	SE._offset     = vertexgrid->PoleSiteForOcoor(Coor);
      }	else {
	SE._offset     = vertexgrid->oIndex(YpCoor);
      }
      acceleratorPut(this->_entries[lexYp],SE);
      ////////////////////////////////////////////////
      // D+ direction
      ////////////////////////////////////////////////
      SE._offset     = vertexgrid->oIndex(DpCoor);
      acceleratorPut(this->_entries[lexDp],SE);
      //////////////////////////////////////////////////
      // Backward one site in time direction
      //////////////////////////////////////////////////
      GetOrthogNbr(grid,Coor,TmCoor,tdim,-delta,permute);
      SE._permute = permute;
      SE._offset = vertexgrid->oIndex(TmCoor);
      acceleratorPut(this->_entries[lexTm],SE);
      ////////////////////////////////////////////////
      // T-X+ hop
      // T-Y+ hop
      // T-D+ hop
      ////////////////////////////////////////////////
      GetNbrForPlusX(grid,TmCoor,TmXpCoor,isPoleX);
      GetNbrForPlusY(grid,TmCoor,TmYpCoor,isPoleY);
      GetNbrForPlusDiagonal(grid,TmCoor,TmDpCoor);
      if ( isPoleX ) {
	SE._offset     = vertexgrid->PoleSiteForOcoor(TmCoor);
      }	else {
	SE._offset     = vertexgrid->oIndex(TmXpCoor);
      }
      acceleratorPut(this->_entries[lexTmXp],SE);

      if ( isPoleY ) {
	SE._offset     = vertexgrid->PoleSiteForOcoor(TmCoor);
      }	else {
	SE._offset     = vertexgrid->oIndex(TmYpCoor);
      }
      acceleratorPut(this->_entries[lexTmYp],SE);

      SE._offset     = vertexgrid->oIndex(TmDpCoor);
      acceleratorPut(this->_entries[lexTmDp],SE);
      
      ////////////////////////////////////////////////
      // Links for the negative XYD staple on the forward T link
      ////////////////////////////////////////////////
      SE._permute = 0;
      SE._offset     = vertexgrid->oIndex(XmCoor);
      SE._polarisation = IcosahedronPatchX;
      if ( XmHemiPatch != HemiPatch && north ) {
	SE._polarisation = IcosahedronPatchDiagonal; 
      }
      acceleratorPut(this->_entries[lexXm],SE);
      GetOrthogNbr(grid,XmCoor,XmTpCoor,tdim,delta,permute);
      SE._permute = permute;
      SE._offset     = vertexgrid->oIndex(XmTpCoor);
      acceleratorPut(this->_entries[lexXmTp],SE);

      SE._permute = 0;
      SE._offset     = vertexgrid->oIndex(YmCoor);
      SE._polarisation = IcosahedronPatchY;
      if ( YmHemiPatch != HemiPatch && south ) {
	SE._polarisation = IcosahedronPatchDiagonal; 
      }
      acceleratorPut(this->_entries[lexYm],SE);
      GetOrthogNbr(grid,YmCoor,YmTpCoor,tdim,delta,permute);
      SE._permute = permute;
      SE._offset     = vertexgrid->oIndex(YmTpCoor);
      acceleratorPut(this->_entries[lexYmTp],SE);

      if ( ! missingLink ) {

	SE._polarisation = IcosahedronPatchDiagonal;   // Basis rotates
	if ( (DmHemiPatch != HemiPatch) && (DmHemisphere==Hemisphere) && south ) {
	  SE._polarisation = IcosahedronPatchX;   // Basis rotates
	}
	if ( (DmHemiPatch != HemiPatch) && (DmHemisphere==Hemisphere) && north ) {
	  SE._polarisation = IcosahedronPatchY;   // Basis rotates
	}
	
	SE._permute = 0;
	SE._offset     = vertexgrid->oIndex(DmCoor);
	acceleratorPut(this->_entries[lexDm],SE);
	
	GetOrthogNbr(grid,DmCoor,DmTpCoor,tdim,delta,permute);
	SE._permute = permute;
	SE._offset     = vertexgrid->oIndex(DmTpCoor);
	acceleratorPut(this->_entries[lexDmTp],SE);

      } else {
	
	SE._permute = 0;
	SE._offset  = 0;
	SE._missing_link = 1;
	acceleratorPut(this->_entries[lexDm],SE);
	acceleratorPut(this->_entries[lexDmTp],SE);

      }
    };

    int ndm1 = grid->Nd()-1;
    int L    = grid->LocalDimensions()[0];
    if ( vertexgrid->ownsSouthPole() ) {
      IcosahedralStencilEntry SE;
      for(uint64_t site=0;site<cart_sites; site ++) { // loops over volume

	Coordinate Coor;
	Coordinate tCoor;
	vertexgrid->oCoorFromOindex(Coor,site);
	int north = Coor[ndm1]/HemiPatches;
	int south = 1-north;
	int tdim = 2;

	if( (Coor[0]==(L-1))&&(Coor[1]==0)&&south ) {

	  int64_t pole_site = vertexgrid->PoleSiteForOcoor(Coor);
	  int64_t lex       = pole_site*np+(Coor[ndm1]%HemiPatches)*2;

	  SE._offset       = site;
	  SE._is_local     = true;
	  SE._polarisation = IcosahedronPatchX;  // ignored
	  SE._adjoint      = false;              // ignored
	  SE._missing_link = false;
	  SE._permute=0;
	  acceleratorPut(this->_entries[lex],SE);
	    
	  lex       = pole_site*np+(Coor[ndm1]%HemiPatches)*2+1;

	  int Rt = vertexgrid->_rdimensions[tdim];
	  int permute;
	  tCoor = Coor;
	  tCoor[tdim] = periAdd(tCoor[tdim],1,Rt,permute);
	  SE._offset  = vertexgrid->oIndex(tCoor);
	  SE._permute= permute;
	  acceleratorPut(this->_entries[lex],SE);
	    
	}
      }
    }
    if ( vertexgrid->ownsNorthPole() ) {
      IcosahedralStencilEntry SE;
      for(uint64_t site=0;site<cart_sites; site ++) {
	Coordinate Coor;
	Coordinate tCoor;
	vertexgrid->oCoorFromOindex(Coor,site);
	int north = Coor[ndm1]/HemiPatches;
	int tdim = 2;
	
	if( (Coor[0]==0)&&(Coor[1]==(L-1))&&north ) {

	  int64_t pole_site = vertexgrid->PoleSiteForOcoor(Coor);
	  int64_t lex       = pole_site*np+(Coor[ndm1]%HemiPatches)*2;

	  SE._offset       = site;
	  SE._is_local     = true;
	  SE._polarisation = IcosahedronPatchY;  // ignored
	  SE._adjoint      = false;              // ignored
	  SE._missing_link = false;
	  SE._permute=0;
	  acceleratorPut(this->_entries[lex],SE);
	    
	  lex       = pole_site*np+(Coor[ndm1]%HemiPatches)*2+1;

	  int Rt = vertexgrid->_rdimensions[tdim];
	  int permute;
	  tCoor = Coor;
	  tCoor[tdim] = periAdd(tCoor[tdim],1,Rt,permute);
	  SE._offset  = vertexgrid->oIndex(tCoor);
	  SE._permute= permute;
	  acceleratorPut(this->_entries[lex],SE);
	    
	}
      }
    }
  }
};
NAMESPACE_END(Grid);
  
