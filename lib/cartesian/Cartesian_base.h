    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cartesian/Cartesian_base.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_CARTESIAN_BASE_H
#define GRID_CARTESIAN_BASE_H

#include <Grid.h>

namespace Grid{

  //////////////////////////////////////////////////////////////////////
  // Commicator provides information on the processor grid
  //////////////////////////////////////////////////////////////////////
  //    unsigned long _ndimension;
  //    std::vector<int> _processors; // processor grid
  //    int              _processor;  // linear processor rank
  //    std::vector<int> _processor_coor;  // linear processor rank
  //////////////////////////////////////////////////////////////////////
  class GridBase : public CartesianCommunicator , public GridThread {

public:

    // Give Lattice access
    template<class object> friend class Lattice;

    GridBase(const std::vector<int> & processor_grid) : CartesianCommunicator(processor_grid) {};


    // Physics Grid information.
    std::vector<int> _simd_layout;// Which dimensions get relayed out over simd lanes.
    std::vector<int> _fdimensions;// Global dimensions of array prior to cb removal
    std::vector<int> _gdimensions;// Global dimensions of array after cb removal
    std::vector<int> _ldimensions;// local dimensions of array with processor images removed
    std::vector<int> _rdimensions;// Reduced local dimensions with simd lane images and processor images removed 
    std::vector<int> _ostride;    // Outer stride for each dimension
    std::vector<int> _istride;    // Inner stride i.e. within simd lane
    int _osites;                  // _isites*_osites = product(dimensions).
    int _isites;
    int _fsites;                  // _isites*_osites = product(dimensions).
    int _gsites;
    std::vector<int> _slice_block;   // subslice information
    std::vector<int> _slice_stride;
    std::vector<int> _slice_nblock;

    // Might need these at some point
    //    std::vector<int> _lstart;     // local start of array in gcoors. _processor_coor[d]*_ldimensions[d]
    //    std::vector<int> _lend;       // local end of array in gcoors    _processor_coor[d]*_ldimensions[d]+_ldimensions_[d]-1

public:

    ////////////////////////////////////////////////////////////////
    // Checkerboarding interface is virtual and overridden by 
    // GridCartesian / GridRedBlackCartesian
    ////////////////////////////////////////////////////////////////
    virtual int CheckerBoarded(int dim)=0;
    virtual int CheckerBoard(std::vector<int> site)=0;
    virtual int CheckerBoardDestination(int source_cb,int shift,int dim)=0;
    virtual int CheckerBoardShift(int source_cb,int dim,int shift,int osite)=0;
    virtual int CheckerBoardShiftForCB(int source_cb,int dim,int shift,int cb)=0;
    int  CheckerBoardFromOindex (int Oindex){
      std::vector<int> ocoor;
      oCoorFromOindex(ocoor,Oindex); 
      return CheckerBoard(ocoor);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Local layout calculations
    //////////////////////////////////////////////////////////////////////////////////////////////
    // These routines are key. Subdivide the linearised cartesian index into
    //      "inner" index identifying which simd lane of object<vFcomplex> is associated with coord
    //      "outer" index identifying which element of _odata in class "Lattice" is associated with coord.
    //
    // Compared to, say, Blitz++ we simply need to store BOTH an inner stride and an outer
    // stride per dimension. The cost of evaluating the indexing information is doubled for an n-dimensional
    // coordinate. Note, however, for data parallel operations the "inner" indexing cost is not paid and all
    // lanes are operated upon simultaneously.
  
    virtual int oIndex(std::vector<int> &coor)
    {
        int idx=0;
	// Works with either global or local coordinates
        for(int d=0;d<_ndimension;d++) idx+=_ostride[d]*(coor[d]%_rdimensions[d]);
        return idx;
    }
    inline int oIndexReduced(std::vector<int> &ocoor)
    {
      int idx=0; 
      // ocoor is already reduced so can eliminate the modulo operation
      // for fast indexing and inline the routine
      for(int d=0;d<_ndimension;d++) idx+=_ostride[d]*ocoor[d];
      return idx;
    }
    inline void oCoorFromOindex (std::vector<int>& coor,int Oindex){
      Lexicographic::CoorFromIndex(coor,Oindex,_rdimensions);
    }


    //////////////////////////////////////////////////////////
    // SIMD lane addressing
    //////////////////////////////////////////////////////////
    inline int iIndex(std::vector<int> &lcoor)
    {
        int idx=0;
        for(int d=0;d<_ndimension;d++) idx+=_istride[d]*(lcoor[d]/_rdimensions[d]);
        return idx;
    }
    inline void iCoorFromIindex(std::vector<int> &coor,int lane)
    {
      Lexicographic::CoorFromIndex(coor,lane,_simd_layout);
    }
    inline int PermuteDim(int dimension){
      return _simd_layout[dimension]>1;
    }
    inline int PermuteType(int dimension){
      int permute_type=0;
      //
      // FIXME:
      //
      // Best way to encode this would be to present a mask 
      // for which simd dimensions are rotated, and the rotation
      // size. If there is only one simd dimension rotated, this is just 
      // a permute. 
      //
      // Cases: PermuteType == 1,2,4,8
      // Distance should be either 0,1,2..
      //
      if ( _simd_layout[dimension] > 2 ) { 
	for(int d=0;d<_ndimension;d++){
	  if ( d != dimension ) assert ( (_simd_layout[d]==1)  );
	}
	permute_type = RotateBit; // How to specify distance; this is not just direction.
	return permute_type;
      }

      for(int d=_ndimension-1;d>dimension;d--){
	if (_simd_layout[d]>1 ) permute_type++;
      }
      return permute_type;
    }
    ////////////////////////////////////////////////////////////////
    // Array sizing queries
    ////////////////////////////////////////////////////////////////

    inline int iSites(void) const { return _isites; };
    inline int Nsimd(void)  const { return _isites; };// Synonymous with iSites
    inline int oSites(void) const { return _osites; };
    inline int lSites(void) const { return _isites*_osites; }; 
    inline int gSites(void) const { return _isites*_osites*_Nprocessors; }; 
    inline int Nd    (void) const { return _ndimension;};

    inline const std::vector<int> &FullDimensions(void)         { return _fdimensions;};
    inline const std::vector<int> &GlobalDimensions(void)       { return _gdimensions;};
    inline const std::vector<int> &LocalDimensions(void)        { return _ldimensions;};
    inline const std::vector<int> &VirtualLocalDimensions(void) { return _ldimensions;};

    ////////////////////////////////////////////////////////////////
    // Global addressing
    ////////////////////////////////////////////////////////////////
    void GlobalIndexToGlobalCoor(int gidx,std::vector<int> &gcoor){
      Lexicographic::CoorFromIndex(gcoor,gidx,_gdimensions);
    }
    void LocalIndexToLocalCoor(int lidx,std::vector<int> &lcoor){
      Lexicographic::CoorFromIndex(lcoor,lidx,_ldimensions);
    }
    void GlobalCoorToGlobalIndex(const std::vector<int> & gcoor,int & gidx){
      gidx=0;
      int mult=1;
      for(int mu=0;mu<_ndimension;mu++) {
	gidx+=mult*gcoor[mu];
	mult*=_gdimensions[mu];
      }
    }
    void GlobalCoorToProcessorCoorLocalCoor(std::vector<int> &pcoor,std::vector<int> &lcoor,const std::vector<int> &gcoor)
    {
      pcoor.resize(_ndimension);
      lcoor.resize(_ndimension);
      for(int mu=0;mu<_ndimension;mu++){
	int _fld  = _fdimensions[mu]/_processors[mu];
	pcoor[mu] = gcoor[mu]/_fld;
	lcoor[mu] = gcoor[mu]%_fld;
      }
    }
    void GlobalCoorToRankIndex(int &rank, int &o_idx, int &i_idx ,const std::vector<int> &gcoor)
    {
      std::vector<int> pcoor;
      std::vector<int> lcoor;
      GlobalCoorToProcessorCoorLocalCoor(pcoor,lcoor,gcoor);
      rank = RankFromProcessorCoor(pcoor);

      std::vector<int> cblcoor(lcoor);
      for(int d=0;d<cblcoor.size();d++){
	if( this->CheckerBoarded(d) ) {
	  cblcoor[d] = lcoor[d]/2;
	}
      }

      i_idx= iIndex(cblcoor);// this does not imply divide by 2 on checker dim
      o_idx= oIndex(lcoor);// this implies divide by 2 on checkerdim
    }

    void RankIndexToGlobalCoor(int rank, int o_idx, int i_idx , std::vector<int> &gcoor)
    {
      gcoor.resize(_ndimension);
      std::vector<int> coor(_ndimension);

      ProcessorCoorFromRank(rank,coor);
      for(int mu=0;mu<_ndimension;mu++) gcoor[mu] = _ldimensions[mu]*coor[mu];

      iCoorFromIindex(coor,i_idx);
      for(int mu=0;mu<_ndimension;mu++) gcoor[mu] += _rdimensions[mu]*coor[mu];

      oCoorFromOindex (coor,o_idx);
      for(int mu=0;mu<_ndimension;mu++) gcoor[mu] += coor[mu];
      
    }
    void RankIndexCbToFullGlobalCoor(int rank, int o_idx, int i_idx, int cb,std::vector<int> &fcoor)
    {
      RankIndexToGlobalCoor(rank,o_idx,i_idx ,fcoor);
      if(CheckerBoarded(0)){
	fcoor[0] = fcoor[0]*2+cb;
      }
    }
    void ProcessorCoorLocalCoorToGlobalCoor(std::vector<int> &Pcoor,std::vector<int> &Lcoor,std::vector<int> &gcoor)
    {
      gcoor.resize(_ndimension);
      for(int mu=0;mu<_ndimension;mu++) gcoor[mu] = Pcoor[mu]*_ldimensions[mu]+Lcoor[mu];
    }
};


}
#endif
