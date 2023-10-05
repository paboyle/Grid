/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cartesian/Cartesian_base.h

    Copyright (C) 2015

    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: paboyle <paboyle@ph.ed.ac.uk>
    Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////
// Commicator provides information on the processor grid
//////////////////////////////////////////////////////////////////////
//    unsigned long _ndimension;
//    Coordinate _processors; // processor grid
//    int              _processor;  // linear processor rank
//    Coordinate _processor_coor;  // linear processor rank
//////////////////////////////////////////////////////////////////////
class GridBase : public CartesianCommunicator , public GridThread {

public:
  int dummy;
  // Give Lattice access
  template<class object> friend class Lattice;

  GridBase(const Coordinate & processor_grid) : CartesianCommunicator(processor_grid) { LocallyPeriodic=0;}; 

  GridBase(const Coordinate & processor_grid,
	   const CartesianCommunicator &parent,
	   int &split_rank) 
    : CartesianCommunicator(processor_grid,parent,split_rank) {LocallyPeriodic=0;};

  GridBase(const Coordinate & processor_grid,
	   const CartesianCommunicator &parent) 
    : CartesianCommunicator(processor_grid,parent,dummy) {LocallyPeriodic=0;};

  virtual ~GridBase() = default;

  // Physics Grid information.
  Coordinate _simd_layout;// Which dimensions get relayed out over simd lanes.
  Coordinate _fdimensions;// (full) Global dimensions of array prior to cb removal
  Coordinate _gdimensions;// Global dimensions of array after cb removal
  Coordinate _ldimensions;// local dimensions of array with processor images removed
  Coordinate _rdimensions;// Reduced local dimensions with simd lane images and processor images removed 
  Coordinate _ostride;    // Outer stride for each dimension
  Coordinate _istride;    // Inner stride i.e. within simd lane
  int _osites;                  // _isites*_osites = product(dimensions).
  int _isites;
  int64_t _fsites;                  // _isites*_osites = product(dimensions).
  int64_t _gsites;
  Coordinate _slice_block;// subslice information
  Coordinate _slice_stride;
  Coordinate _slice_nblock;

  Coordinate _lstart;     // local start of array in gcoors _processor_coor[d]*_ldimensions[d]
  Coordinate _lend  ;     // local end of array in gcoors   _processor_coor[d]*_ldimensions[d]+_ldimensions_[d]-1

  bool _isCheckerBoarded; 
  int        LocallyPeriodic;
  Coordinate _checker_dim_mask;

public:

  ////////////////////////////////////////////////////////////////
  // Checkerboarding interface is virtual and overridden by 
  // GridCartesian / GridRedBlackCartesian
  ////////////////////////////////////////////////////////////////
  virtual int CheckerBoarded(int dim)=0;
  virtual int CheckerBoard(const Coordinate &site)=0;
  virtual int CheckerBoardDestination(int source_cb,int shift,int dim)=0;
  virtual int CheckerBoardShift(int source_cb,int dim,int shift,int osite)=0;
  virtual int CheckerBoardShiftForCB(int source_cb,int dim,int shift,int cb)=0;
  virtual int CheckerBoardFromOindex (int Oindex)=0;
  virtual int CheckerBoardFromOindexTable (int Oindex)=0;

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
  
  virtual int oIndex(Coordinate &coor)
  {
    int idx=0;
    // Works with either global or local coordinates
    for(int d=0;d<_ndimension;d++) idx+=_ostride[d]*(coor[d]%_rdimensions[d]);
    return idx;
  }
  virtual int iIndex(Coordinate &lcoor)
  {
    int idx=0;
    for(int d=0;d<_ndimension;d++) idx+=_istride[d]*(lcoor[d]/_rdimensions[d]);
    return idx;
  }
  inline int oIndexReduced(Coordinate &ocoor)
  {
    int idx=0; 
    // ocoor is already reduced so can eliminate the modulo operation
    // for fast indexing and inline the routine
    for(int d=0;d<_ndimension;d++) idx+=_ostride[d]*ocoor[d];
    return idx;
  }
  inline void oCoorFromOindex (Coordinate& coor,int Oindex){
    Lexicographic::CoorFromIndex(coor,Oindex,_rdimensions);
  }

  inline void InOutCoorToLocalCoor (Coordinate &ocoor, Coordinate &icoor, Coordinate &lcoor) {
    lcoor.resize(_ndimension);
    for (int d = 0; d < _ndimension; d++)
      lcoor[d] = ocoor[d] + _rdimensions[d] * icoor[d];
  }

  //////////////////////////////////////////////////////////
  // SIMD lane addressing
  //////////////////////////////////////////////////////////
  inline void iCoorFromIindex(Coordinate &coor,int lane)
  {
    Lexicographic::CoorFromIndex(coor,lane,_simd_layout);
  }

  inline int PermuteDim(int dimension){
    return _simd_layout[dimension]>1;
  }
  inline int PermuteType(int dimension){
    int permute_type=0;
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
  inline int64_t gSites(void) const { return (int64_t)_isites*(int64_t)_osites*(int64_t)_Nprocessors; }; 
  inline int Nd    (void) const { return _ndimension;};

  inline const Coordinate LocalStarts(void)             { return _lstart;    };
  inline const Coordinate &FullDimensions(void)         { return _fdimensions;};
  inline const Coordinate &GlobalDimensions(void)       { return _gdimensions;};
  inline const Coordinate &LocalDimensions(void)        { return _ldimensions;};
  inline const Coordinate &VirtualLocalDimensions(void) { return _ldimensions;};

  ////////////////////////////////////////////////////////////////
  // Utility to print the full decomposition details 
  ////////////////////////////////////////////////////////////////

  void show_decomposition(){
    std::cout << GridLogMessage << "\tFull Dimensions    : " << _fdimensions << std::endl;
    std::cout << GridLogMessage << "\tSIMD layout        : " << _simd_layout << std::endl;
    std::cout << GridLogMessage << "\tGlobal Dimensions  : " << _gdimensions << std::endl;
    std::cout << GridLogMessage << "\tLocal Dimensions   : " << _ldimensions << std::endl;
    std::cout << GridLogMessage << "\tReduced Dimensions : " << _rdimensions << std::endl;
    std::cout << GridLogMessage << "\tOuter strides      : " << _ostride << std::endl;
    std::cout << GridLogMessage << "\tInner strides      : " << _istride << std::endl;
    std::cout << GridLogMessage << "\tiSites             : " << _isites << std::endl;
    std::cout << GridLogMessage << "\toSites             : " << _osites << std::endl;
    std::cout << GridLogMessage << "\tlSites             : " << lSites() << std::endl;        
    std::cout << GridLogMessage << "\tgSites             : " << gSites() << std::endl;
    std::cout << GridLogMessage << "\tNd                 : " << _ndimension << std::endl;             
  } 

  ////////////////////////////////////////////////////////////////
  // Global addressing
  ////////////////////////////////////////////////////////////////
  void GlobalIndexToGlobalCoor(int64_t gidx,Coordinate &gcoor){
    assert(gidx< gSites());
    Lexicographic::CoorFromIndex(gcoor,gidx,_gdimensions);
  }
  void LocalIndexToLocalCoor(int lidx,Coordinate &lcoor){
    assert(lidx<lSites());
    Lexicographic::CoorFromIndex(lcoor,lidx,_ldimensions);
  }
  void GlobalCoorToGlobalIndex(const Coordinate & gcoor,int64_t & gidx){
    gidx=0;
    int mult=1;
    for(int mu=0;mu<_ndimension;mu++) {
      gidx+=mult*gcoor[mu];
      mult*=_gdimensions[mu];
    }
  }
  void GlobalCoorToProcessorCoorLocalCoor(Coordinate &pcoor,Coordinate &lcoor,const Coordinate &gcoor)
  {
    pcoor.resize(_ndimension);
    lcoor.resize(_ndimension);
    for(int mu=0;mu<_ndimension;mu++){
      int _fld  = _fdimensions[mu]/_processors[mu];
      pcoor[mu] = gcoor[mu]/_fld;
      lcoor[mu] = gcoor[mu]%_fld;
    }
  }
  void GlobalCoorToRankIndex(int &rank, int &o_idx, int &i_idx ,const Coordinate &gcoor)
  {
    Coordinate pcoor;
    Coordinate lcoor;
    GlobalCoorToProcessorCoorLocalCoor(pcoor,lcoor,gcoor);
    rank = RankFromProcessorCoor(pcoor);
    /*
      Coordinate cblcoor(lcoor);
      for(int d=0;d<cblcoor.size();d++){
      if( this->CheckerBoarded(d) ) {
      cblcoor[d] = lcoor[d]/2;
      }
      }
    */
    i_idx= iIndex(lcoor);
    o_idx= oIndex(lcoor);
  }

  void RankIndexToGlobalCoor(int rank, int o_idx, int i_idx , Coordinate &gcoor)
  {
    gcoor.resize(_ndimension);
    Coordinate coor(_ndimension);

    ProcessorCoorFromRank(rank,coor);
    for(int mu=0;mu<_ndimension;mu++) gcoor[mu] = _ldimensions[mu]*coor[mu];

    iCoorFromIindex(coor,i_idx);
    for(int mu=0;mu<_ndimension;mu++) gcoor[mu] += _rdimensions[mu]*coor[mu];

    oCoorFromOindex (coor,o_idx);
    for(int mu=0;mu<_ndimension;mu++) gcoor[mu] += coor[mu];
      
  }
  void RankIndexCbToFullGlobalCoor(int rank, int o_idx, int i_idx, int cb,Coordinate &fcoor)
  {
    RankIndexToGlobalCoor(rank,o_idx,i_idx ,fcoor);
    if(CheckerBoarded(0)){
      fcoor[0] = fcoor[0]*2+cb;
    }
  }
  void ProcessorCoorLocalCoorToGlobalCoor(Coordinate &Pcoor,Coordinate &Lcoor,Coordinate &gcoor)
  {
    gcoor.resize(_ndimension);
    for(int mu=0;mu<_ndimension;mu++) gcoor[mu] = Pcoor[mu]*_ldimensions[mu]+Lcoor[mu];
  }
};

NAMESPACE_END(Grid);
#endif
