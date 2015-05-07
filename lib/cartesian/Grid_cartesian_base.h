#ifndef GRID_CARTESIAN_BASE_H
#define GRID_CARTESIAN_BASE_H

#include <Grid.h>
#include <Grid_communicator.h>

namespace Grid{

  //////////////////////////////////////////////////////////////////////
  // Commicator provides information on the processor grid
  //////////////////////////////////////////////////////////////////////
  //    unsigned long _ndimension;
  //    std::vector<int> _processors; // processor grid
  //    int              _processor;  // linear processor rank
  //    std::vector<int> _processor_coor;  // linear processor rank
  //////////////////////////////////////////////////////////////////////
  class GridBase : public CartesianCommunicator {

public:

    // Give Lattice access
    template<class object> friend class Lattice;

    GridBase(std::vector<int> & processor_grid) : CartesianCommunicator(processor_grid) {};
            
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
    virtual int CheckerBoardDestination(int source_cb,int shift)=0;
    virtual int CheckerBoardShift(int source_cb,int dim,int shift,int osite)=0;
    inline int  CheckerBoardFromOindex (int Oindex){
      std::vector<int> ocoor;
      oCoorFromOindex(ocoor,Oindex); 
      int ss=0;
      for(int d=0;d<_ndimension;d++){
	ss=ss+ocoor[d];
      }      
      return ss&0x1;
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
    static inline void CoorFromIndex (std::vector<int>& coor,int index,std::vector<int> &dims){
      int nd= dims.size();
      coor.resize(nd);
      for(int d=0;d<nd;d++){
	coor[d] = index % dims[d];
	index   = index / dims[d];
      }
    }
    inline void oCoorFromOindex (std::vector<int>& coor,int Oindex){
      CoorFromIndex(coor,Oindex,_rdimensions);
    }
    static inline void IndexFromCoor (std::vector<int>& coor,int &index,std::vector<int> &dims){
      int nd=dims.size();
      int stride=1;
      index=0;
      for(int d=0;d<nd;d++){
	index = index+stride*coor[d];
	stride=stride*dims[d];
      }
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
      CoorFromIndex(coor,lane,_simd_layout);
    }
    inline int PermuteDim(int dimension){
      return _simd_layout[dimension]>1;
    }
    inline int PermuteType(int dimension){
      int permute_type=0;
      for(int d=_ndimension-1;d>dimension;d--){
	if (_simd_layout[d]>1 ) permute_type++;
      }
      return permute_type;
    }
    ////////////////////////////////////////////////////////////////
    // Array sizing queries
    ////////////////////////////////////////////////////////////////

    inline int iSites(void) { return _isites; };
    inline int Nsimd(void)  { return _isites; };// Synonymous with iSites
    inline int oSites(void) { return _osites; };
    inline int lSites(void) { return _isites*_osites; }; 
    inline int gSites(void) { return _isites*_osites*_Nprocessors; }; 
    inline int Nd    (void) { return _ndimension;};

    inline const std::vector<int> &FullDimensions(void)         { return _fdimensions;};
    inline const std::vector<int> &GlobalDimensions(void)       { return _gdimensions;};
    inline const std::vector<int> &LocalDimensions(void)        { return _ldimensions;};
    inline const std::vector<int> &VirtualLocalDimensions(void) { return _ldimensions;};

    ////////////////////////////////////////////////////////////////
    // Global addressing
    ////////////////////////////////////////////////////////////////
    void GlobalIndexToGlobalCoor(int gidx,std::vector<int> &gcoor){
      CoorFromIndex(gcoor,gidx,_gdimensions);
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
	pcoor[mu] = gcoor[mu]/_ldimensions[mu];
	lcoor[mu] = gcoor[mu]%_ldimensions[mu];
      }
    }
    void GlobalCoorToRankIndex(int &rank, int &o_idx, int &i_idx ,const std::vector<int> &gcoor)
    {
      std::vector<int> pcoor;
      std::vector<int> lcoor;
      GlobalCoorToProcessorCoorLocalCoor(pcoor,lcoor,gcoor);
      rank = RankFromProcessorCoor(pcoor);
      i_idx= iIndex(lcoor);
      o_idx= oIndex(lcoor);
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
