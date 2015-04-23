#ifndef GRID_SUMMATION_H
#define GRID_SUMMATION_H

void subdivides(GridBase *coarse,GridBase *fine)
{
  assert(coarse->_ndimension == fine->_ndimension);

  int _ndimension = coarse->_ndimension;

  // local and global volumes subdivide cleanly after SIMDization
  for(int d=0;d<_ndimension;d++){
    assert((fine->_rdimensions[d] / coarse->_rdimensions[d])* coarse->_rdimensions[d]==fine->_rdimensions[d]); 
    assert((fine->_ldimensions[d] / coarse->_ldimensions[d])* coarse->_ldimensions[d]==fine->_ldimensions[d]); 
    assert((fine->_gdimensions[d] / coarse->_gdimensions[d])* coarse->_gdimensions[d]==fine->_gdimensions[d]); 
    assert((fine->_fdimensions[d] / coarse->_fdimensions[d])* coarse->_fdimensions[d]==fine->_fdimensions[d]); 
  }
}
// Generic name : Coarsen?
//              : SubMeshSum?
//
template<class vobj>
inline void sumBlocks(Lattice<vobj> &coarseData,const Lattice<vobj> &fineData)
{
  GridBase * fine  = findData._grid;
  GridBase * coarse= findData._grid;

  subdivides(coars,fine); // require they map

  int _ndimension = coarse->_ndimension;
  
  std::vector<bool> replicated(_ndimension,false);
  std::vector<int>  block_r   (_dimension);
  std::vector<int>  block_f   (_dimension);

  ///////////////////////////////////////////////////////////
  // Detect whether the result is replicated in dimension d
  ///////////////////////////////////////////////////////////
  for(int d=0 ; d<_ndimension;d++){
    if ( (_fdimensions[d] == 1) && (coarse->_processors[d]>1) ) {
      replicated[d]=true;
    }
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
    block_l[d] = fine->_ldimensions[d] / coarse->_ldimensions[d];
    block_f[d] = fine->_fdimensions[d] / coarse->_fdimensions[d];
  }

  coaseData=zero;

  //FIXME Bagel's strategy: loop over fine sites
  // identify corresponding coarse site, but coarse sites are 
  // divided across threads. Not so easy to do in openmp but
  // there must be a way
  for(int sf=0;sf<fine->oSites();sf++){
    
    int sc;
    vobj sum=zero;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);

    GridBase::CoorFromIndex(coor_f,sf,fine->_rdimensions);

    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/fine->_rdimensions;

    GridBase::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

    coarseData._odata[sc]=coarseData._odata[sc]+fineData._odata[sf];

  }
  return;
}
#endif
