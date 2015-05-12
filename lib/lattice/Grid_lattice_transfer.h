#ifndef GRID_LATTICE_TRANSFER_H
#define GRID_LATTICE_TRANSFER_H

namespace Grid {

inline void subdivides(GridBase *coarse,GridBase *fine)
{
  assert(coarse->_ndimension == fine->_ndimension);

  int _ndimension = coarse->_ndimension;

  // local and global volumes subdivide cleanly after SIMDization
  for(int d=0;d<_ndimension;d++){
    assert(coarse->_processors[d]  == fine->_processors[d]);
    assert(coarse->_simd_layout[d] == fine->_simd_layout[d]);
    assert((fine->_rdimensions[d] / coarse->_rdimensions[d])* coarse->_rdimensions[d]==fine->_rdimensions[d]); 
  }
}

  ////////////////////////////////////////////////////////////////////////////////////////////
  // remove and insert a half checkerboard
  ////////////////////////////////////////////////////////////////////////////////////////////
  template<class vobj> inline void pickCheckerboard(int cb,Lattice<vobj> &half,const Lattice<vobj> &full){
    half.checkerboard = cb;
    int ssh=0;
PARALLEL_FOR_LOOP
    for(int ss=0;ss<full._grid->oSites();ss++){
      std::vector<int> coor;
      int cbos;
      
      full._grid->oCoorFromOindex(coor,ss);
      cbos=half._grid->CheckerBoard(coor);
      
      if (cbos==cb) {
	
	half._odata[ssh] = full._odata[ss];
	ssh++;
      }
    }
  }
  template<class vobj> inline void setCheckerboard(Lattice<vobj> &full,const Lattice<vobj> &half){
    int cb = half.checkerboard;
    int ssh=0;
PARALLEL_FOR_LOOP
    for(int ss=0;ss<full._grid->oSites();ss++){
      std::vector<int> coor;
      int cbos;
      
      full._grid->oCoorFromOindex(coor,ss);
      cbos=half._grid->CheckerBoard(coor);
      
      if (cbos==cb) {
	full._odata[ss]=half._odata[ssh];
	ssh++;
      }
    }
  }
  

template<class vobj,int nbasis>
inline void projectBlockBasis(Lattice<iVector<vComplex,nbasis > > &coarseData,
			      const             Lattice<vobj>   &fineData,
			      const std::vector<Lattice<vobj> > &Basis)
{
  GridBase * fine  = fineData._grid;
  GridBase * coarse= coarseData._grid;
  int  _ndimension = coarse->_ndimension;

  // checks
  assert( nbasis == Basis.size() );
  subdivides(coarse,fine); 
  for(int i=0;i<nbasis;i++){
    conformable(Basis,fineData);
  }

  std::vector<int>  block_r      (_ndimension);
  
  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
  }

  coarseData=zero;

  // Loop with a cache friendly loop ordering
  for(int sf=0;sf<fine->oSites();sf++){

    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);
    GridBase::CoorFromIndex(coor_f,sf,fine->_rdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    GridBase::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

    for(int i=0;i<nbasis;i++) {
      
      coarseData._odata[sc][i]=coarseData._odata[sc][i]
	+ innerProduct(Basis[i]._odata[sf],fineData._odata[sf]);

    }
  }
  return;
  
}


template<class vobj,int nbasis>
inline void promoteBlockBasis(const Lattice<iVector<vComplex,nbasis > > &coarseData,
			                        Lattice<vobj>   &fineData,
			      const std::vector<Lattice<vobj> > &Basis)
{
  GridBase * fine  = fineData._grid;
  GridBase * coarse= coarseData._grid;
  int  _ndimension = coarse->_ndimension;

  // checks
  assert( nbasis == Basis.size() );
  subdivides(coarse,fine); 
  for(int i=0;i<nbasis;i++){
    conformable(Basis,fineData);
  }

  std::vector<int>  block_r      (_ndimension);
  
  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
  }

  // Loop with a cache friendly loop ordering
  for(int sf=0;sf<fine->oSites();sf++){

    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);

    GridBase::CoorFromIndex(coor_f,sf,fine->_rdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    GridBase::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

    for(int i=0;i<nbasis;i++) {

      if(i==0) fineData._odata[sf]=                    coarseData._odata[sc][i]*Basis[i]._odata[sf];
      else     fineData._odata[sf]=fineData._odata[sf]+coarseData._odata[sc][i]*Basis[i]._odata[sf];

    }
  }
  return;
  
}

// useful in multigrid project;
// Generic name : Coarsen?
template<class vobj>
inline void sumBlocks(Lattice<vobj> &coarseData,const Lattice<vobj> &fineData)
{
  GridBase * fine  = fineData._grid;
  GridBase * coarse= coarseData._grid;

  subdivides(coarse,fine); // require they map

  int _ndimension = coarse->_ndimension;
  
  std::vector<int>  block_r      (_ndimension);
  
  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
  }

  coarseData=zero;
  for(int sf=0;sf<fine->oSites();sf++){
    
    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);

    GridBase::CoorFromIndex(coor_f,sf,fine->_rdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    GridBase::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

    coarseData._odata[sc]=coarseData._odata[sc]+fineData._odata[sf];

  }
  return;
}

}
#endif
