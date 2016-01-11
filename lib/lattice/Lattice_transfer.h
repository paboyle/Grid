    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_transfer.h

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
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
    //PARALLEL_FOR_LOOP
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
    //PARALLEL_FOR_LOOP
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
  

template<class vobj,class CComplex,int nbasis>
inline void blockProject(Lattice<iVector<CComplex,nbasis > > &coarseData,
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
    conformable(Basis[i],fineData);
  }

  std::vector<int>  block_r      (_ndimension);
  
  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
    assert(block_r[d]*coarse->_rdimensions[d] == fine->_rdimensions[d]);
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
      
      coarseData._odata[sc](i)=coarseData._odata[sc](i)
	+ innerProduct(Basis[i]._odata[sf],fineData._odata[sf]);

    }
  }
  return;
}

template<class vobj,class CComplex>
inline void blockZAXPY(Lattice<vobj> &fineZ,
		       const Lattice<CComplex> &coarseA,
		       const Lattice<vobj> &fineX,
		       const Lattice<vobj> &fineY)
{
  GridBase * fine  = fineZ._grid;
  GridBase * coarse= coarseA._grid;

  fineZ.checkerboard=fineX.checkerboard;
  subdivides(coarse,fine); // require they map
  conformable(fineX,fineY);
  conformable(fineX,fineZ);

  int _ndimension = coarse->_ndimension;
  
  std::vector<int>  block_r      (_ndimension);

  // FIXME merge with subdivide checking routine as this is redundant
  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
    assert(block_r[d]*coarse->_rdimensions[d]==fine->_rdimensions[d]);
  }

PARALLEL_FOR_LOOP
  for(int sf=0;sf<fine->oSites();sf++){
    
    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);

    GridBase::CoorFromIndex(coor_f,sf,fine->_rdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    GridBase::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

    // z = A x + y
    fineZ._odata[sf]=coarseA._odata[sc]*fineX._odata[sf]+fineY._odata[sf];

  }

  return;
}
template<class vobj,class CComplex>
  inline void blockInnerProduct(Lattice<CComplex> &CoarseInner,
				const Lattice<vobj> &fineX,
				const Lattice<vobj> &fineY)
{
  typedef decltype(innerProduct(fineX._odata[0],fineY._odata[0])) dotp;

  GridBase *coarse(CoarseInner._grid);
  GridBase *fine  (fineX._grid);

  Lattice<dotp> fine_inner(fine);
  Lattice<dotp> coarse_inner(coarse);

  fine_inner = localInnerProduct(fineX,fineY);
  blockSum(coarse_inner,fine_inner);
PARALLEL_FOR_LOOP
  for(int ss=0;ss<coarse->oSites();ss++){
    CoarseInner._odata[ss] = coarse_inner._odata[ss];
  }
}
template<class vobj,class CComplex>
inline void blockNormalise(Lattice<CComplex> &ip,Lattice<vobj> &fineX)
{
  GridBase *coarse = ip._grid;
  Lattice<vobj> zz(fineX._grid); zz=zero;
  blockInnerProduct(ip,fineX,fineX);
  ip = pow(ip,-0.5);
  blockZAXPY(fineX,ip,fineX,zz);
}
// useful in multigrid project;
// Generic name : Coarsen?
template<class vobj>
inline void blockSum(Lattice<vobj> &coarseData,const Lattice<vobj> &fineData)
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

template<class vobj>
inline void blockPick(GridBase *coarse,const Lattice<vobj> &unpicked,Lattice<vobj> &picked,std::vector<int> coor)
{
  GridBase * fine = unpicked._grid;

  Lattice<vobj> zz(fine);
  Lattice<iScalar<vInteger> > fcoor(fine);

  zz = zero;

  picked = unpicked;
  for(int d=0;d<fine->_ndimension;d++){
    LatticeCoordinate(fcoor,d);
    int block= fine->_rdimensions[d] / coarse->_rdimensions[d];
    int lo   = (coor[d])*block;
    int hi   = (coor[d]+1)*block;
    picked = where( (fcoor<hi) , picked, zz);
    picked = where( (fcoor>=lo), picked, zz);
  }
}

template<class vobj,class CComplex>
inline void blockOrthogonalise(Lattice<CComplex> &ip,std::vector<Lattice<vobj> > &Basis)
{
  GridBase *coarse = ip._grid;
  GridBase *fine   = Basis[0]._grid;

  int       nbasis = Basis.size() ;
  int  _ndimension = coarse->_ndimension;

  // checks
  subdivides(coarse,fine); 
  for(int i=0;i<nbasis;i++){
    conformable(Basis[i]._grid,fine);
  }

  for(int v=0;v<nbasis;v++) {
    for(int u=0;u<v;u++) {
      //Inner product & remove component 
      blockInnerProduct(ip,Basis[u],Basis[v]);
      ip = -ip;
      blockZAXPY<vobj,CComplex> (Basis[v],ip,Basis[u],Basis[v]);
    }
    blockNormalise(ip,Basis[v]);
  }
}

template<class vobj,class CComplex,int nbasis>
inline void blockPromote(const Lattice<iVector<CComplex,nbasis > > &coarseData,
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
    conformable(Basis[i]._grid,fine);
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
      if(i==0) fineData._odata[sf]=coarseData._odata[sc](i) * Basis[i]._odata[sf];
      else     fineData._odata[sf]=fineData._odata[sf]+coarseData._odata[sc](i)*Basis[i]._odata[sf];

    }
  }
  return;
  
}


template<class vobj>
void Replicate(Lattice<vobj> &coarse,Lattice<vobj> & fine)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *cg = coarse._grid;
  GridBase *fg =   fine._grid;

  int nd = cg->_ndimension;

  subdivides(cg,fg); 

  assert(cg->_ndimension==fg->_ndimension);

  std::vector<int> ratio(cg->_ndimension);

  for(int d=0;d<cg->_ndimension;d++){
    ratio[d] = fg->_fdimensions[d]/cg->_fdimensions[d];
  }

  std::vector<int> fcoor(nd);
  std::vector<int> ccoor(nd);
  for(int g=0;g<fg->gSites();g++){

    fg->GlobalIndexToGlobalCoor(g,fcoor);
    for(int d=0;d<nd;d++){
      ccoor[d] = fcoor[d]%cg->_gdimensions[d];
    }
    
    sobj tmp;
    peekSite(tmp,coarse,ccoor);
    pokeSite(tmp,fine,fcoor);
  }

}


}
#endif
