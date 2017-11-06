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
    //parallel_for
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
    //parallel_for
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

  // Loop over coars parallel, and then loop over fine associated with coarse.
  parallel_for(int sf=0;sf<fine->oSites();sf++){

    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);
    Lexicographic::CoorFromIndex(coor_f,sf,fine->_rdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    Lexicographic::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

PARALLEL_CRITICAL
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
  assert(fineX.checkerboard==fineY.checkerboard);
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

  parallel_for(int sf=0;sf<fine->oSites();sf++){
    
    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);

    Lexicographic::CoorFromIndex(coor_f,sf,fine->_rdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    Lexicographic::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

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

  Lattice<dotp> fine_inner(fine); fine_inner.checkerboard = fineX.checkerboard;
  Lattice<dotp> coarse_inner(coarse);

  // Precision promotion?
  fine_inner = localInnerProduct(fineX,fineY);
  blockSum(coarse_inner,fine_inner);
  parallel_for(int ss=0;ss<coarse->oSites();ss++){
    CoarseInner._odata[ss] = coarse_inner._odata[ss];
  }
}
template<class vobj,class CComplex>
inline void blockNormalise(Lattice<CComplex> &ip,Lattice<vobj> &fineX)
{
  GridBase *coarse = ip._grid;
  Lattice<vobj> zz(fineX._grid); zz=zero; zz.checkerboard=fineX.checkerboard;
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

  // Turn this around to loop threaded over sc and interior loop 
  // over sf would thread better
  coarseData=zero;
  parallel_region {

    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);

    parallel_for_internal(int sf=0;sf<fine->oSites();sf++){
    
      Lexicographic::CoorFromIndex(coor_f,sf,fine->_rdimensions);
      for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
      Lexicographic::IndexFromCoor(coor_c,sc,coarse->_rdimensions);
      
PARALLEL_CRITICAL
      coarseData._odata[sc]=coarseData._odata[sc]+fineData._odata[sf];

    }
  }
  return;
}

template<class vobj>
inline void blockPick(GridBase *coarse,const Lattice<vobj> &unpicked,Lattice<vobj> &picked,std::vector<int> coor)
{
  GridBase * fine = unpicked._grid;

  Lattice<vobj> zz(fine); zz.checkerboard = unpicked.checkerboard;
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
  parallel_region {
    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);

    parallel_for_internal(int sf=0;sf<fine->oSites();sf++){

      Lexicographic::CoorFromIndex(coor_f,sf,fine->_rdimensions);
      for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
      Lexicographic::IndexFromCoor(coor_c,sc,coarse->_rdimensions);
      
      for(int i=0;i<nbasis;i++) {
	if(i==0) fineData._odata[sf]=coarseData._odata[sc](i) * Basis[i]._odata[sf];
	else     fineData._odata[sf]=fineData._odata[sf]+coarseData._odata[sc](i)*Basis[i]._odata[sf];
      }
    }
  }
  return;
  
}

// Useful for precision conversion, or indeed anything where an operator= does a conversion on scalars.
// Simd layouts need not match since we use peek/poke Local
template<class vobj,class vvobj>
void localConvert(const Lattice<vobj> &in,Lattice<vvobj> &out)
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vvobj::scalar_object ssobj;

  GridBase *ig = in._grid;
  GridBase *og = out._grid;

  int ni = ig->_ndimension;
  int no = og->_ndimension;

  assert(ni == no);

  for(int d=0;d<no;d++){
    assert(ig->_processors[d]  == og->_processors[d]);
    assert(ig->_ldimensions[d] == og->_ldimensions[d]);
    assert(ig->lSites() == og->lSites());
  }

  parallel_for(int idx=0;idx<ig->lSites();idx++){
    sobj s;
    ssobj ss;

    std::vector<int> lcoor(ni);
    ig->LocalIndexToLocalCoor(idx,lcoor);
    peekLocalSite(s,in,lcoor);
    ss=s;
    pokeLocalSite(ss,out,lcoor);
  }
}


template<class vobj>
void InsertSlice(const Lattice<vobj> &lowDim,Lattice<vobj> & higherDim,int slice, int orthog)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim._grid;
  GridBase *hg = higherDim._grid;
  int nl = lg->_ndimension;
  int nh = hg->_ndimension;

  assert(nl+1 == nh);
  assert(orthog<nh);
  assert(orthog>=0);
  assert(hg->_processors[orthog]==1);

  int dl; dl = 0;
  for(int d=0;d<nh;d++){
    if ( d != orthog) {
      assert(lg->_processors[dl]  == hg->_processors[d]);
      assert(lg->_ldimensions[dl] == hg->_ldimensions[d]);
      dl++;
    }
  }

  // the above should guarantee that the operations are local
  parallel_for(int idx=0;idx<lg->lSites();idx++){
    sobj s;
    std::vector<int> lcoor(nl);
    std::vector<int> hcoor(nh);
    lg->LocalIndexToLocalCoor(idx,lcoor);
    int ddl=0;
    hcoor[orthog] = slice;
    for(int d=0;d<nh;d++){
      if ( d!=orthog ) { 
	hcoor[d]=lcoor[ddl++];
      }
    }
    peekLocalSite(s,lowDim,lcoor);
    pokeLocalSite(s,higherDim,hcoor);
  }
}

template<class vobj>
void ExtractSlice(Lattice<vobj> &lowDim,const Lattice<vobj> & higherDim,int slice, int orthog)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim._grid;
  GridBase *hg = higherDim._grid;
  int nl = lg->_ndimension;
  int nh = hg->_ndimension;

  assert(nl+1 == nh);
  assert(orthog<nh);
  assert(orthog>=0);
  assert(hg->_processors[orthog]==1);

  int dl; dl = 0;
    for(int d=0;d<nh;d++){
      if ( d != orthog) {
	assert(lg->_processors[dl]  == hg->_processors[d]);
	assert(lg->_ldimensions[dl] == hg->_ldimensions[d]);
	dl++;
    }
  }
  // the above should guarantee that the operations are local
  parallel_for(int idx=0;idx<lg->lSites();idx++){
    sobj s;
    std::vector<int> lcoor(nl);
    std::vector<int> hcoor(nh);
    lg->LocalIndexToLocalCoor(idx,lcoor);
    int ddl=0;
    hcoor[orthog] = slice;
    for(int d=0;d<nh;d++){
      if ( d!=orthog ) { 
	hcoor[d]=lcoor[ddl++];
      }
    }
    peekLocalSite(s,higherDim,hcoor);
    pokeLocalSite(s,lowDim,lcoor);
  }

}


template<class vobj>
void InsertSliceLocal(const Lattice<vobj> &lowDim, Lattice<vobj> & higherDim,int slice_lo,int slice_hi, int orthog)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim._grid;
  GridBase *hg = higherDim._grid;
  int nl = lg->_ndimension;
  int nh = hg->_ndimension;

  assert(nl == nh);
  assert(orthog<nh);
  assert(orthog>=0);

  for(int d=0;d<nh;d++){
    assert(lg->_processors[d]  == hg->_processors[d]);
    assert(lg->_ldimensions[d] == hg->_ldimensions[d]);
  }

  // the above should guarantee that the operations are local
  parallel_for(int idx=0;idx<lg->lSites();idx++){
    sobj s;
    std::vector<int> lcoor(nl);
    std::vector<int> hcoor(nh);
    lg->LocalIndexToLocalCoor(idx,lcoor);
    if( lcoor[orthog] == slice_lo ) { 
      hcoor=lcoor;
      hcoor[orthog] = slice_hi;
      peekLocalSite(s,lowDim,lcoor);
      pokeLocalSite(s,higherDim,hcoor);
    }
  }
}


template<class vobj>
void ExtractSliceLocal(Lattice<vobj> &lowDim, Lattice<vobj> & higherDim,int slice_lo,int slice_hi, int orthog)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim._grid;
  GridBase *hg = higherDim._grid;
  int nl = lg->_ndimension;
  int nh = hg->_ndimension;

  assert(nl == nh);
  assert(orthog<nh);
  assert(orthog>=0);

  for(int d=0;d<nh;d++){
    assert(lg->_processors[d]  == hg->_processors[d]);
    assert(lg->_ldimensions[d] == hg->_ldimensions[d]);
  }

  // the above should guarantee that the operations are local
  parallel_for(int idx=0;idx<lg->lSites();idx++){
    sobj s;
    std::vector<int> lcoor(nl);
    std::vector<int> hcoor(nh);
    lg->LocalIndexToLocalCoor(idx,lcoor);
    if( lcoor[orthog] == slice_lo ) { 
      hcoor=lcoor;
      hcoor[orthog] = slice_hi;
      peekLocalSite(s,higherDim,hcoor);
      pokeLocalSite(s,lowDim,lcoor);
    }
  }
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

//Copy SIMD-vectorized lattice to array of scalar objects in lexicographic order
template<typename vobj, typename sobj>
typename std::enable_if<isSIMDvectorized<vobj>::value && !isSIMDvectorized<sobj>::value, void>::type 
unvectorizeToLexOrdArray(std::vector<sobj> &out, const Lattice<vobj> &in)
{

  typedef typename vobj::vector_type vtype;
  
  GridBase* in_grid = in._grid;
  out.resize(in_grid->lSites());
  
  int ndim = in_grid->Nd();
  int in_nsimd = vtype::Nsimd();

  std::vector<std::vector<int> > in_icoor(in_nsimd);
      
  for(int lane=0; lane < in_nsimd; lane++){
    in_icoor[lane].resize(ndim);
    in_grid->iCoorFromIindex(in_icoor[lane], lane);
  }
  
  parallel_for(int in_oidx = 0; in_oidx < in_grid->oSites(); in_oidx++){ //loop over outer index
    //Assemble vector of pointers to output elements
    std::vector<sobj*> out_ptrs(in_nsimd);

    std::vector<int> in_ocoor(ndim);
    in_grid->oCoorFromOindex(in_ocoor, in_oidx);

    std::vector<int> lcoor(in_grid->Nd());
      
    for(int lane=0; lane < in_nsimd; lane++){
      for(int mu=0;mu<ndim;mu++)
	lcoor[mu] = in_ocoor[mu] + in_grid->_rdimensions[mu]*in_icoor[lane][mu];

      int lex;
      Lexicographic::IndexFromCoor(lcoor, lex, in_grid->_ldimensions);
      out_ptrs[lane] = &out[lex];
    }
    
    //Unpack into those ptrs
    const vobj & in_vobj = in._odata[in_oidx];
    extract1(in_vobj, out_ptrs, 0);
  }
}
//Copy SIMD-vectorized lattice to array of scalar objects in lexicographic order
template<typename vobj, typename sobj>
typename std::enable_if<isSIMDvectorized<vobj>::value 
                    && !isSIMDvectorized<sobj>::value, void>::type 
vectorizeFromLexOrdArray( std::vector<sobj> &in, Lattice<vobj> &out)
{

  typedef typename vobj::vector_type vtype;
  
  GridBase* grid = out._grid;
  assert(in.size()==grid->lSites());
  
  int ndim     = grid->Nd();
  int nsimd    = vtype::Nsimd();

  std::vector<std::vector<int> > icoor(nsimd);
      
  for(int lane=0; lane < nsimd; lane++){
    icoor[lane].resize(ndim);
    grid->iCoorFromIindex(icoor[lane],lane);
  }
  
  parallel_for(uint64_t oidx = 0; oidx < grid->oSites(); oidx++){ //loop over outer index
    //Assemble vector of pointers to output elements
    std::vector<sobj*> ptrs(nsimd);

    std::vector<int> ocoor(ndim);
    grid->oCoorFromOindex(ocoor, oidx);

    std::vector<int> lcoor(grid->Nd());
      
    for(int lane=0; lane < nsimd; lane++){

      for(int mu=0;mu<ndim;mu++){
	lcoor[mu] = ocoor[mu] + grid->_rdimensions[mu]*icoor[lane][mu];
      }

      int lex;
      Lexicographic::IndexFromCoor(lcoor, lex, grid->_ldimensions);
      ptrs[lane] = &in[lex];
    }
    
    //pack from those ptrs
    vobj vecobj;
    merge1(vecobj, ptrs, 0);
    out._odata[oidx] = vecobj; 
  }
}

//Convert a Lattice from one precision to another
template<class VobjOut, class VobjIn>
void precisionChange(Lattice<VobjOut> &out, const Lattice<VobjIn> &in){
  assert(out._grid->Nd() == in._grid->Nd());
  out.checkerboard = in.checkerboard;
  GridBase *in_grid=in._grid;
  GridBase *out_grid = out._grid;

  typedef typename VobjOut::scalar_object SobjOut;
  typedef typename VobjIn::scalar_object SobjIn;

  int ndim = out._grid->Nd();
  int out_nsimd = out_grid->Nsimd();
    
  std::vector<std::vector<int> > out_icoor(out_nsimd);
      
  for(int lane=0; lane < out_nsimd; lane++){
    out_icoor[lane].resize(ndim);
    out_grid->iCoorFromIindex(out_icoor[lane], lane);
  }
        
  std::vector<SobjOut> in_slex_conv(in_grid->lSites());
  unvectorizeToLexOrdArray(in_slex_conv, in);
    
  parallel_for(uint64_t out_oidx=0;out_oidx<out_grid->oSites();out_oidx++){
    std::vector<int> out_ocoor(ndim);
    out_grid->oCoorFromOindex(out_ocoor, out_oidx);

    std::vector<SobjOut*> ptrs(out_nsimd);      

    std::vector<int> lcoor(out_grid->Nd());
      
    for(int lane=0; lane < out_nsimd; lane++){
      for(int mu=0;mu<ndim;mu++)
	lcoor[mu] = out_ocoor[mu] + out_grid->_rdimensions[mu]*out_icoor[lane][mu];
	
      int llex; Lexicographic::IndexFromCoor(lcoor, llex, out_grid->_ldimensions);
      ptrs[lane] = &in_slex_conv[llex];
    }
    merge(out._odata[out_oidx], ptrs, 0);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Communicate between grids
////////////////////////////////////////////////////////////////////////////////
//
// All to all plan
//
// Subvolume on fine grid is v.    Vectors a,b,c,d 
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// SIMPLEST CASE:
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mesh of nodes (2) ; subdivide to  1 subdivisions
//
// Lex ord:   
//          N0 va0 vb0  N1 va1 vb1 
//
// For each dimension do an all to all
//
// full AllToAll(0)
//          N0 va0 va1    N1 vb0 vb1
//
// REARRANGE
//          N0 va01       N1 vb01
//
// Must also rearrange data to get into the NEW lex order of grid at each stage. Some kind of "insert/extract".
// NB: Easiest to programme if keep in lex order.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// SIMPLE CASE:
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Mesh of nodes (2x2) ; subdivide to  1x1 subdivisions
//
// Lex ord:   
//          N0 va0 vb0 vc0 vd0       N1 va1 vb1 vc1 vd1  
//          N2 va2 vb2 vc2 vd2       N3 va3 vb3 vc3 vd3 
//
// Ratio = full[dim] / split[dim]
//
// For each dimension do an all to all; get Nvec -> Nvec / ratio
//                                          Ldim -> Ldim * ratio
//                                          LocalVol -> LocalVol * ratio
// full AllToAll(0)
//          N0 va0 vb0 va1 vb1       N1 vc0 vd0 vc1 vd1   
//          N2 va2 vb2 va3 vb3       N3 vc2 vd2 vc3 vd3 
//
// REARRANGE
//          N0 va01 vb01      N1 vc01 vd01
//          N2 va23 vb23      N3 vc23 vd23
//
// full AllToAll(1)           // Not what is wanted. FIXME
//          N0 va01 va23      N1 vc01 vc23 
//          N2 vb01 vb23      N3 vd01 vd23
// 
// REARRANGE
//          N0 va0123      N1 vc0123
//          N2 vb0123      N3 vd0123
//
// Must also rearrange data to get into the NEW lex order of grid at each stage. Some kind of "insert/extract".
// NB: Easiest to programme if keep in lex order.
//
/////////////////////////////////////////////////////////

template<class Vobj>
void Grid_split(std::vector<Lattice<Vobj> > & full,Lattice<Vobj>   & split)
{
  typedef typename Vobj::scalar_object Sobj;

  int full_vecs   = full.size();

  assert(full_vecs>=1);

  GridBase * full_grid = full[0]._grid;
  GridBase *split_grid = split._grid;

  int       ndim  = full_grid->_ndimension;
  int  full_nproc = full_grid->_Nprocessors;
  int split_nproc =split_grid->_Nprocessors;

  ////////////////////////////////
  // Checkerboard management
  ////////////////////////////////
  int cb = full[0].checkerboard;
  split.checkerboard = cb;

  //////////////////////////////
  // Checks
  //////////////////////////////
  assert(full_grid->_ndimension==split_grid->_ndimension);
  for(int n=0;n<full_vecs;n++){
    assert(full[n].checkerboard == cb);
    for(int d=0;d<ndim;d++){
      assert(full[n]._grid->_gdimensions[d]==split._grid->_gdimensions[d]);
      assert(full[n]._grid->_fdimensions[d]==split._grid->_fdimensions[d]);
    }
  }

  int   nvector   =full_nproc/split_nproc; 
  assert(nvector*split_nproc==full_nproc);
  assert(nvector == full_vecs);

  std::vector<int> ratio(ndim);
  for(int d=0;d<ndim;d++){
    ratio[d] = full_grid->_processors[d]/ split_grid->_processors[d];
  }

  uint64_t lsites = full_grid->lSites();
  uint64_t     sz = lsites * nvector;
  std::vector<Sobj> tmpdata(sz);
  std::vector<Sobj> alldata(sz);
  std::vector<Sobj> scalardata(lsites); 

  for(int v=0;v<nvector;v++){
    unvectorizeToLexOrdArray(scalardata,full[v]);    
    parallel_for(int site=0;site<lsites;site++){
      alldata[v*lsites+site] = scalardata[site];
    }
  }

  int nvec = nvector; // Counts down to 1 as we collapse dims
  std::vector<int> ldims = full_grid->_ldimensions;
  std::vector<int> lcoor(ndim);

  for(int d=ndim-1;d>=0;d--){

    if ( ratio[d] != 1 ) {

      full_grid ->AllToAll(d,alldata,tmpdata);
      //      std::cout << GridLogMessage << "Grid_split: dim " <<d<<" ratio "<<ratio[d]<<" nvec "<<nvec<<" procs "<<split_grid->_processors[d]<<std::endl;
      //      for(int v=0;v<nvec;v++){
      //	std::cout << "Grid_split: alldata["<<v<<"] " << alldata[v] <<std::endl;
      //	std::cout << "Grid_split: tmpdata["<<v<<"] " << tmpdata[v] <<std::endl;
      //      }
      //////////////////////////////////////////
      //Local volume for this dimension is expanded by ratio of processor extents
      // Number of vectors is decreased by same factor
      // Rearrange to lexico for bigger volume
      //////////////////////////////////////////
      nvec    /= ratio[d];

      auto rdims = ldims; rdims[d]  *=   ratio[d];
      auto rsites= lsites*ratio[d];
      for(int v=0;v<nvec;v++){

	// For loop over each site within old subvol
	for(int lsite=0;lsite<lsites;lsite++){

	  Lexicographic::CoorFromIndex(lcoor, lsite, ldims);	  

	  for(int r=0;r<ratio[d];r++){ // ratio*nvec terms

	    auto rcoor = lcoor;	    rcoor[d]  += r*ldims[d];

	    int rsite; Lexicographic::IndexFromCoor(rcoor, rsite, rdims);	  
	    rsite += v * rsites;

	    int rmul=nvec*lsites;
	    int vmul=     lsites;
	    alldata[rsite] = tmpdata[lsite+r*rmul+v*vmul];
	    //	    if ( lsite==0 ) {
	    //	      std::cout << "Grid_split: grow alldata["<<rsite<<"] " << alldata[rsite] << " <- tmpdata["<< lsite+r*rmul+v*vmul<<"] "<<tmpdata[lsite+r*rmul+v*vmul]  <<std::endl;
	    //	    }	      
	  }
	}
      }
      ldims[d]*= ratio[d];
      lsites  *= ratio[d];

      if ( split_grid->_processors[d] > 1 ) {
	tmpdata = alldata;
	split_grid->AllToAll(d,tmpdata,alldata);
      }
    }
  }
  vectorizeFromLexOrdArray(alldata,split);    
}

template<class Vobj>
void Grid_split(Lattice<Vobj> &full,Lattice<Vobj>   & split)
{
  int nvector = full._grid->_Nprocessors / split._grid->_Nprocessors;
  std::vector<Lattice<Vobj> > full_v(nvector,full._grid);
  for(int n=0;n<nvector;n++){
    full_v[n] = full;
  }
  Grid_split(full_v,split);
}

template<class Vobj>
void Grid_unsplit(std::vector<Lattice<Vobj> > & full,Lattice<Vobj>   & split)
{
  typedef typename Vobj::scalar_object Sobj;

  int full_vecs   = full.size();

  assert(full_vecs>=1);

  GridBase * full_grid = full[0]._grid;
  GridBase *split_grid = split._grid;

  int       ndim  = full_grid->_ndimension;
  int  full_nproc = full_grid->_Nprocessors;
  int split_nproc =split_grid->_Nprocessors;

  ////////////////////////////////
  // Checkerboard management
  ////////////////////////////////
  int cb = full[0].checkerboard;
  split.checkerboard = cb;

  //////////////////////////////
  // Checks
  //////////////////////////////
  assert(full_grid->_ndimension==split_grid->_ndimension);
  for(int n=0;n<full_vecs;n++){
    assert(full[n].checkerboard == cb);
    for(int d=0;d<ndim;d++){
      assert(full[n]._grid->_gdimensions[d]==split._grid->_gdimensions[d]);
      assert(full[n]._grid->_fdimensions[d]==split._grid->_fdimensions[d]);
    }
  }

  int   nvector   =full_nproc/split_nproc; 
  assert(nvector*split_nproc==full_nproc);
  assert(nvector == full_vecs);

  std::vector<int> ratio(ndim);
  for(int d=0;d<ndim;d++){
    ratio[d] = full_grid->_processors[d]/ split_grid->_processors[d];
  }

  uint64_t lsites = full_grid->lSites();
  uint64_t     sz = lsites * nvector;
  std::vector<Sobj> tmpdata(sz);
  std::vector<Sobj> alldata(sz);
  std::vector<Sobj> scalardata(lsites); 

  unvectorizeToLexOrdArray(alldata,split);    

  /////////////////////////////////////////////////////////////////
  // Start from split grid and work towards full grid
  /////////////////////////////////////////////////////////////////
  std::vector<int> lcoor(ndim);
  std::vector<int> rcoor(ndim);

  int nvec = 1;
  lsites = split_grid->lSites();
  std::vector<int> ldims = split_grid->_ldimensions;

  //  for(int d=ndim-1;d>=0;d--){
  for(int d=0;d<ndim;d++){

    if ( ratio[d] != 1 ) {


      if ( split_grid->_processors[d] > 1 ) {
	tmpdata = alldata;
	split_grid->AllToAll(d,tmpdata,alldata);
      }

      //////////////////////////////////////////
      //Local volume for this dimension is expanded by ratio of processor extents
      // Number of vectors is decreased by same factor
      // Rearrange to lexico for bigger volume
      //////////////////////////////////////////
      auto rsites= lsites/ratio[d];
      auto rdims = ldims; rdims[d]/=ratio[d];

      for(int v=0;v<nvec;v++){

	// rsite, rcoor --> smaller local volume
	// lsite, lcoor --> bigger original (single node?) volume
	// For loop over each site within smaller subvol
	for(int rsite=0;rsite<rsites;rsite++){

	  Lexicographic::CoorFromIndex(rcoor, rsite, rdims);	  
	  int lsite;

	  for(int r=0;r<ratio[d];r++){ 

	    lcoor = rcoor; lcoor[d] += r*rdims[d];
	    Lexicographic::IndexFromCoor(lcoor, lsite, ldims); lsite += v * lsites;

	    int rmul=nvec*rsites;
	    int vmul=     rsites;
	    tmpdata[rsite+r*rmul+v*vmul]=alldata[lsite];

	  }
	}
      }
      nvec   *= ratio[d];
      ldims[d]=rdims[d];
      lsites  =rsites;

      full_grid ->AllToAll(d,tmpdata,alldata);
    }
  }

  lsites = full_grid->lSites();
  for(int v=0;v<nvector;v++){
    assert(v<full.size());
    parallel_for(int site=0;site<lsites;site++){
      scalardata[site] = alldata[v*lsites+site];
    }
    vectorizeFromLexOrdArray(scalardata,full[v]);    
  }
}

 
}
#endif
