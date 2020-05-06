/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_transfer.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Christoph Lehner <christoph@lhnr.de>

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
  half.Checkerboard() = cb;

  auto half_v = half.View();
  auto full_v = full.View();
  thread_for(ss, full.Grid()->oSites(),{
    int cbos;
    Coordinate coor;
    full.Grid()->oCoorFromOindex(coor,ss);
    cbos=half.Grid()->CheckerBoard(coor);

    if (cbos==cb) {
      int ssh=half.Grid()->oIndex(coor);
      half_v[ssh] = full_v[ss];
    }
  });
}

template<class vobj> inline void setCheckerboard(Lattice<vobj> &full,const Lattice<vobj> &half){
  int cb = half.Checkerboard();
  auto half_v = half.View();
  auto full_v = full.View();
  thread_for(ss,full.Grid()->oSites(),{

    Coordinate coor;
    int cbos;

    full.Grid()->oCoorFromOindex(coor,ss);
    cbos=half.Grid()->CheckerBoard(coor);
      
    if (cbos==cb) {
      int ssh=half.Grid()->oIndex(coor);
      full_v[ss]=half_v[ssh];
    }
  });
}

////////////////////////////////////////////////////////////////////////////////////////////
// Flexible Type Conversion for internal promotion to double as well as graceful
// treatment of scalar-compatible types
////////////////////////////////////////////////////////////////////////////////////////////
accelerator_inline void convertType(ComplexD & out, const std::complex<double> & in) {
  out = in;
}

accelerator_inline void convertType(ComplexF & out, const std::complex<float> & in) {
  out = in;
}

#ifdef __CUDA_ARCH__
accelerator_inline void convertType(vComplexF & out, const ComplexF & in) {
  ((ComplexF*)&out)[SIMTlane(vComplexF::Nsimd())] = in;
}
accelerator_inline void convertType(vComplexD & out, const ComplexD & in) {
  ((ComplexD*)&out)[SIMTlane(vComplexD::Nsimd())] = in;
}
accelerator_inline void convertType(vComplexD2 & out, const ComplexD & in) {
  ((ComplexD*)&out)[SIMTlane(vComplexD::Nsimd()*2)] = in;
}
#endif

accelerator_inline void convertType(vComplexF & out, const vComplexD2 & in) {
  out.v = Optimization::PrecisionChange::DtoS(in._internal[0].v,in._internal[1].v);
}

accelerator_inline void convertType(vComplexD2 & out, const vComplexF & in) {
  Optimization::PrecisionChange::StoD(in.v,out._internal[0].v,out._internal[1].v);
}

template<typename T1,typename T2,int N>
  accelerator_inline void convertType(iMatrix<T1,N> & out, const iMatrix<T2,N> & in);
template<typename T1,typename T2,int N>
  accelerator_inline void convertType(iVector<T1,N> & out, const iVector<T2,N> & in);

template<typename T1,typename T2, typename std::enable_if<!isGridScalar<T1>::value, T1>::type* = nullptr>
accelerator_inline void convertType(T1 & out, const iScalar<T2> & in) {
  convertType(out,in._internal);
}

template<typename T1,typename T2>
accelerator_inline void convertType(iScalar<T1> & out, const T2 & in) {
  convertType(out._internal,in);
}

template<typename T1,typename T2,int N>
accelerator_inline void convertType(iMatrix<T1,N> & out, const iMatrix<T2,N> & in) {
  for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      convertType(out._internal[i][j],in._internal[i][j]);
}

template<typename T1,typename T2,int N>
accelerator_inline void convertType(iVector<T1,N> & out, const iVector<T2,N> & in) {
  for (int i=0;i<N;i++)
    convertType(out._internal[i],in._internal[i]);
}

template<typename T, typename std::enable_if<isGridFundamental<T>::value, T>::type* = nullptr>
accelerator_inline void convertType(T & out, const T & in) {
  out = in;
}

template<typename T1,typename T2>
accelerator_inline void convertType(Lattice<T1> & out, const Lattice<T2> & in) {
  auto out_v = out.AcceleratorView(ViewWrite);
  auto in_v  = in.AcceleratorView(ViewRead);

  accelerator_for(ss,out_v.size(),T1::Nsimd(),{
      convertType(out_v[ss],in_v(ss));
    });
}

////////////////////////////////////////////////////////////////////////////////////////////
// precision-promoted local inner product
////////////////////////////////////////////////////////////////////////////////////////////
template<class vobj>
inline auto localInnerProductD(const Lattice<vobj> &lhs,const Lattice<vobj> &rhs)
-> Lattice<iScalar<decltype(TensorRemove(innerProductD2(lhs.View()[0],rhs.View()[0])))>>
{
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  auto rhs_v = rhs.AcceleratorView(ViewRead);

  typedef decltype(TensorRemove(innerProductD2(lhs_v[0],rhs_v[0]))) t_inner;
  Lattice<iScalar<t_inner>> ret(lhs.Grid());
  auto ret_v = ret.AcceleratorView(ViewWrite);

  accelerator_for(ss,rhs_v.size(),vobj::Nsimd(),{
      convertType(ret_v[ss],innerProductD2(lhs_v(ss),rhs_v(ss)));
    });

  return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////
// block routines
////////////////////////////////////////////////////////////////////////////////////////////
template<class vobj,class CComplex,int nbasis,class VLattice>
inline void blockProject(Lattice<iVector<CComplex,nbasis > > &coarseData,
			   const             Lattice<vobj>   &fineData,
			   const VLattice &Basis)
{
  GridBase * fine  = fineData.Grid();
  GridBase * coarse= coarseData.Grid();

  Lattice<iScalar<CComplex>> ip(coarse);
  Lattice<vobj>     fineDataRed = fineData;

  //  auto fineData_   = fineData.View();
  auto coarseData_ = coarseData.AcceleratorView(ViewWrite);
  auto ip_         = ip.AcceleratorView(ViewReadWrite);
  for(int v=0;v<nbasis;v++) {
    blockInnerProductD(ip,Basis[v],fineDataRed); // ip = <basis|fine>
    accelerator_for( sc, coarse->oSites(), vobj::Nsimd(), {
	convertType(coarseData_[sc](v),ip_[sc]);
      });

    // improve numerical stability of projection
    // |fine> = |fine> - <basis|fine> |basis>
    ip=-ip;
    blockZAXPY(fineDataRed,ip,Basis[v],fineDataRed); 
  }
}

template<class vobj,class CComplex,int nbasis>
inline void blockProject1(Lattice<iVector<CComplex,nbasis > > &coarseData,
			 const             Lattice<vobj>   &fineData,
			 const std::vector<Lattice<vobj> > &Basis)
{
  typedef iVector<CComplex,nbasis > coarseSiteData;
  coarseSiteData elide;
  typedef decltype(coalescedRead(elide)) ScalarComplex;
  GridBase * fine  = fineData.Grid();
  GridBase * coarse= coarseData.Grid();
  int  _ndimension = coarse->_ndimension;

  // checks
  assert( nbasis == Basis.size() );
  subdivides(coarse,fine); 
  for(int i=0;i<nbasis;i++){
    conformable(Basis[i],fineData);
  }

  Coordinate block_r      (_ndimension);
  
  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
    assert(block_r[d]*coarse->_rdimensions[d] == fine->_rdimensions[d]);
  }
  int blockVol = fine->oSites()/coarse->oSites();

  coarseData=Zero();

  auto fineData_   = fineData.View();
  auto coarseData_ = coarseData.View();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // To make this lock free, loop over coars parallel, and then loop over fine associated with coarse.
  // Otherwise do fine inner product per site, and make the update atomic
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  accelerator_for( sci, nbasis*coarse->oSites(), vobj::Nsimd(), {

    auto sc=sci/nbasis;
    auto i=sci%nbasis;
    auto Basis_      = Basis[i].View();

    Coordinate coor_c(_ndimension);
    Lexicographic::CoorFromIndex(coor_c,sc,coarse->_rdimensions);  // Block coordinate

    int sf;
    decltype(innerProduct(Basis_(sf),fineData_(sf))) reduce=Zero();

    for(int sb=0;sb<blockVol;sb++){

      Coordinate coor_b(_ndimension);
      Coordinate coor_f(_ndimension);

      Lexicographic::CoorFromIndex(coor_b,sb,block_r);
      for(int d=0;d<_ndimension;d++) coor_f[d]=coor_c[d]*block_r[d]+coor_b[d];
      Lexicographic::IndexFromCoor(coor_f,sf,fine->_rdimensions);
      
      reduce=reduce+innerProduct(Basis_(sf),fineData_(sf));
    }
    coalescedWrite(coarseData_[sc](i),reduce);
  });
  return;
}

template<class vobj,class vobj2,class CComplex>
  inline void blockZAXPY(Lattice<vobj> &fineZ,
			 const Lattice<CComplex> &coarseA,
			 const Lattice<vobj2> &fineX,
			 const Lattice<vobj> &fineY)
{
  GridBase * fine  = fineZ.Grid();
  GridBase * coarse= coarseA.Grid();

  fineZ.Checkerboard()=fineX.Checkerboard();
  assert(fineX.Checkerboard()==fineY.Checkerboard());
  subdivides(coarse,fine); // require they map
  conformable(fineX,fineY);
  conformable(fineX,fineZ);

  int _ndimension = coarse->_ndimension;

  Coordinate  block_r      (_ndimension);

  // FIXME merge with subdivide checking routine as this is redundant
  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
    assert(block_r[d]*coarse->_rdimensions[d]==fine->_rdimensions[d]);
  }

  auto fineZ_  = fineZ.AcceleratorView(ViewWrite);
  auto fineX_  = fineX.AcceleratorView(ViewRead);
  auto fineY_  = fineY.AcceleratorView(ViewRead);
  auto coarseA_= coarseA.AcceleratorView(ViewRead);

  accelerator_for(sf, fine->oSites(), CComplex::Nsimd(), {

      int sc;
      Coordinate coor_c(_ndimension);
      Coordinate coor_f(_ndimension);

      Lexicographic::CoorFromIndex(coor_f,sf,fine->_rdimensions);
      for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
      Lexicographic::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

      // z = A x + y
#ifdef __CUDA_ARCH__
      typename vobj2::tensor_reduced::scalar_object cA;
      typename vobj::scalar_object cAx;
#else
      typename vobj2::tensor_reduced cA;
      vobj cAx;
#endif
      convertType(cA,TensorRemove(coarseA_(sc)));
      auto prod = cA*fineX_(sf);
      convertType(cAx,prod);
      coalescedWrite(fineZ_[sf],cAx+fineY_(sf));

    });

  return;
}

template<class vobj,class CComplex>
  inline void blockInnerProductD(Lattice<CComplex> &CoarseInner,
				 const Lattice<vobj> &fineX,
				 const Lattice<vobj> &fineY)
{
  typedef iScalar<decltype(TensorRemove(innerProductD2(vobj(),vobj())))> dotp;

  GridBase *coarse(CoarseInner.Grid());
  GridBase *fine  (fineX.Grid());

  Lattice<dotp> fine_inner(fine); fine_inner.Checkerboard() = fineX.Checkerboard();
  Lattice<dotp> coarse_inner(coarse);

  auto CoarseInner_  = CoarseInner.AcceleratorView(ViewWrite);
  auto coarse_inner_ = coarse_inner.AcceleratorView(ViewReadWrite);

  // Precision promotion
  fine_inner = localInnerProductD(fineX,fineY);
  blockSum(coarse_inner,fine_inner);
  accelerator_for(ss, coarse->oSites(), 1, {
      convertType(CoarseInner_[ss], TensorRemove(coarse_inner_[ss]));
    });
 
}

template<class vobj,class CComplex> // deprecate
inline void blockInnerProduct(Lattice<CComplex> &CoarseInner,
			      const Lattice<vobj> &fineX,
			      const Lattice<vobj> &fineY)
{
  typedef decltype(innerProduct(vobj(),vobj())) dotp;

  GridBase *coarse(CoarseInner.Grid());
  GridBase *fine  (fineX.Grid());

  Lattice<dotp> fine_inner(fine); fine_inner.Checkerboard() = fineX.Checkerboard();
  Lattice<dotp> coarse_inner(coarse);

  // Precision promotion?
  auto CoarseInner_  = CoarseInner.AcceleratorView(ViewWrite);
  auto coarse_inner_ = coarse_inner.AcceleratorView(ViewReadWrite);

  fine_inner = localInnerProduct(fineX,fineY);
  blockSum(coarse_inner,fine_inner);
  accelerator_for(ss, coarse->oSites(), 1, {
    CoarseInner_[ss] = coarse_inner_[ss];
  });
}

template<class vobj,class CComplex>
inline void blockNormalise(Lattice<CComplex> &ip,Lattice<vobj> &fineX)
{
  GridBase *coarse = ip.Grid();
  Lattice<vobj> zz(fineX.Grid()); zz=Zero(); zz.Checkerboard()=fineX.Checkerboard();
  blockInnerProduct(ip,fineX,fineX);
  ip = pow(ip,-0.5);
  blockZAXPY(fineX,ip,fineX,zz);
}
// useful in multigrid project;
// Generic name : Coarsen?
template<class vobj>
inline void blockSum(Lattice<vobj> &coarseData,const Lattice<vobj> &fineData) 
{
  GridBase * fine  = fineData.Grid();
  GridBase * coarse= coarseData.Grid();

  subdivides(coarse,fine); // require they map

  int _ndimension = coarse->_ndimension;

  Coordinate  block_r      (_ndimension);

  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
  }
  int blockVol = fine->oSites()/coarse->oSites();

  auto coarseData_ = coarseData.AcceleratorView(ViewReadWrite);
  auto fineData_   = fineData.AcceleratorView(ViewRead);

  accelerator_for(sc,coarse->oSites(),1,{

      // One thread per sub block
      Coordinate coor_c(_ndimension);
      Lexicographic::CoorFromIndex(coor_c,sc,coarse->_rdimensions);  // Block coordinate
      coarseData_[sc]=Zero();

      for(int sb=0;sb<blockVol;sb++){

	int sf;
	Coordinate coor_b(_ndimension);
	Coordinate coor_f(_ndimension);
	Lexicographic::CoorFromIndex(coor_b,sb,block_r);               // Block sub coordinate
	for(int d=0;d<_ndimension;d++) coor_f[d]=coor_c[d]*block_r[d] + coor_b[d];
	Lexicographic::IndexFromCoor(coor_f,sf,fine->_rdimensions);

	coarseData_[sc]=coarseData_[sc]+fineData_[sf];
      }

    });
  return;
}


template<class vobj>
inline void blockPick(GridBase *coarse,const Lattice<vobj> &unpicked,Lattice<vobj> &picked,Coordinate coor)
{
  GridBase * fine = unpicked.Grid();

  Lattice<vobj> zz(fine); zz.Checkerboard() = unpicked.Checkerboard();
  Lattice<iScalar<vInteger> > fcoor(fine);

  zz = Zero();

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

template<class CComplex,class VLattice>
inline void blockOrthonormalize(Lattice<CComplex> &ip,VLattice &Basis)
{
  GridBase *coarse = ip.Grid();
  GridBase *fine   = Basis[0].Grid();

  int       nbasis = Basis.size() ;

  // checks
  subdivides(coarse,fine);
  for(int i=0;i<nbasis;i++){
    conformable(Basis[i].Grid(),fine);
  }

  for(int v=0;v<nbasis;v++) {
    for(int u=0;u<v;u++) {
      //Inner product & remove component
      blockInnerProductD(ip,Basis[u],Basis[v]);
      ip = -ip;
      blockZAXPY(Basis[v],ip,Basis[u],Basis[v]);
    }
    blockNormalise(ip,Basis[v]);
  }
}

template<class vobj,class CComplex>
inline void blockOrthogonalise(Lattice<CComplex> &ip,std::vector<Lattice<vobj> > &Basis) // deprecated inaccurate naming
{
  blockOrthonormalize(ip,Basis);
}

#if 0
// TODO: CPU optimized version here
template<class vobj,class CComplex,int nbasis>
inline void blockPromote(const Lattice<iVector<CComplex,nbasis > > &coarseData,
			 Lattice<vobj>   &fineData,
			 const std::vector<Lattice<vobj> > &Basis)
{
  GridBase * fine  = fineData.Grid();
  GridBase * coarse= coarseData.Grid();
  int  _ndimension = coarse->_ndimension;

  // checks
  assert( nbasis == Basis.size() );
  subdivides(coarse,fine); 
  for(int i=0;i<nbasis;i++){
    conformable(Basis[i].Grid(),fine);
  }

  Coordinate  block_r      (_ndimension);
  
  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
  }
  auto fineData_   = fineData.View();
  auto coarseData_ = coarseData.View();

  // Loop with a cache friendly loop ordering
  accelerator_for(sf,fine->oSites(),1,{
    int sc;
    Coordinate coor_c(_ndimension);
    Coordinate coor_f(_ndimension);

    Lexicographic::CoorFromIndex(coor_f,sf,fine->_rdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    Lexicographic::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

    for(int i=0;i<nbasis;i++) {
      auto basis_ = Basis[i].View();
      if(i==0) fineData_[sf]=coarseData_[sc](i) *basis_[sf]);
      else     fineData_[sf]=fineData_[sf]+coarseData_[sc](i)*basis_[sf]);
    }
  });
  return;
  
}
#else
template<class vobj,class CComplex,int nbasis,class VLattice>
inline void blockPromote(const Lattice<iVector<CComplex,nbasis > > &coarseData,
			 Lattice<vobj>   &fineData,
			 const VLattice &Basis)
{
  GridBase * fine  = fineData.Grid();
  GridBase * coarse= coarseData.Grid();
  fineData=Zero();
  for(int i=0;i<nbasis;i++) {
    Lattice<iScalar<CComplex> > ip = PeekIndex<0>(coarseData,i);
    auto  ip_ =  ip.AcceleratorView(ViewRead);
    blockZAXPY(fineData,ip,Basis[i],fineData);
  }
}
#endif

// Useful for precision conversion, or indeed anything where an operator= does a conversion on scalars.
// Simd layouts need not match since we use peek/poke Local
template<class vobj,class vvobj>
void localConvert(const Lattice<vobj> &in,Lattice<vvobj> &out)
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vvobj::scalar_object ssobj;

  GridBase *ig = in.Grid();
  GridBase *og = out.Grid();

  int ni = ig->_ndimension;
  int no = og->_ndimension;

  assert(ni == no);

  for(int d=0;d<no;d++){
    assert(ig->_processors[d]  == og->_processors[d]);
    assert(ig->_ldimensions[d] == og->_ldimensions[d]);
    assert(ig->lSites() == og->lSites());
  }

  thread_for(idx, ig->lSites(),{
    sobj s;
    ssobj ss;

    Coordinate lcoor(ni);
    ig->LocalIndexToLocalCoor(idx,lcoor);
    peekLocalSite(s,in,lcoor);
    ss=s;
    pokeLocalSite(ss,out,lcoor);
  });
}

template<class vobj>
void localCopyRegion(const Lattice<vobj> &From,Lattice<vobj> & To,Coordinate FromLowerLeft, Coordinate ToLowerLeft, Coordinate RegionSize)
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  static const int words=sizeof(vobj)/sizeof(vector_type);

  GridBase *Fg = From.Grid();
  GridBase *Tg = To.Grid();
  assert(!Fg->_isCheckerBoarded);
  assert(!Tg->_isCheckerBoarded);
  int Nsimd = Fg->Nsimd();
  int nF = Fg->_ndimension;
  int nT = Tg->_ndimension;
  int nd = nF;
  assert(nF == nT);

  for(int d=0;d<nd;d++){
    assert(Fg->_processors[d]  == Tg->_processors[d]);
  }

  // the above should guarantee that the operations are local
  Coordinate ldf = Fg->_ldimensions;
  Coordinate rdf = Fg->_rdimensions;
  Coordinate isf = Fg->_istride;
  Coordinate osf = Fg->_ostride;
  Coordinate rdt = Tg->_rdimensions;
  Coordinate ist = Tg->_istride;
  Coordinate ost = Tg->_ostride;
  auto t_v = To.AcceleratorView(ViewWrite);
  auto f_v = From.AcceleratorView(ViewRead);
  accelerator_for(idx,Fg->lSites(),1,{
    sobj s;
    Coordinate Fcoor(nd);
    Coordinate Tcoor(nd);
    Lexicographic::CoorFromIndex(Fcoor,idx,ldf);
    int in_region=1;
    for(int d=0;d<nd;d++){
      if ( (Fcoor[d] < FromLowerLeft[d]) || (Fcoor[d]>=FromLowerLeft[d]+RegionSize[d]) ){ 
	in_region=0;
      }
      Tcoor[d] = ToLowerLeft[d]+ Fcoor[d]-FromLowerLeft[d];
    }
    if (in_region) {
      Integer idx_f = 0; for(int d=0;d<nd;d++) idx_f+=isf[d]*(Fcoor[d]/rdf[d]);
      Integer idx_t = 0; for(int d=0;d<nd;d++) idx_t+=ist[d]*(Tcoor[d]/rdt[d]);
      Integer odx_f = 0; for(int d=0;d<nd;d++) odx_f+=osf[d]*(Fcoor[d]%rdf[d]);
      Integer odx_t = 0; for(int d=0;d<nd;d++) odx_t+=ost[d]*(Tcoor[d]%rdt[d]);
      scalar_type * fp = (scalar_type *)&f_v[odx_f];
      scalar_type * tp = (scalar_type *)&t_v[odx_t];
      for(int w=0;w<words;w++){
	tp[idx_t+w*Nsimd] = fp[idx_f+w*Nsimd];  // FIXME IF RRII layout, type pun no worke
      }
      //      peekLocalSite(s,From,Fcoor);
      //      pokeLocalSite(s,To  ,Tcoor);
    }
  });
}


template<class vobj>
void InsertSlice(const Lattice<vobj> &lowDim,Lattice<vobj> & higherDim,int slice, int orthog)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim.Grid();
  GridBase *hg = higherDim.Grid();
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
  thread_for(idx,lg->lSites(),{
    sobj s;
    Coordinate lcoor(nl);
    Coordinate hcoor(nh);
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
  });
}

template<class vobj>
void ExtractSlice(Lattice<vobj> &lowDim,const Lattice<vobj> & higherDim,int slice, int orthog)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim.Grid();
  GridBase *hg = higherDim.Grid();
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
  thread_for(idx,lg->lSites(),{
    sobj s;
    Coordinate lcoor(nl);
    Coordinate hcoor(nh);
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
  });

}


template<class vobj>
void InsertSliceLocal(const Lattice<vobj> &lowDim, Lattice<vobj> & higherDim,int slice_lo,int slice_hi, int orthog)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim.Grid();
  GridBase *hg = higherDim.Grid();
  int nl = lg->_ndimension;
  int nh = hg->_ndimension;

  assert(nl == nh);
  assert(orthog<nh);
  assert(orthog>=0);

  for(int d=0;d<nh;d++){
    if ( d!=orthog ) {
    assert(lg->_processors[d]  == hg->_processors[d]);
    assert(lg->_ldimensions[d] == hg->_ldimensions[d]);
  }
  }

  // the above should guarantee that the operations are local
  thread_for(idx,lg->lSites(),{
    sobj s;
    Coordinate lcoor(nl);
    Coordinate hcoor(nh);
    lg->LocalIndexToLocalCoor(idx,lcoor);
    if( lcoor[orthog] == slice_lo ) { 
      hcoor=lcoor;
      hcoor[orthog] = slice_hi;
      peekLocalSite(s,lowDim,lcoor);
      pokeLocalSite(s,higherDim,hcoor);
    }
  });
}


template<class vobj>
void ExtractSliceLocal(Lattice<vobj> &lowDim,const Lattice<vobj> & higherDim,int slice_lo,int slice_hi, int orthog)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim.Grid();
  GridBase *hg = higherDim.Grid();
  int nl = lg->_ndimension;
  int nh = hg->_ndimension;

  assert(nl == nh);
  assert(orthog<nh);
  assert(orthog>=0);

  for(int d=0;d<nh;d++){
    if ( d!=orthog ) {
    assert(lg->_processors[d]  == hg->_processors[d]);
    assert(lg->_ldimensions[d] == hg->_ldimensions[d]);
  }
  }

  // the above should guarantee that the operations are local
  thread_for(idx,lg->lSites(),{
    sobj s;
    Coordinate lcoor(nl);
    Coordinate hcoor(nh);
    lg->LocalIndexToLocalCoor(idx,lcoor);
    if( lcoor[orthog] == slice_lo ) { 
      hcoor=lcoor;
      hcoor[orthog] = slice_hi;
      peekLocalSite(s,higherDim,hcoor);
      pokeLocalSite(s,lowDim,lcoor);
    }
  });
}


template<class vobj>
void Replicate(Lattice<vobj> &coarse,Lattice<vobj> & fine)
{
  typedef typename vobj::scalar_object sobj;

  GridBase *cg = coarse.Grid();
  GridBase *fg =   fine.Grid();

  int nd = cg->_ndimension;

  subdivides(cg,fg); 

  assert(cg->_ndimension==fg->_ndimension);

  Coordinate ratio(cg->_ndimension);

  for(int d=0;d<cg->_ndimension;d++){
    ratio[d] = fg->_fdimensions[d]/cg->_fdimensions[d];
  }

  Coordinate fcoor(nd);
  Coordinate ccoor(nd);
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
  
  GridBase* in_grid = in.Grid();
  out.resize(in_grid->lSites());
  
  int ndim = in_grid->Nd();
  int in_nsimd = vtype::Nsimd();

  std::vector<Coordinate > in_icoor(in_nsimd);
      
  for(int lane=0; lane < in_nsimd; lane++){
    in_icoor[lane].resize(ndim);
    in_grid->iCoorFromIindex(in_icoor[lane], lane);
  }

  //loop over outer index
  auto in_v  = in.View();
  thread_for(in_oidx,in_grid->oSites(),{
    //Assemble vector of pointers to output elements
    ExtractPointerArray<sobj> out_ptrs(in_nsimd);

    Coordinate in_ocoor(ndim);
    in_grid->oCoorFromOindex(in_ocoor, in_oidx);

    Coordinate lcoor(in_grid->Nd());
      
    for(int lane=0; lane < in_nsimd; lane++){

      for(int mu=0;mu<ndim;mu++){
	lcoor[mu] = in_ocoor[mu] + in_grid->_rdimensions[mu]*in_icoor[lane][mu];
      }

      int lex;
      Lexicographic::IndexFromCoor(lcoor, lex, in_grid->_ldimensions);
      assert(lex < out.size());
      out_ptrs[lane] = &out[lex];
    }
    
    //Unpack into those ptrs
    const vobj & in_vobj = in_v[in_oidx];
    extract(in_vobj, out_ptrs, 0);
  });
}

template<typename vobj, typename sobj>
typename std::enable_if<isSIMDvectorized<vobj>::value && !isSIMDvectorized<sobj>::value, void>::type 
unvectorizeToRevLexOrdArray(std::vector<sobj> &out, const Lattice<vobj> &in)
{

  typedef typename vobj::vector_type vtype;
  
  GridBase* in_grid = in._grid;
  out.resize(in_grid->lSites());
  
  int ndim = in_grid->Nd();
  int in_nsimd = vtype::Nsimd();

  std::vector<Coordinate > in_icoor(in_nsimd);
      
  for(int lane=0; lane < in_nsimd; lane++){
    in_icoor[lane].resize(ndim);
    in_grid->iCoorFromIindex(in_icoor[lane], lane);
  }
  
  thread_for(in_oidx, in_grid->oSites(),{
    //Assemble vector of pointers to output elements
    std::vector<sobj*> out_ptrs(in_nsimd);

    Coordinate in_ocoor(ndim);
    in_grid->oCoorFromOindex(in_ocoor, in_oidx);

    Coordinate lcoor(in_grid->Nd());
      
    for(int lane=0; lane < in_nsimd; lane++){
      for(int mu=0;mu<ndim;mu++)
	lcoor[mu] = in_ocoor[mu] + in_grid->_rdimensions[mu]*in_icoor[lane][mu];

      int lex;
      Lexicographic::IndexFromCoorReversed(lcoor, lex, in_grid->_ldimensions);
      out_ptrs[lane] = &out[lex];
    }
    
    //Unpack into those ptrs
    const vobj & in_vobj = in._odata[in_oidx];
    extract1(in_vobj, out_ptrs, 0);
  });
}

//Copy SIMD-vectorized lattice to array of scalar objects in lexicographic order
template<typename vobj, typename sobj>
typename std::enable_if<isSIMDvectorized<vobj>::value 
			&& !isSIMDvectorized<sobj>::value, void>::type 
vectorizeFromLexOrdArray( std::vector<sobj> &in, Lattice<vobj> &out)
{

  typedef typename vobj::vector_type vtype;
  
  GridBase* grid = out.Grid();
  assert(in.size()==grid->lSites());
  
  const int ndim     = grid->Nd();
  constexpr int nsimd    = vtype::Nsimd();

  std::vector<Coordinate > icoor(nsimd);
      
  for(int lane=0; lane < nsimd; lane++){
    icoor[lane].resize(ndim);
    grid->iCoorFromIindex(icoor[lane],lane);
  }
  auto out_v = out.View();
  thread_for(oidx, grid->oSites(),{
    //Assemble vector of pointers to output elements
    ExtractPointerArray<sobj> ptrs(nsimd);

    Coordinate ocoor(ndim);
    Coordinate lcoor(ndim);
    grid->oCoorFromOindex(ocoor, oidx);
      
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
    merge(vecobj, ptrs, 0);
    out_v[oidx] = vecobj; 
  });
}

template<typename vobj, typename sobj>
typename std::enable_if<isSIMDvectorized<vobj>::value 
                    && !isSIMDvectorized<sobj>::value, void>::type 
vectorizeFromRevLexOrdArray( std::vector<sobj> &in, Lattice<vobj> &out)
{

  typedef typename vobj::vector_type vtype;
  
  GridBase* grid = out._grid;
  assert(in.size()==grid->lSites());
  
  int ndim     = grid->Nd();
  int nsimd    = vtype::Nsimd();

  std::vector<Coordinate > icoor(nsimd);
      
  for(int lane=0; lane < nsimd; lane++){
    icoor[lane].resize(ndim);
    grid->iCoorFromIindex(icoor[lane],lane);
  }
  
  thread_for(oidx, grid->oSites(), {
    //Assemble vector of pointers to output elements
    std::vector<sobj*> ptrs(nsimd);

    Coordinate ocoor(ndim);
    grid->oCoorFromOindex(ocoor, oidx);

    Coordinate lcoor(grid->Nd());
      
    for(int lane=0; lane < nsimd; lane++){

      for(int mu=0;mu<ndim;mu++){
	lcoor[mu] = ocoor[mu] + grid->_rdimensions[mu]*icoor[lane][mu];
      }

      int lex;
      Lexicographic::IndexFromCoorReversed(lcoor, lex, grid->_ldimensions);
      ptrs[lane] = &in[lex];
    }
    
    //pack from those ptrs
    vobj vecobj;
    merge1(vecobj, ptrs, 0);
    out._odata[oidx] = vecobj; 
  });
}

//Convert a Lattice from one precision to another
template<class VobjOut, class VobjIn>
void precisionChange(Lattice<VobjOut> &out, const Lattice<VobjIn> &in)
{
  assert(out.Grid()->Nd() == in.Grid()->Nd());
  for(int d=0;d<out.Grid()->Nd();d++){
    assert(out.Grid()->FullDimensions()[d] == in.Grid()->FullDimensions()[d]);
  }
  out.Checkerboard() = in.Checkerboard();
  GridBase *in_grid=in.Grid();
  GridBase *out_grid = out.Grid();

  typedef typename VobjOut::scalar_object SobjOut;
  typedef typename VobjIn::scalar_object SobjIn;

  int ndim = out.Grid()->Nd();
  int out_nsimd = out_grid->Nsimd();
    
  std::vector<Coordinate > out_icoor(out_nsimd);
      
  for(int lane=0; lane < out_nsimd; lane++){
    out_icoor[lane].resize(ndim);
    out_grid->iCoorFromIindex(out_icoor[lane], lane);
  }
        
  std::vector<SobjOut> in_slex_conv(in_grid->lSites());
  unvectorizeToLexOrdArray(in_slex_conv, in);
    
  auto out_v = out.View();
  thread_for(out_oidx,out_grid->oSites(),{
    Coordinate out_ocoor(ndim);
    out_grid->oCoorFromOindex(out_ocoor, out_oidx);

    ExtractPointerArray<SobjOut> ptrs(out_nsimd);      

    Coordinate lcoor(out_grid->Nd());
      
    for(int lane=0; lane < out_nsimd; lane++){
      for(int mu=0;mu<ndim;mu++)
	lcoor[mu] = out_ocoor[mu] + out_grid->_rdimensions[mu]*out_icoor[lane][mu];
	
      int llex; Lexicographic::IndexFromCoor(lcoor, llex, out_grid->_ldimensions);
      ptrs[lane] = &in_slex_conv[llex];
    }
    merge(out_v[out_oidx], ptrs, 0);
  });
}

////////////////////////////////////////////////////////////////////////////////
// Communicate between grids
////////////////////////////////////////////////////////////////////////////////
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
/*
 *  Let chunk = (fvol*nvec)/sP be size of a chunk.         ( Divide lexico vol * nvec into fP/sP = M chunks )
 *  
 *  2nd A2A (over sP nodes; subdivide the fP into sP chunks of M)
 * 
 *     node 0     1st chunk of node 0M..(1M-1); 2nd chunk of node 0M..(1M-1)..   data chunk x M x sP = fL / sP * M * sP = fL * M growth
 *     node 1     1st chunk of node 1M..(2M-1); 2nd chunk of node 1M..(2M-1)..
 *     node 2     1st chunk of node 2M..(3M-1); 2nd chunk of node 2M..(3M-1)..
 *     node 3     1st chunk of node 3M..(3M-1); 2nd chunk of node 2M..(3M-1)..
 *  etc...
 */
template<class Vobj>
void Grid_split(std::vector<Lattice<Vobj> > & full,Lattice<Vobj>   & split)
{
  typedef typename Vobj::scalar_object Sobj;

  int full_vecs   = full.size();

  assert(full_vecs>=1);

  GridBase * full_grid = full[0].Grid();
  GridBase *split_grid = split.Grid();

  int       ndim  = full_grid->_ndimension;
  int  full_nproc = full_grid->_Nprocessors;
  int split_nproc =split_grid->_Nprocessors;

  ////////////////////////////////
  // Checkerboard management
  ////////////////////////////////
  int cb = full[0].Checkerboard();
  split.Checkerboard() = cb;

  //////////////////////////////
  // Checks
  //////////////////////////////
  assert(full_grid->_ndimension==split_grid->_ndimension);
  for(int n=0;n<full_vecs;n++){
    assert(full[n].Checkerboard() == cb);
    for(int d=0;d<ndim;d++){
      assert(full[n].Grid()->_gdimensions[d]==split.Grid()->_gdimensions[d]);
      assert(full[n].Grid()->_fdimensions[d]==split.Grid()->_fdimensions[d]);
    }
  }

  int   nvector   =full_nproc/split_nproc; 
  assert(nvector*split_nproc==full_nproc);
  assert(nvector == full_vecs);

  Coordinate ratio(ndim);
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
    thread_for(site,lsites,{
      alldata[v*lsites+site] = scalardata[site];
    });
  }

  int nvec = nvector; // Counts down to 1 as we collapse dims
  Coordinate ldims = full_grid->_ldimensions;

  for(int d=ndim-1;d>=0;d--){

    if ( ratio[d] != 1 ) {

      full_grid ->AllToAll(d,alldata,tmpdata);
      if ( split_grid->_processors[d] > 1 ) {
	alldata=tmpdata;
	split_grid->AllToAll(d,alldata,tmpdata);
      }

      auto rdims = ldims; 
      auto     M = ratio[d];
      auto rsites= lsites*M;// increases rsites by M
      nvec      /= M;       // Reduce nvec by subdivision factor
      rdims[d]  *= M;       // increase local dim by same factor

      int sP =   split_grid->_processors[d];
      int fP =    full_grid->_processors[d];

      int fvol   = lsites;
      
      int chunk  = (nvec*fvol)/sP;          assert(chunk*sP == nvec*fvol);

      // Loop over reordered data post A2A
      thread_for(c, chunk, {
	Coordinate coor(ndim);
	for(int m=0;m<M;m++){
	  for(int s=0;s<sP;s++){
	    
	    // addressing; use lexico
	    int lex_r;
	    uint64_t lex_c        = c+chunk*m+chunk*M*s;
	    uint64_t lex_fvol_vec = c+chunk*s;
	    uint64_t lex_fvol     = lex_fvol_vec%fvol;
	    uint64_t lex_vec      = lex_fvol_vec/fvol;

	    // which node sets an adder to the coordinate
	    Lexicographic::CoorFromIndex(coor, lex_fvol, ldims);	  
	    coor[d] += m*ldims[d];
	    Lexicographic::IndexFromCoor(coor, lex_r, rdims);	  
	    lex_r += lex_vec * rsites;

	    // LexicoFind coordinate & vector number within split lattice
	    alldata[lex_r] = tmpdata[lex_c];

	  }
	}
      });
      ldims[d]*= ratio[d];
      lsites  *= ratio[d];

    }
  }
  vectorizeFromLexOrdArray(alldata,split);    
}

template<class Vobj>
void Grid_split(Lattice<Vobj> &full,Lattice<Vobj>   & split)
{
  int nvector = full.Grid()->_Nprocessors / split.Grid()->_Nprocessors;
  std::vector<Lattice<Vobj> > full_v(nvector,full.Grid());
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

  GridBase * full_grid = full[0].Grid();
  GridBase *split_grid = split.Grid();

  int       ndim  = full_grid->_ndimension;
  int  full_nproc = full_grid->_Nprocessors;
  int split_nproc =split_grid->_Nprocessors;

  ////////////////////////////////
  // Checkerboard management
  ////////////////////////////////
  int cb = full[0].Checkerboard();
  split.Checkerboard() = cb;

  //////////////////////////////
  // Checks
  //////////////////////////////
  assert(full_grid->_ndimension==split_grid->_ndimension);
  for(int n=0;n<full_vecs;n++){
    assert(full[n].Checkerboard() == cb);
    for(int d=0;d<ndim;d++){
      assert(full[n].Grid()->_gdimensions[d]==split.Grid()->_gdimensions[d]);
      assert(full[n].Grid()->_fdimensions[d]==split.Grid()->_fdimensions[d]);
    }
  }

  int   nvector   =full_nproc/split_nproc; 
  assert(nvector*split_nproc==full_nproc);
  assert(nvector == full_vecs);

  Coordinate ratio(ndim);
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

  int nvec = 1;
  uint64_t rsites        = split_grid->lSites();
  Coordinate rdims = split_grid->_ldimensions;

  for(int d=0;d<ndim;d++){

    if ( ratio[d] != 1 ) {

      auto     M = ratio[d];

      int sP =   split_grid->_processors[d];
      int fP =    full_grid->_processors[d];
      
      auto ldims = rdims;  ldims[d]  /= M;  // Decrease local dims by same factor
      auto lsites= rsites/M;                // Decreases rsites by M
      
      int fvol   = lsites;
      int chunk  = (nvec*fvol)/sP;          assert(chunk*sP == nvec*fvol);
	
      {
	// Loop over reordered data post A2A
	thread_for(c, chunk,{
	  Coordinate coor(ndim);
	  for(int m=0;m<M;m++){
	    for(int s=0;s<sP;s++){

	      // addressing; use lexico
	      int lex_r;
	      uint64_t lex_c = c+chunk*m+chunk*M*s;
	      uint64_t lex_fvol_vec = c+chunk*s;
	      uint64_t lex_fvol     = lex_fvol_vec%fvol;
	      uint64_t lex_vec      = lex_fvol_vec/fvol;
	      
	      // which node sets an adder to the coordinate
	      Lexicographic::CoorFromIndex(coor, lex_fvol, ldims);	  
	      coor[d] += m*ldims[d];
	      Lexicographic::IndexFromCoor(coor, lex_r, rdims);	  
	      lex_r += lex_vec * rsites;
	      
	      // LexicoFind coordinate & vector number within split lattice
	      tmpdata[lex_c] = alldata[lex_r];
	    }
	  }
        });
      }

      if ( split_grid->_processors[d] > 1 ) {
	split_grid->AllToAll(d,tmpdata,alldata);
	tmpdata=alldata;
      }
      full_grid ->AllToAll(d,tmpdata,alldata);
      rdims[d]/= M;
      rsites  /= M;
      nvec    *= M;       // Increase nvec by subdivision factor
    }
  }

  lsites = full_grid->lSites();
  for(int v=0;v<nvector;v++){
    thread_for(site, lsites,{
      scalardata[site] = alldata[v*lsites+site];
    });
    vectorizeFromLexOrdArray(scalardata,full[v]);    
  }
}

NAMESPACE_END(Grid);

