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
template<class vobj> inline void pickCheckerboard(int cb,Lattice<vobj> &half,const Lattice<vobj> &full)
{
  half.Checkerboard() = cb;

  autoView( half_v, half, CpuWrite);
  autoView( full_v, full, CpuRead);
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
template<class vobj> inline void setCheckerboard(Lattice<vobj> &full,const Lattice<vobj> &half)
{
  int cb = half.Checkerboard();
  autoView( half_v , half, CpuRead);
  autoView( full_v , full, CpuWrite);
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

template<class vobj> inline void acceleratorPickCheckerboard(int cb,Lattice<vobj> &half,const Lattice<vobj> &full, int checker_dim_half=0)
{
  half.Checkerboard() = cb;
  autoView(half_v, half, AcceleratorWrite);
  autoView(full_v, full, AcceleratorRead);
  Coordinate rdim_full             = full.Grid()->_rdimensions;
  Coordinate rdim_half             = half.Grid()->_rdimensions;
  unsigned long ndim_half          = half.Grid()->_ndimension;
  Coordinate checker_dim_mask_half = half.Grid()->_checker_dim_mask;
  Coordinate ostride_half          = half.Grid()->_ostride;
  accelerator_for(ss, full.Grid()->oSites(),full.Grid()->Nsimd(),{
    
    Coordinate coor;
    int cbos;
    int linear=0;

    Lexicographic::CoorFromIndex(coor,ss,rdim_full);
    assert(coor.size()==ndim_half);

    for(int d=0;d<ndim_half;d++){ 
      if(checker_dim_mask_half[d]) linear += coor[d];
    }
    cbos = (linear&0x1);

    if (cbos==cb) {
      int ssh=0;
      for(int d=0;d<ndim_half;d++) {
        if (d == checker_dim_half) ssh += ostride_half[d] * ((coor[d] / 2) % rdim_half[d]);
        else ssh += ostride_half[d] * (coor[d] % rdim_half[d]);
      }
      coalescedWrite(half_v[ssh],full_v(ss));
    }
  });
}
template<class vobj> inline void acceleratorSetCheckerboard(Lattice<vobj> &full,const Lattice<vobj> &half, int checker_dim_half=0)
{
  int cb = half.Checkerboard();
  autoView(half_v , half, AcceleratorRead);
  autoView(full_v , full, AcceleratorWrite);
  Coordinate rdim_full             = full.Grid()->_rdimensions;
  Coordinate rdim_half             = half.Grid()->_rdimensions;
  unsigned long ndim_half          = half.Grid()->_ndimension;
  Coordinate checker_dim_mask_half = half.Grid()->_checker_dim_mask;
  Coordinate ostride_half          = half.Grid()->_ostride;
  accelerator_for(ss,full.Grid()->oSites(),full.Grid()->Nsimd(),{

    Coordinate coor;
    int cbos;
    int linear=0;
  
    Lexicographic::CoorFromIndex(coor,ss,rdim_full);
    assert(coor.size()==ndim_half);

    for(int d=0;d<ndim_half;d++){ 
      if(checker_dim_mask_half[d]) linear += coor[d];
    }
    cbos = (linear&0x1);

    if (cbos==cb) {
      int ssh=0;
      for(int d=0;d<ndim_half;d++){
        if (d == checker_dim_half) ssh += ostride_half[d] * ((coor[d] / 2) % rdim_half[d]);
        else ssh += ostride_half[d] * (coor[d] % rdim_half[d]);
      }
      coalescedWrite(full_v[ss],half_v(ssh));
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

template<typename T>
accelerator_inline EnableIf<isGridFundamental<T>> convertType(T & out, const T & in) {
  out = in;
}

// This would allow for conversions between GridFundamental types, but is not strictly needed as yet
/*template<typename T1, typename T2>
accelerator_inline typename std::enable_if<isGridFundamental<T1>::value && isGridFundamental<T2>::value>::type
// Or to make this very broad, conversions between anything that's not a GridTensor could be allowed
//accelerator_inline typename std::enable_if<!isGridTensor<T1>::value && !isGridTensor<T2>::value>::type
convertType(T1 & out, const T2 & in) {
  out = in;
}*/

#ifdef GRID_SIMT
accelerator_inline void convertType(vComplexF & out, const ComplexF & in) {
  ((ComplexF*)&out)[acceleratorSIMTlane(vComplexF::Nsimd())] = in;
}
accelerator_inline void convertType(vComplexD & out, const ComplexD & in) {
  ((ComplexD*)&out)[acceleratorSIMTlane(vComplexD::Nsimd())] = in;
}
accelerator_inline void convertType(vComplexD2 & out, const ComplexD & in) {
  ((ComplexD*)&out)[acceleratorSIMTlane(vComplexD::Nsimd()*2)] = in;
}
#endif

accelerator_inline void convertType(vComplexF & out, const vComplexD2 & in) {
  precisionChange(out,in);
}

accelerator_inline void convertType(vComplexD2 & out, const vComplexF & in) {
  precisionChange(out,in);
}

template<typename T1,typename T2>
accelerator_inline void convertType(iScalar<T1> & out, const iScalar<T2> & in) {
  convertType(out._internal,in._internal);
}

template<typename T1,typename T2>
accelerator_inline NotEnableIf<isGridScalar<T1>> convertType(T1 & out, const iScalar<T2> & in) {
  convertType(out,in._internal);
}

template<typename T1,typename T2>
accelerator_inline NotEnableIf<isGridScalar<T2>> convertType(iScalar<T1> & out, const T2 & in) {
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

template<typename T1,typename T2>
accelerator_inline void convertType(Lattice<T1> & out, const Lattice<T2> & in) {
  autoView( out_v , out,AcceleratorWrite);
  autoView( in_v  , in ,AcceleratorRead);
  accelerator_for(ss,out_v.size(),T1::Nsimd(),{
      convertType(out_v[ss],in_v(ss));
  });
}

////////////////////////////////////////////////////////////////////////////////////////////
// precision-promoted local inner product
////////////////////////////////////////////////////////////////////////////////////////////
template<class vobj>
inline auto localInnerProductD(const Lattice<vobj> &lhs,const Lattice<vobj> &rhs)
-> Lattice<iScalar<decltype(TensorRemove(innerProductD2(lhs.View(CpuRead)[0],rhs.View(CpuRead)[0])))>>
{
  autoView( lhs_v , lhs, AcceleratorRead);
  autoView( rhs_v , rhs, AcceleratorRead);

  typedef decltype(TensorRemove(innerProductD2(lhs_v[0],rhs_v[0]))) t_inner;
  Lattice<iScalar<t_inner>> ret(lhs.Grid());

  {
    autoView(ret_v, ret,AcceleratorWrite);
    accelerator_for(ss,rhs_v.size(),vobj::Nsimd(),{
      convertType(ret_v[ss],innerProductD2(lhs_v(ss),rhs_v(ss)));
    });
  }
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

  autoView( coarseData_ , coarseData, AcceleratorWrite);
  autoView( ip_         , ip,         AcceleratorWrite);
  RealD t_IP=0;
  RealD t_co=0;
  RealD t_za=0;
  for(int v=0;v<nbasis;v++) {
    t_IP-=usecond();
    blockInnerProductD(ip,Basis[v],fineDataRed); // ip = <basis|fine>
    t_IP+=usecond();
    t_co-=usecond();
    accelerator_for( sc, coarse->oSites(), vobj::Nsimd(), {
	convertType(coarseData_[sc](v),ip_[sc]);
    });
    t_co+=usecond();

    // improve numerical stability of projection
    // |fine> = |fine> - <basis|fine> |basis>
    ip=-ip;
    t_za-=usecond();
    blockZAXPY(fineDataRed,ip,Basis[v],fineDataRed); 
    t_za+=usecond();
  }
  //  std::cout << GridLogPerformance << " blockProject : blockInnerProduct :  "<<t_IP<<" us"<<std::endl;
  //  std::cout << GridLogPerformance << " blockProject : conv              :  "<<t_co<<" us"<<std::endl;
  //  std::cout << GridLogPerformance << " blockProject : blockZaxpy        :  "<<t_za<<" us"<<std::endl;
}
// This only minimises data motion from CPU to GPU
// there is chance of better implementation that does a vxk loop of inner products to data share
// at the GPU thread level
template<class vobj,class CComplex,int nbasis,class VLattice>
inline void batchBlockProject(std::vector<Lattice<iVector<CComplex,nbasis>>> &coarseData,
                               const std::vector<Lattice<vobj>> &fineData,
                               const VLattice &Basis)
{
  int NBatch = fineData.size();
  assert(coarseData.size() == NBatch);

  GridBase * fine  = fineData[0].Grid();
  GridBase * coarse= coarseData[0].Grid();

  Lattice<iScalar<CComplex>> ip(coarse);
  std::vector<Lattice<vobj>> fineDataCopy = fineData;

  autoView(ip_, ip, AcceleratorWrite);
  for(int v=0;v<nbasis;v++) {
    for (int k=0; k<NBatch; k++) {
      autoView( coarseData_ , coarseData[k], AcceleratorWrite);
      blockInnerProductD(ip,Basis[v],fineDataCopy[k]); // ip = <basis|fine>
      accelerator_for( sc, coarse->oSites(), vobj::Nsimd(), {
        convertType(coarseData_[sc](v),ip_[sc]);
      });

      // improve numerical stability of projection
      // |fine> = |fine> - <basis|fine> |basis>
      ip=-ip;
      blockZAXPY(fineDataCopy[k],ip,Basis[v],fineDataCopy[k]); 
    }
  }
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

  autoView( fineZ_  , fineZ, AcceleratorWrite);
  autoView( fineX_  , fineX, AcceleratorRead);
  autoView( fineY_  , fineY, AcceleratorRead);
  autoView( coarseA_, coarseA, AcceleratorRead);
  Coordinate fine_rdimensions = fine->_rdimensions;
  Coordinate coarse_rdimensions = coarse->_rdimensions;

  accelerator_for(sf, fine->oSites(), CComplex::Nsimd(), {

      int sc;
      Coordinate coor_c(_ndimension);
      Coordinate coor_f(_ndimension);

      Lexicographic::CoorFromIndex(coor_f,sf,fine_rdimensions);
      for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
      Lexicographic::IndexFromCoor(coor_c,sc,coarse_rdimensions);

      // z = A x + y
#ifdef GRID_SIMT
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

  // Precision promotion
  RealD t;
  t=-usecond();
  fine_inner = localInnerProductD<vobj>(fineX,fineY);
  //  t+=usecond(); std::cout << GridLogPerformance << " blockInnerProduct : localInnerProductD "<<t<<" us"<<std::endl;
  
  t=-usecond();
  blockSum(coarse_inner,fine_inner);
  //  t+=usecond(); std::cout << GridLogPerformance << " blockInnerProduct : blockSum "<<t<<" us"<<std::endl;
  t=-usecond();
  {
    autoView( CoarseInner_  , CoarseInner,AcceleratorWrite);
    autoView( coarse_inner_ , coarse_inner,AcceleratorRead);
    accelerator_for(ss, coarse->oSites(), 1, {
      convertType(CoarseInner_[ss], TensorRemove(coarse_inner_[ss]));
    });
  }
  //  t+=usecond(); std::cout << GridLogPerformance << " blockInnerProduct : convertType "<<t<<" us"<<std::endl;
 
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
  fine_inner = localInnerProduct(fineX,fineY);
  blockSum(coarse_inner,fine_inner);
  {
    autoView( CoarseInner_  , CoarseInner, AcceleratorWrite);
    autoView( coarse_inner_ , coarse_inner, AcceleratorRead);
    accelerator_for(ss, coarse->oSites(), 1, {
	CoarseInner_[ss] = coarse_inner_[ss];
    });
  }
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
  const int maxsubsec=256;
  typedef iVector<vobj,maxsubsec> vSubsec;

  GridBase * fine  = fineData.Grid();
  GridBase * coarse= coarseData.Grid();

  subdivides(coarse,fine); // require they map

  int _ndimension = coarse->_ndimension;

  Coordinate  block_r      (_ndimension);

  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
  }
  int blockVol = fine->oSites()/coarse->oSites();

  // Turn this around to loop threaded over sc and interior loop 
  // over sf would thread better
  autoView( coarseData_ , coarseData, AcceleratorWrite);
  autoView( fineData_   , fineData, AcceleratorRead);

  auto coarseData_p  = &coarseData_[0];
  auto fineData_p    = &fineData_[0];
  
  Coordinate fine_rdimensions = fine->_rdimensions;
  Coordinate coarse_rdimensions = coarse->_rdimensions;

  vobj zz = Zero();

  // Somewhat lazy calculation
  // Find the biggest power of two subsection divisor less than or equal to maxsubsec
  int subsec=maxsubsec;
  int subvol;
  subvol=blockVol/subsec;
  while(subvol*subsec!=blockVol){
    subsec = subsec/2;
    subvol=blockVol/subsec;
  };

  Lattice<vSubsec> coarseTmp(coarse);
  autoView( coarseTmp_, coarseTmp, AcceleratorWriteDiscard);
  auto coarseTmp_p= &coarseTmp_[0];
  
  // Sum within subsecs in a first kernel
  accelerator_for(sce,subsec*coarse->oSites(),vobj::Nsimd(),{

      int sc=sce/subsec;
      int e=sce%subsec;
      
      // One thread per sub block
      Coordinate coor_c(_ndimension);
      Lexicographic::CoorFromIndex(coor_c,sc,coarse_rdimensions);  // Block coordinate

      auto cd = coalescedRead(zz);
      for(int sb=e*subvol;sb<MIN((e+1)*subvol,blockVol);sb++){
	int sf;
	Coordinate coor_b(_ndimension);
	Coordinate coor_f(_ndimension);
	Lexicographic::CoorFromIndex(coor_b,sb,block_r);               // Block sub coordinate
	for(int d=0;d<_ndimension;d++) coor_f[d]=coor_c[d]*block_r[d] + coor_b[d];
	Lexicographic::IndexFromCoor(coor_f,sf,fine_rdimensions);
	
	cd=cd+coalescedRead(fineData_p[sf]);
      }

      coalescedWrite(coarseTmp_[sc](e),cd);

    });
   // Sum across subsecs in a second kernel
   accelerator_for(sc,coarse->oSites(),vobj::Nsimd(),{
      auto cd = coalescedRead(coarseTmp_p[sc](0));
      for(int e=1;e<subsec;e++){
	cd=cd+coalescedRead(coarseTmp_p[sc](e));
      }
      coalescedWrite(coarseData_p[sc],cd);
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

#ifdef GRID_ACCELERATED
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
  autoView( fineData_   , fineData, AcceleratorWrite);
  autoView( coarseData_ , coarseData, AcceleratorRead);

  typedef LatticeView<vobj> Vview;
  std::vector<Vview> AcceleratorVecViewContainer_h; 
  for(int v=0;v<nbasis;v++) {
    AcceleratorVecViewContainer_h.push_back(Basis[v].View(AcceleratorRead));
  }
  static deviceVector<Vview> AcceleratorVecViewContainer; AcceleratorVecViewContainer.resize(nbasis); 
  acceleratorCopyToDevice(&AcceleratorVecViewContainer_h[0],&AcceleratorVecViewContainer[0],nbasis *sizeof(Vview));
  auto Basis_p = &AcceleratorVecViewContainer[0];
  // Loop with a cache friendly loop ordering
  Coordinate frdimensions=fine->_rdimensions;
  Coordinate crdimensions=coarse->_rdimensions;
  accelerator_for(sf,fine->oSites(),vobj::Nsimd(),{
    int sc;
    Coordinate coor_c(_ndimension);
    Coordinate coor_f(_ndimension);

    Lexicographic::CoorFromIndex(coor_f,sf,frdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    Lexicographic::IndexFromCoor(coor_c,sc,crdimensions);

    auto sum= coarseData_(sc)(0) *Basis_p[0](sf);
    for(int i=1;i<nbasis;i++) sum = sum + coarseData_(sc)(i)*Basis_p[i](sf);
    coalescedWrite(fineData_[sf],sum);
  });
  for(int v=0;v<nbasis;v++) {
    AcceleratorVecViewContainer_h[v].ViewClose();
  }
  return;
}
#else
// CPU version
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

    //Lattice<CComplex> cip(coarse);
    //autoView( cip_ , cip, AcceleratorWrite);
    //autoView(  ip_ ,  ip, AcceleratorRead);
    //accelerator_forNB(sc,coarse->oSites(),CComplex::Nsimd(),{
    //	coalescedWrite(cip_[sc], ip_(sc)());
    //  });
    //blockZAXPY<vobj,CComplex >(fineData,cip,Basis[i],fineData);
    blockZAXPY(fineData,ip,Basis[i],fineData);
  }
}
#endif

template<class vobj,class CComplex,int nbasis,class VLattice>
inline void batchBlockPromote(const std::vector<Lattice<iVector<CComplex,nbasis>>> &coarseData,
                               std::vector<Lattice<vobj>> &fineData,
                               const VLattice &Basis)
{
  int NBatch = coarseData.size();
  assert(fineData.size() == NBatch);

  GridBase * fine   = fineData[0].Grid();
  GridBase * coarse = coarseData[0].Grid();
  for (int k=0; k<NBatch; k++)
    fineData[k]=Zero();
  for (int i=0;i<nbasis;i++) {
    for (int k=0; k<NBatch; k++) {
      Lattice<iScalar<CComplex>> ip = PeekIndex<0>(coarseData[k],i);
      blockZAXPY(fineData[k],ip,Basis[i],fineData[k]);
    }
  }
}

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

  autoView(in_v,in,CpuRead);
  autoView(out_v,out,CpuWrite);
  thread_for(idx, ig->lSites(),{
    sobj s;
    ssobj ss;

    Coordinate lcoor(ni);
    ig->LocalIndexToLocalCoor(idx,lcoor);
    peekLocalSite(s,in_v,lcoor);
    ss=s;
    pokeLocalSite(ss,out_v,lcoor);
  });
}

template<class vobj>
void localCopyRegion(const Lattice<vobj> &From,Lattice<vobj> & To,Coordinate FromLowerLeft, Coordinate ToLowerLeft, Coordinate RegionSize)
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  const int words=sizeof(vobj)/sizeof(vector_type);

  //////////////////////////////////////////////////////////////////////////////////////////
  // checks should guarantee that the operations are local
  //////////////////////////////////////////////////////////////////////////////////////////

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

  ///////////////////////////////////////////////////////////
  // do the index calc on the GPU
  ///////////////////////////////////////////////////////////
  Coordinate f_ostride = Fg->_ostride;
  Coordinate f_istride = Fg->_istride;
  Coordinate f_rdimensions = Fg->_rdimensions;
  Coordinate t_ostride = Tg->_ostride;
  Coordinate t_istride = Tg->_istride;
  Coordinate t_rdimensions = Tg->_rdimensions;

  size_t nsite = 1;
  for(int i=0;i<nd;i++) nsite *= RegionSize[i];

  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  autoView(from_v,From,AcceleratorRead);
  autoView(to_v,To,AcceleratorWrite);

  accelerator_for(idx,nsite,1,{

      Coordinate from_coor, to_coor, base;
      Lexicographic::CoorFromIndex(base,idx,RegionSize);
      for(int i=0;i<nd;i++){
	from_coor[i] = base[i] + FromLowerLeft[i];
	to_coor[i] = base[i] + ToLowerLeft[i];
      }
      int from_oidx = 0; for(int d=0;d<nd;d++) from_oidx+=f_ostride[d]*(from_coor[d]%f_rdimensions[d]);
      int from_lane = 0; for(int d=0;d<nd;d++) from_lane+=f_istride[d]*(from_coor[d]/f_rdimensions[d]);
      int to_oidx   = 0; for(int d=0;d<nd;d++) to_oidx+=t_ostride[d]*(to_coor[d]%t_rdimensions[d]);
      int to_lane   = 0; for(int d=0;d<nd;d++) to_lane+=t_istride[d]*(to_coor[d]/t_rdimensions[d]);

      const vector_type* from = (const vector_type *)&from_v[from_oidx];
      vector_type* to = (vector_type *)&to_v[to_oidx];
      
      scalar_type stmp;
      for(int w=0;w<words;w++){
	stmp = getlane(from[w], from_lane);
	putlane(to[w], stmp, to_lane);
      }
  });
}

template<class vobj>
void InsertSliceFast(const Lattice<vobj> &From,Lattice<vobj> & To,int slice, int orthog)
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  const int words=sizeof(vobj)/sizeof(vector_type);

  //////////////////////////////////////////////////////////////////////////////////////////
  // checks should guarantee that the operations are local
  //////////////////////////////////////////////////////////////////////////////////////////
  GridBase *Fg = From.Grid();
  GridBase *Tg = To.Grid();
  assert(!Fg->_isCheckerBoarded);
  assert(!Tg->_isCheckerBoarded);
  int Nsimd = Fg->Nsimd();
  int nF = Fg->_ndimension;
  int nT = Tg->_ndimension;
  assert(nF+1 == nT);

  ///////////////////////////////////////////////////////////
  // do the index calc on the GPU
  ///////////////////////////////////////////////////////////
  Coordinate f_ostride = Fg->_ostride;
  Coordinate f_istride = Fg->_istride;
  Coordinate f_rdimensions = Fg->_rdimensions;
  Coordinate t_ostride = Tg->_ostride;
  Coordinate t_istride = Tg->_istride;
  Coordinate t_rdimensions = Tg->_rdimensions;
  Coordinate RegionSize = Fg->_ldimensions;
  size_t nsite = 1;
  for(int i=0;i<nF;i++) nsite *= RegionSize[i]; // whole volume of lower dim grid

  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  autoView(from_v,From,AcceleratorRead);
  autoView(to_v,To,AcceleratorWrite);

  accelerator_for(idx,nsite,1,{

      Coordinate from_coor(nF), to_coor(nT);
      Lexicographic::CoorFromIndex(from_coor,idx,RegionSize);
      int j=0;
      for(int i=0;i<nT;i++){
	if ( i!=orthog ) { 
	  to_coor[i] = from_coor[j];
	  j++;
	} else {
	  to_coor[i] = slice;
	}
      }
      int from_oidx = 0; for(int d=0;d<nF;d++) from_oidx+=f_ostride[d]*(from_coor[d]%f_rdimensions[d]);
      int from_lane = 0; for(int d=0;d<nF;d++) from_lane+=f_istride[d]*(from_coor[d]/f_rdimensions[d]);
      int to_oidx   = 0; for(int d=0;d<nT;d++) to_oidx+=t_ostride[d]*(to_coor[d]%t_rdimensions[d]);
      int to_lane   = 0; for(int d=0;d<nT;d++) to_lane+=t_istride[d]*(to_coor[d]/t_rdimensions[d]);

      const vector_type* from = (const vector_type *)&from_v[from_oidx];
      vector_type* to = (vector_type *)&to_v[to_oidx];
      
      scalar_type stmp;
      for(int w=0;w<words;w++){
	stmp = getlane(from[w], from_lane);
	putlane(to[w], stmp, to_lane);
      }
  });
}

template<class vobj>
void ExtractSliceFast(Lattice<vobj> &To,const Lattice<vobj> & From,int slice, int orthog)
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  const int words=sizeof(vobj)/sizeof(vector_type);

  //////////////////////////////////////////////////////////////////////////////////////////
  // checks should guarantee that the operations are local
  //////////////////////////////////////////////////////////////////////////////////////////
  GridBase *Fg = From.Grid();
  GridBase *Tg = To.Grid();
  assert(!Fg->_isCheckerBoarded);
  assert(!Tg->_isCheckerBoarded);
  int Nsimd = Fg->Nsimd();
  int nF = Fg->_ndimension;
  int nT = Tg->_ndimension;
  assert(nT+1 == nF);

  ///////////////////////////////////////////////////////////
  // do the index calc on the GPU
  ///////////////////////////////////////////////////////////
  Coordinate f_ostride = Fg->_ostride;
  Coordinate f_istride = Fg->_istride;
  Coordinate f_rdimensions = Fg->_rdimensions;
  Coordinate t_ostride = Tg->_ostride;
  Coordinate t_istride = Tg->_istride;
  Coordinate t_rdimensions = Tg->_rdimensions;
  Coordinate RegionSize = Tg->_ldimensions;
  size_t nsite = 1;
  for(int i=0;i<nT;i++) nsite *= RegionSize[i]; // whole volume of lower dim grid

  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;

  autoView(from_v,From,AcceleratorRead);
  autoView(to_v,To,AcceleratorWrite);

  accelerator_for(idx,nsite,1,{

      Coordinate from_coor(nF), to_coor(nT);
      Lexicographic::CoorFromIndex(to_coor,idx,RegionSize);
      int j=0;
      for(int i=0;i<nF;i++){
	if ( i!=orthog ) { 
	  from_coor[i] = to_coor[j];
	  j++;
	} else {
	  from_coor[i] = slice;
	}
      }
      int from_oidx = 0; for(int d=0;d<nF;d++) from_oidx+=f_ostride[d]*(from_coor[d]%f_rdimensions[d]);
      int from_lane = 0; for(int d=0;d<nF;d++) from_lane+=f_istride[d]*(from_coor[d]/f_rdimensions[d]);
      int to_oidx   = 0; for(int d=0;d<nT;d++) to_oidx+=t_ostride[d]*(to_coor[d]%t_rdimensions[d]);
      int to_lane   = 0; for(int d=0;d<nT;d++) to_lane+=t_istride[d]*(to_coor[d]/t_rdimensions[d]);

      const vector_type* from = (const vector_type *)&from_v[from_oidx];
      vector_type* to = (vector_type *)&to_v[to_oidx];
      
      scalar_type stmp;
      for(int w=0;w<words;w++){
	stmp = getlane(from[w], from_lane);
	putlane(to[w], stmp, to_lane);
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
  autoView(lowDimv,lowDim,CpuRead);
  autoView(higherDimv,higherDim,CpuWrite);
  thread_for(idx,lg->lSites(),{
    sobj s;
    Coordinate lcoor(nl);
    Coordinate hcoor(nh);
    lg->LocalIndexToLocalCoor(idx,lcoor);
    int ddl=0;
    hcoor[orthog] = slice;
    for(int d=0;d<nh;d++){
      if ( d!=orthog ) { 
	hcoor[d]=lcoor[ddl];
	if ( hg->_checker_dim == d ) {
	  hcoor[d]=hcoor[d]*2; // factor in the full coor for peekLocalSite
	  lcoor[ddl]=lcoor[ddl]*2; // factor in the full coor for peekLocalSite
	}
	ddl++;
      }
      
    }
    peekLocalSite(s,lowDimv,lcoor);
    pokeLocalSite(s,higherDimv,hcoor);
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
  lowDim.Checkerboard() = higherDim.Checkerboard();

  int dl; dl = 0;
  for(int d=0;d<nh;d++){
    if ( d != orthog) {
      assert(lg->_processors[dl]  == hg->_processors[d]);
      assert(lg->_ldimensions[dl] == hg->_ldimensions[d]);
      dl++;
    }
  }
  // the above should guarantee that the operations are local
  autoView(lowDimv,lowDim,CpuWrite);
  autoView(higherDimv,higherDim,CpuRead);
  thread_for(idx,lg->lSites(),{
    sobj s;
    Coordinate lcoor(nl);
    Coordinate hcoor(nh);
    lg->LocalIndexToLocalCoor(idx,lcoor);
    hcoor[orthog] = slice;
    int ddl=0;
    for(int d=0;d<nh;d++){
      if ( d!=orthog ) { 
	hcoor[d]=lcoor[ddl];
	if ( hg->_checker_dim == d ) {
	  hcoor[d]=hcoor[d]*2;     // factor in the full gridd coor for peekLocalSite
	  lcoor[ddl]=lcoor[ddl]*2; // factor in the full coor for peekLocalSite
	}
	ddl++;
      }
    }
    peekLocalSite(s,higherDimv,hcoor);
    pokeLocalSite(s,lowDimv,lcoor);
  });

}

//Can I implement with local copyregion??
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
  Coordinate sz = lg->_ldimensions;
  sz[orthog]=1;
  Coordinate f_ll(nl,0); f_ll[orthog]=slice_lo;
  Coordinate t_ll(nh,0); t_ll[orthog]=slice_hi;
  localCopyRegion(lowDim,higherDim,f_ll,t_ll,sz);
}


template<class vobj>
void ExtractSliceLocal(Lattice<vobj> &lowDim,const Lattice<vobj> & higherDim,int slice_lo,int slice_hi, int orthog)
{
  InsertSliceLocal(higherDim,lowDim,slice_hi,slice_lo,orthog);
}


template<class vobj>
void Replicate(const Lattice<vobj> &coarse,Lattice<vobj> & fine)
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
  for(int64_t g=0;g<fg->gSites();g++){

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
  autoView( in_v  , in, CpuRead);
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
  autoView( out_v , out, CpuWrite);
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

//Very fast precision change. Requires in/out objects to reside on same Grid (e.g. by using double2 for the double-precision field)
template<class VobjOut, class VobjIn>
void precisionChangeFast(Lattice<VobjOut> &out, const Lattice<VobjIn> &in)
{
  typedef typename VobjOut::vector_type Vout;
  typedef typename VobjIn::vector_type Vin;
  const int N = sizeof(VobjOut)/sizeof(Vout);
  conformable(out.Grid(),in.Grid());
  out.Checkerboard() = in.Checkerboard();
  int nsimd = out.Grid()->Nsimd();
  autoView( out_v  , out, AcceleratorWrite);
  autoView(  in_v ,   in, AcceleratorRead);
  accelerator_for(idx,out.Grid()->oSites(),1,{
      Vout *vout = (Vout *)&out_v[idx];
      Vin  *vin  = (Vin  *)&in_v[idx];
      precisionChange(vout,vin,N);
  });
}
//Convert a Lattice from one precision to another (original, slow implementation)
template<class VobjOut, class VobjIn>
void precisionChangeOrig(Lattice<VobjOut> &out, const Lattice<VobjIn> &in)
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
  int in_nsimd = in_grid->Nsimd();
  std::vector<Coordinate > out_icoor(out_nsimd);
      
  for(int lane=0; lane < out_nsimd; lane++){
    out_icoor[lane].resize(ndim);
    out_grid->iCoorFromIindex(out_icoor[lane], lane);
  }
        
  std::vector<SobjOut> in_slex_conv(in_grid->lSites());
  unvectorizeToLexOrdArray(in_slex_conv, in);
    
  autoView( out_v , out, CpuWrite);
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

//The workspace for a precision change operation allowing for the reuse of the mapping to save time on subsequent calls
class precisionChangeWorkspace{
  std::pair<Integer,Integer>* fmap_device; //device pointer
  //maintain grids for checking
  GridBase* _out_grid;
  GridBase* _in_grid;
public:
  precisionChangeWorkspace(GridBase *out_grid, GridBase *in_grid): _out_grid(out_grid), _in_grid(in_grid){
    //Build a map between the sites and lanes of the output field and the input field as we cannot use the Grids on the device
    assert(out_grid->Nd() == in_grid->Nd());
    for(int d=0;d<out_grid->Nd();d++){
      assert(out_grid->FullDimensions()[d] == in_grid->FullDimensions()[d]);
    }
    int Nsimd_out = out_grid->Nsimd();

    std::vector<Coordinate> out_icorrs(out_grid->Nsimd()); //reuse these
    for(int lane=0; lane < out_grid->Nsimd(); lane++)
      out_grid->iCoorFromIindex(out_icorrs[lane], lane);
  
    std::vector<std::pair<Integer,Integer> > fmap_host(out_grid->lSites()); //lsites = osites*Nsimd
    thread_for(out_oidx,out_grid->oSites(),{
	Coordinate out_ocorr; 
	out_grid->oCoorFromOindex(out_ocorr, out_oidx);
      
	Coordinate lcorr; //the local coordinate (common to both in and out as full coordinate)
	for(int out_lane=0; out_lane < Nsimd_out; out_lane++){
	  out_grid->InOutCoorToLocalCoor(out_ocorr, out_icorrs[out_lane], lcorr);
	
	  //int in_oidx = in_grid->oIndex(lcorr), in_lane = in_grid->iIndex(lcorr);
	  //Note oIndex and OcorrFromOindex (and same for iIndex) are not inverse for checkerboarded lattice, the former coordinates being defined on the full lattice and the latter on the reduced lattice
	  //Until this is fixed we need to circumvent the problem locally. Here I will use the coordinates defined on the reduced lattice for simplicity
	  int in_oidx = 0, in_lane = 0;
	  for(int d=0;d<in_grid->_ndimension;d++){
	    in_oidx += in_grid->_ostride[d] * ( lcorr[d] % in_grid->_rdimensions[d] );
	    in_lane += in_grid->_istride[d] * ( lcorr[d] / in_grid->_rdimensions[d] );
	  }
	  fmap_host[out_lane + Nsimd_out*out_oidx] = std::pair<Integer,Integer>( in_oidx, in_lane );
	}
      });

    //Copy the map to the device (if we had a way to tell if an accelerator is in use we could avoid this copy for CPU-only machines)
    size_t fmap_bytes = out_grid->lSites() * sizeof(std::pair<Integer,Integer>);
    fmap_device = (std::pair<Integer,Integer>*)acceleratorAllocDevice(fmap_bytes);
    acceleratorCopyToDevice(fmap_host.data(), fmap_device, fmap_bytes); 
  }

  //Prevent moving or copying
  precisionChangeWorkspace(const precisionChangeWorkspace &r) = delete;
  precisionChangeWorkspace(precisionChangeWorkspace &&r) = delete;
  precisionChangeWorkspace &operator=(const precisionChangeWorkspace &r) = delete;
  precisionChangeWorkspace &operator=(precisionChangeWorkspace &&r) = delete;
  
  std::pair<Integer,Integer> const* getMap() const{ return fmap_device; }

  void checkGrids(GridBase* out, GridBase* in) const{
    conformable(out, _out_grid);
    conformable(in, _in_grid);
  }
  
  ~precisionChangeWorkspace(){
    acceleratorFreeDevice(fmap_device);
  }
};


//We would like to use precisionChangeFast when possible. However usage of this requires the Grids to be the same (runtime check)
//*and* the precisionChange(VobjOut::vector_type, VobjIn, int) function to be defined for the types; this requires an extra compile-time check which we do using some SFINAE trickery
template<class VobjOut, class VobjIn>
auto _precisionChangeFastWrap(Lattice<VobjOut> &out, const Lattice<VobjIn> &in, int dummy)->decltype( precisionChange( ((typename VobjOut::vector_type*)0), ((typename VobjIn::vector_type*)0), 1), int()){
  if(out.Grid() == in.Grid()){
    precisionChangeFast(out,in);
    return 1;
  }else{
    return 0;
  }
}
template<class VobjOut, class VobjIn>
int _precisionChangeFastWrap(Lattice<VobjOut> &out, const Lattice<VobjIn> &in, long dummy){ //note long here is intentional; it means the above is preferred if available
  return 0;
}


//Convert a lattice of one precision to another. Much faster than original implementation but requires a pregenerated workspace
//which contains the mapping data.
template<class VobjOut, class VobjIn>
void precisionChange(Lattice<VobjOut> &out, const Lattice<VobjIn> &in, const precisionChangeWorkspace &workspace){
  if(_precisionChangeFastWrap(out,in,0)) return;
  
  static_assert( std::is_same<typename VobjOut::scalar_typeD, typename VobjIn::scalar_typeD>::value == 1, "precisionChange: tensor types must be the same" ); //if tensor types are same the DoublePrecision type must be the same

  out.Checkerboard() = in.Checkerboard();
  constexpr int Nsimd_out = VobjOut::Nsimd();

  workspace.checkGrids(out.Grid(),in.Grid());
  std::pair<Integer,Integer> const* fmap_device = workspace.getMap();

  //Do the copy/precision change
  autoView( out_v , out, AcceleratorWrite);
  autoView( in_v , in, AcceleratorRead);

  accelerator_for(out_oidx, out.Grid()->oSites(), 1,{
      std::pair<Integer,Integer> const* fmap_osite = fmap_device + out_oidx*Nsimd_out;
      for(int out_lane=0; out_lane < Nsimd_out; out_lane++){      
	int in_oidx = fmap_osite[out_lane].first;
	int in_lane = fmap_osite[out_lane].second;
	copyLane(out_v[out_oidx], out_lane, in_v[in_oidx], in_lane);
      }
    });
}

//Convert a Lattice from one precision to another. Much faster than original implementation but slower than precisionChangeFast
//or precisionChange called with pregenerated workspace, as it needs to internally generate the workspace on the host and copy to device
template<class VobjOut, class VobjIn>
void precisionChange(Lattice<VobjOut> &out, const Lattice<VobjIn> &in){
  if(_precisionChangeFastWrap(out,in,0)) return;   
  precisionChangeWorkspace workspace(out.Grid(), in.Grid());
  precisionChange(out, in, workspace);
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

//////////////////////////////////////////////////////
// Faster but less accurate blockProject
//////////////////////////////////////////////////////
template<class vobj,class CComplex,int nbasis,class VLattice>
inline void blockProjectFast(Lattice<iVector<CComplex,nbasis > > &coarseData,
			     const             Lattice<vobj>   &fineData,
			     const VLattice &Basis)
{
  GridBase * fine  = fineData.Grid();
  GridBase * coarse= coarseData.Grid();

  Lattice<iScalar<CComplex> > ip(coarse);

  autoView( coarseData_ , coarseData, AcceleratorWrite);
  autoView( ip_         , ip,         AcceleratorWrite);
  RealD t_IP=0;
  RealD t_co=0;
  for(int v=0;v<nbasis;v++) {
    t_IP-=usecond();
    blockInnerProductD(ip,Basis[v],fineData); 
    t_IP+=usecond();
    t_co-=usecond();
    accelerator_for( sc, coarse->oSites(), vobj::Nsimd(), {
	convertType(coarseData_[sc](v),ip_[sc]);
      });
    t_co+=usecond();
  }
}


NAMESPACE_END(Grid);

