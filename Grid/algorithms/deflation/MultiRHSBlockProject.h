/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: MultiRHSDeflation.h

    Copyright (C) 2023

Author: Peter Boyle <pboyle@bnl.gov>

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


/* 
   MultiRHS block projection

   Import basis -> nblock x nbasis x  (block x internal) 
   Import vector of fine lattice objects -> nblock x nrhs x (block x internal) 

   => coarse_(nrhs x nbasis )^block = via batched GEMM

//template<class vobj,class CComplex,int nbasis,class VLattice>
//inline void blockProject(Lattice<iVector<CComplex,nbasis > > &coarseData,
//			   const VLattice &fineData,
//			   const VLattice &Basis)
*/

template<class Field>
class MultiRHSBlockProject
{
public:

  typedef typename Field::scalar_type   scalar;
  typedef typename Field::scalar_object scalar_object;
  typedef Field Fermion;

  int nbasis;
  GridBase *coarse_grid;
  GridBase *fine_grid;
  uint64_t block_vol;
  uint64_t fine_vol;
  uint64_t coarse_vol;
  uint64_t words;

  ////////////////////////////////////////////////////////////////////////////
  // Blocking geometry in full local coordinates. Addressing in lSites rather
  // than (lane,oSite) lets the fine and coarse spaces carry different SIMD
  // layouts, so an unvectorised coarse space may block a vectorised fine one.
  ////////////////////////////////////////////////////////////////////////////
  Coordinate fine_ldimensions;
  Coordinate coarse_ldimensions;
  Coordinate block_ldimensions;
  Coordinate fine_simd;
  Coordinate coarse_simd;

  // Row major layout "C" order:
  // BLAS_V[coarse_vol][nbasis][block_vol][words]
  // BLAS_F[coarse_vol][nrhs][block_vol][words]
  // BLAS_C[coarse_vol][nrhs][nbasis]
  /*
   * in Fortran column major notation (cuBlas order)
   *
   * Vxb = [v1(x)][..][vn(x)] ... x coarse vol
   *
   * Fxr = [r1(x)][..][rm(x)] ... x coarse vol
   *
   * Block project:
   * C_br = V^dag F x coarse vol
   *
   * Block promote:
   * F_xr = Vxb Cbr x coarse_vol
   */  
  deviceVector<scalar> BLAS_V;      // words * block_vol * nbasis x coarse_vol 
  deviceVector<scalar> BLAS_F;      // nrhs x fine_vol * words   -- the sources
  deviceVector<scalar> BLAS_C;      // nrhs x coarse_vol * nbasis -- the coarse coeffs

  RealD blasNorm2(deviceVector<scalar> &blas)
  {
    scalar ss(0.0);
    std::vector<scalar> tmp(blas.size());
    acceleratorCopyFromDevice(&blas[0],&tmp[0],blas.size()*sizeof(scalar));
    for(int64_t s=0;s<blas.size();s++){
      ss=ss+tmp[s]*adj(tmp[s]);
    }
    coarse_grid->GlobalSum(ss);
    return real(ss);
  }
  
  MultiRHSBlockProject(){};
 ~MultiRHSBlockProject(){ Deallocate(); };
  
  void Deallocate(void)
  {
    nbasis=0;
    coarse_grid=nullptr;
    fine_grid=nullptr;
    fine_vol=0;
    block_vol=0;
    coarse_vol=0;
    words=0;
    BLAS_V.resize(0);
    BLAS_F.resize(0);
    BLAS_C.resize(0);
  }
  void Allocate(int _nbasis,GridBase *_fgrid,GridBase *_cgrid)
  {
    nbasis=_nbasis;

    fine_grid=_fgrid;
    coarse_grid=_cgrid;

    fine_vol   = fine_grid->lSites();
    coarse_vol = coarse_grid->lSites();
    block_vol = fine_vol/coarse_vol;

    words = sizeof(scalar_object)/sizeof(scalar);

    int nd = coarse_grid->_ndimension;
    GRID_ASSERT(fine_grid->_ndimension == nd);

    fine_ldimensions.resize(nd);
    coarse_ldimensions.resize(nd);
    block_ldimensions.resize(nd);
    fine_simd   = fine_grid->_simd_layout;
    coarse_simd = coarse_grid->_simd_layout;
    for(int d=0;d<nd;d++){
      GRID_ASSERT(fine_grid->_processors[d] == coarse_grid->_processors[d]);
      fine_ldimensions  [d] = fine_grid->_rdimensions  [d]*fine_grid->_simd_layout  [d];
      coarse_ldimensions[d] = coarse_grid->_rdimensions[d]*coarse_grid->_simd_layout[d];
      block_ldimensions [d] = fine_ldimensions[d]/coarse_ldimensions[d];
      GRID_ASSERT(block_ldimensions[d]*coarse_ldimensions[d] == fine_ldimensions[d]);
    }
    GRID_ASSERT(block_vol == fine_vol/coarse_vol);

    BLAS_V.resize (fine_vol * words * nbasis );
  }
  void ImportFineGridVectors(std::vector <Field > &vecs, deviceVector<scalar> &blas)
  {
    GRID_TRACE("ImportFineGridVectors");
    int nvec = vecs.size();
    typedef typename Field::vector_object vobj;
    //    std::cout << GridLogMessage <<" BlockProjector importing "<<nvec<< " fine grid vectors" <<std::endl;

    GRID_ASSERT(vecs[0].Grid()==fine_grid);

    int _ndimension = coarse_grid->_ndimension;

    uint64_t sz = blas.size();

    acceleratorMemSet(&blas[0],0,blas.size()*sizeof(scalar));

    Coordinate fine_rdimensions = fine_grid->_rdimensions;
    Coordinate coarse_l = coarse_ldimensions;
    Coordinate block_l  = block_ldimensions;
    Coordinate fsimd    = fine_simd;
    int64_t bv= block_vol;
    for(int v=0;v<vecs.size();v++){

      //      std::cout << " BlockProjector importing vector"<<v<<" "<<norm2(vecs[v])<<std::endl;
      autoView( fineData   , vecs[v], AcceleratorRead);

      auto blasData_p  = &blas[0];
      auto fineData_p  = &fineData[0];

      int64_t osites = fine_grid->oSites();

      // loop over fine sites
      const int Nsimd = vobj::Nsimd();
      GRID_ASSERT(sz == coarse_vol * block_vol * nvec * words);
      uint64_t lwords= words; // local variable for copy in to GPU
      accelerator_for(sf,osites,Nsimd,{
#ifdef GRID_SIMT
        {
	  int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
	  for(int lane=0;lane<Nsimd;lane++) {
#endif
	  // One thread per fine site
	  Coordinate coor_f(_ndimension);
	  Coordinate coor_l(_ndimension);
	  Coordinate coor_b(_ndimension);
	  Coordinate coor_c(_ndimension);

	  // Fine (oSite,lane) to full local coor
	  Lexicographic::CoorFromIndex(coor_f,sf,fine_rdimensions);
	  Lexicographic::CoorFromIndex(coor_l,lane,fsimd);
	  for(int d=0;d<_ndimension;d++) coor_f[d] += fine_rdimensions[d]*coor_l[d];

	  for(int d=0;d<_ndimension;d++) coor_b[d] = coor_f[d]%block_l[d];
	  for(int d=0;d<_ndimension;d++) coor_c[d] = coor_f[d]/block_l[d];

	  int sc;// coarse site
	  int sb;// block site
	  Lexicographic::IndexFromCoor(coor_c,sc,coarse_l);
	  Lexicographic::IndexFromCoor(coor_b,sb,block_l);

          scalar_object data = extractLane(lane,fineData[sf]);

	  // BLAS_F[coarse_vol][nvec][block_vol][words]
	  int64_t site = (sc*nvec + v)*bv
	               + sb;

	  //	  GRID_ASSERT(site*lwords<sz);

	  scalar_object * ptr = (scalar_object *)&blasData_p[site*lwords];

	  *ptr = data;
#ifdef GRID_SIMT
	}
#else
	}
#endif
      });
      //      std::cout << " import fine Blas norm "<<blasNorm2(blas)<<std::endl;
      //      std::cout << " BlockProjector imported vector"<<v<<std::endl;
    }
  }
  ////////////////////////////////////////////////////////////////////////////
  // Import direct from multiRHS fine order, avoiding the unpack to a vector
  // of single RHS fields.
  //
  //   fine mrhs grid  : rhs is dimension 0, undistributed, unvectorised
  //   Grid  order     : F[fine_vol][nrhs][words]
  //   BLAS  order     : BLAS_F[lane][coarse_vol][nrhs][block_vol][words]
  //
  // The gather of block_vol from fine_vol, and the transpose of nrhs against
  // block_vol, are the irreducible part: this is not the identity.
  ////////////////////////////////////////////////////////////////////////////
  void ImportFineGridMrhsVectors(Field &vec_mrhs, deviceVector<scalar> &blas)
  {
    typedef typename Field::vector_object vobj;

    GridBase *fine_mrhs_grid = vec_mrhs.Grid();
    int _ndimension = coarse_grid->_ndimension;

    GRID_ASSERT(fine_mrhs_grid->_ndimension == _ndimension+1);
    GRID_ASSERT(fine_mrhs_grid->_simd_layout[0] == 1);
    GRID_ASSERT(fine_mrhs_grid->_processors[0] == 1);
    for(int d=0;d<_ndimension;d++){
      GRID_ASSERT(fine_mrhs_grid->_rdimensions[d+1] == fine_grid->_rdimensions[d]);
      GRID_ASSERT(fine_mrhs_grid->_simd_layout[d+1] == fine_grid->_simd_layout[d]);
    }
    int nvec = fine_mrhs_grid->_rdimensions[0];   // nrhs

    uint64_t sz = blas.size();
    acceleratorMemSet(&blas[0],0,blas.size()*sizeof(scalar));

    Coordinate fine_mrhs_rdimensions = fine_mrhs_grid->_rdimensions;
    Coordinate fine_rdimensions      = fine_grid->_rdimensions;
    Coordinate coarse_l = coarse_ldimensions;
    Coordinate block_l  = block_ldimensions;
    Coordinate fsimd    = fine_simd;
    int64_t bv= block_vol;

    autoView( fineData   , vec_mrhs, AcceleratorRead);
    auto blasData_p  = &blas[0];
    auto fineData_p  = &fineData[0];

    int64_t osites     = fine_grid->oSites();       // D dimensional
    int64_t osites_hi  = fine_mrhs_grid->oSites();  // nvec * osites

    const int Nsimd = vobj::Nsimd();
    GRID_ASSERT(sz == coarse_vol * block_vol * nvec * words);
    uint64_t lwords= words;
    int64_t lnvec = nvec;

    accelerator_for(sfr,osites_hi,Nsimd,{
#ifdef GRID_SIMT
      {
	int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
	for(int lane=0;lane<Nsimd;lane++) {
#endif
	  Coordinate coor_hi(_ndimension+1);
	  Coordinate coor_f(_ndimension);
	  Coordinate coor_l(_ndimension);
	  Coordinate coor_b(_ndimension);
	  Coordinate coor_c(_ndimension);

	  // rhs is dimension 0 of the D+1 grid
	  Lexicographic::CoorFromIndex(coor_hi,sfr,fine_mrhs_rdimensions);
	  int v = coor_hi[0];
	  for(int d=0;d<_ndimension;d++) coor_f[d] = coor_hi[d+1];

	  Lexicographic::CoorFromIndex(coor_l,lane,fsimd);
	  for(int d=0;d<_ndimension;d++) coor_f[d] += fine_rdimensions[d]*coor_l[d];

	  for(int d=0;d<_ndimension;d++) coor_b[d] = coor_f[d]%block_l[d];
	  for(int d=0;d<_ndimension;d++) coor_c[d] = coor_f[d]/block_l[d];

	  int sc;// coarse site
	  int sb;// block site
	  Lexicographic::IndexFromCoor(coor_c,sc,coarse_l);
	  Lexicographic::IndexFromCoor(coor_b,sb,block_l);

	  scalar_object data = extractLane(lane,fineData[sfr]);

	  int64_t site = (sc*lnvec + v)*bv
	               + sb;

	  scalar_object * ptr = (scalar_object *)&blasData_p[site*lwords];
	  *ptr = data;
#ifdef GRID_SIMT
      }
#else
      }
#endif
    });
  }

  ////////////////////////////////////////////////////////////////////////////
  // Export direct to multiRHS coarse order.
  //
  //   Grid  order : C[coarse_vol][nrhs][nbasis]
  //   BLAS  order : BLAS_C[lane][coarse_vol][nrhs][nbasis]
  //
  // At Nsimd()==1 these are the same sequence of addresses.
  ////////////////////////////////////////////////////////////////////////////
  template<class vobj>
  void ExportCoarseGridMrhsVectors(Lattice<vobj> &vec_mrhs, deviceVector<scalar> &blas)
  {
    typedef typename vobj::scalar_object coarse_scalar_object;

    GridBase *coarse_mrhs_grid = vec_mrhs.Grid();
    int _ndimension = coarse_grid->_ndimension;

    GRID_ASSERT(coarse_mrhs_grid->_ndimension == _ndimension+1);
    GRID_ASSERT(coarse_mrhs_grid->_simd_layout[0] == 1);
    GRID_ASSERT(coarse_mrhs_grid->_processors[0] == 1);
    for(int d=0;d<_ndimension;d++){
      GRID_ASSERT(coarse_mrhs_grid->_rdimensions[d+1] == coarse_grid->_rdimensions[d]);
      GRID_ASSERT(coarse_mrhs_grid->_simd_layout[d+1] == coarse_grid->_simd_layout[d]);
    }
    int nvec = coarse_mrhs_grid->_rdimensions[0];   // nrhs

    Coordinate coarse_mrhs_rdimensions = coarse_mrhs_grid->_rdimensions;
    Coordinate coarse_rdimensions      = coarse_grid->_rdimensions;
    Coordinate coarse_l = coarse_ldimensions;
    Coordinate csimd    = coarse_simd;

    autoView( coarseData   , vec_mrhs, AcceleratorWrite);
    auto blasData_p    = &blas[0];
    auto coarseData_p  = &coarseData[0];

    int64_t osites    = coarse_grid->oSites();       // D dimensional
    int64_t osites_hi = coarse_mrhs_grid->oSites();  // nvec * osites

    const int Nsimd = vobj::Nsimd();
    uint64_t cwords=sizeof(typename vobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(cwords==nbasis);
    int64_t lnvec = nvec;

    accelerator_for(scr,osites_hi,Nsimd,{
#ifdef GRID_SIMT
      {
	int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
	for(int lane=0;lane<Nsimd;lane++) {
#endif
	  Coordinate coor_hi(_ndimension+1);
	  Coordinate coor_l(_ndimension);
	  Coordinate coor_c(_ndimension);

	  Lexicographic::CoorFromIndex(coor_hi,scr,coarse_mrhs_rdimensions);
	  int v = coor_hi[0];
	  for(int d=0;d<_ndimension;d++) coor_c[d] = coor_hi[d+1];

	  Lexicographic::CoorFromIndex(coor_l,lane,csimd);
	  for(int d=0;d<_ndimension;d++) coor_c[d] += coarse_rdimensions[d]*coor_l[d];

	  int sc;
	  Lexicographic::IndexFromCoor(coor_c,sc,coarse_l);

	  int64_t blas_site = (sc*lnvec + v)*cwords;
	  coarse_scalar_object * ptr = (coarse_scalar_object *)&blasData_p[blas_site];
	  coarse_scalar_object data = *ptr;
	  insertLane(lane,coarseData[scr],data);
#ifdef GRID_SIMT
      }
#else
      }
#endif
    });
  }

  ////////////////////////////////////////////////////////////////////////////
  // Reverse directions: BLAS_F -> fine mrhs field, coarse mrhs field -> BLAS_C
  ////////////////////////////////////////////////////////////////////////////
  void ExportFineGridMrhsVectors(Field &vec_mrhs, deviceVector<scalar> &blas)
  {
    typedef typename Field::vector_object vobj;

    GridBase *fine_mrhs_grid = vec_mrhs.Grid();
    int _ndimension = coarse_grid->_ndimension;

    GRID_ASSERT(fine_mrhs_grid->_ndimension == _ndimension+1);
    GRID_ASSERT(fine_mrhs_grid->_simd_layout[0] == 1);
    for(int d=0;d<_ndimension;d++){
      GRID_ASSERT(fine_mrhs_grid->_rdimensions[d+1] == fine_grid->_rdimensions[d]);
      GRID_ASSERT(fine_mrhs_grid->_simd_layout[d+1] == fine_grid->_simd_layout[d]);
    }
    int nvec = fine_mrhs_grid->_rdimensions[0];

    Coordinate fine_mrhs_rdimensions = fine_mrhs_grid->_rdimensions;
    Coordinate fine_rdimensions      = fine_grid->_rdimensions;
    Coordinate coarse_l = coarse_ldimensions;
    Coordinate block_l  = block_ldimensions;
    Coordinate fsimd    = fine_simd;
    int64_t bv= block_vol;

    autoView( fineData   , vec_mrhs, AcceleratorWrite);
    auto blasData_p  = &blas[0];
    auto fineData_p  = &fineData[0];

    int64_t osites     = fine_grid->oSites();
    int64_t osites_hi  = fine_mrhs_grid->oSites();

    const int Nsimd = vobj::Nsimd();
    uint64_t lwords= words;
    int64_t lnvec = nvec;

    accelerator_for(sfr,osites_hi,Nsimd,{
#ifdef GRID_SIMT
      {
	int lane=acceleratorSIMTlane(Nsimd);
#else
	for(int lane=0;lane<Nsimd;lane++) {
#endif
	  Coordinate coor_hi(_ndimension+1);
	  Coordinate coor_f(_ndimension);
	  Coordinate coor_l(_ndimension);
	  Coordinate coor_b(_ndimension);
	  Coordinate coor_c(_ndimension);

	  Lexicographic::CoorFromIndex(coor_hi,sfr,fine_mrhs_rdimensions);
	  int v = coor_hi[0];
	  for(int d=0;d<_ndimension;d++) coor_f[d] = coor_hi[d+1];

	  Lexicographic::CoorFromIndex(coor_l,lane,fsimd);
	  for(int d=0;d<_ndimension;d++) coor_f[d] += fine_rdimensions[d]*coor_l[d];

	  for(int d=0;d<_ndimension;d++) coor_b[d] = coor_f[d]%block_l[d];
	  for(int d=0;d<_ndimension;d++) coor_c[d] = coor_f[d]/block_l[d];

	  int sc,sb;
	  Lexicographic::IndexFromCoor(coor_c,sc,coarse_l);
	  Lexicographic::IndexFromCoor(coor_b,sb,block_l);

	  int64_t site = (sc*lnvec + v)*bv
	               + sb;

	  scalar_object * ptr = (scalar_object *)&blasData_p[site*lwords];
	  scalar_object data = *ptr;
	  insertLane(lane,fineData[sfr],data);
#ifdef GRID_SIMT
      }
#else
      }
#endif
    });
  }

  template<class vobj>
  void ImportCoarseGridMrhsVectors(Lattice<vobj> &vec_mrhs, deviceVector<scalar> &blas)
  {
    typedef typename vobj::scalar_object coarse_scalar_object;

    GridBase *coarse_mrhs_grid = vec_mrhs.Grid();
    int _ndimension = coarse_grid->_ndimension;

    GRID_ASSERT(coarse_mrhs_grid->_ndimension == _ndimension+1);
    GRID_ASSERT(coarse_mrhs_grid->_simd_layout[0] == 1);
    for(int d=0;d<_ndimension;d++){
      GRID_ASSERT(coarse_mrhs_grid->_rdimensions[d+1] == coarse_grid->_rdimensions[d]);
      GRID_ASSERT(coarse_mrhs_grid->_simd_layout[d+1] == coarse_grid->_simd_layout[d]);
    }
    int nvec = coarse_mrhs_grid->_rdimensions[0];

    Coordinate coarse_mrhs_rdimensions = coarse_mrhs_grid->_rdimensions;
    Coordinate coarse_rdimensions      = coarse_grid->_rdimensions;
    Coordinate coarse_l = coarse_ldimensions;
    Coordinate csimd    = coarse_simd;

    autoView( coarseData   , vec_mrhs, AcceleratorRead);
    auto blasData_p    = &blas[0];
    auto coarseData_p  = &coarseData[0];

    int64_t osites    = coarse_grid->oSites();
    int64_t osites_hi = coarse_mrhs_grid->oSites();

    const int Nsimd = vobj::Nsimd();
    uint64_t cwords=sizeof(typename vobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(cwords==nbasis);
    int64_t lnvec = nvec;

    accelerator_for(scr,osites_hi,Nsimd,{
#ifdef GRID_SIMT
      {
	int lane=acceleratorSIMTlane(Nsimd);
#else
	for(int lane=0;lane<Nsimd;lane++) {
#endif
	  Coordinate coor_hi(_ndimension+1);
	  Coordinate coor_l(_ndimension);
	  Coordinate coor_c(_ndimension);

	  Lexicographic::CoorFromIndex(coor_hi,scr,coarse_mrhs_rdimensions);
	  int v = coor_hi[0];
	  for(int d=0;d<_ndimension;d++) coor_c[d] = coor_hi[d+1];

	  Lexicographic::CoorFromIndex(coor_l,lane,csimd);
	  for(int d=0;d<_ndimension;d++) coor_c[d] += coarse_rdimensions[d]*coor_l[d];

	  int sc;
	  Lexicographic::IndexFromCoor(coor_c,sc,coarse_l);

	  coarse_scalar_object data = extractLane(lane,coarseData[scr]);

	  int64_t blas_site = (sc*lnvec + v)*cwords;
	  coarse_scalar_object * ptr = (coarse_scalar_object *)&blasData_p[blas_site];
	  *ptr = data;
#ifdef GRID_SIMT
      }
#else
      }
#endif
    });
  }

  void ExportFineGridVectors(std::vector <Field> &vecs, deviceVector<scalar> &blas)
  {
    GRID_TRACE("ExportFineGridVectors");
    typedef typename Field::vector_object vobj;

    int nvec = vecs.size();

    GRID_ASSERT(vecs[0].Grid()==fine_grid);

    int _ndimension = coarse_grid->_ndimension;

    Coordinate fine_rdimensions = fine_grid->_rdimensions;
    Coordinate coarse_l = coarse_ldimensions;
    Coordinate block_l  = block_ldimensions;
    Coordinate fsimd    = fine_simd;

    //    std::cout << " export fine Blas norm "<<blasNorm2(blas)<<std::endl;

    int64_t bv= block_vol;
    for(int v=0;v<vecs.size();v++){

      autoView( fineData   , vecs[v], AcceleratorWrite);

      auto blasData_p  = &blas[0];
      auto fineData_p    = &fineData[0];

      int64_t osites = fine_grid->oSites();
      uint64_t lwords = words;
      //      std::cout << " Nsimd is "<<vobj::Nsimd() << std::endl;
      //      std::cout << " lwords is "<<lwords << std::endl;
      //      std::cout << " sizeof(scalar_object) is "<<sizeof(scalar_object) << std::endl;
      // loop over fine sites
      accelerator_for(sf,osites,vobj::Nsimd(),{
      
#ifdef GRID_SIMT
        {
	  int lane=acceleratorSIMTlane(vobj::Nsimd()); // buffer lane
#else
	  for(int lane=0;lane<vobj::Nsimd();lane++) {
#endif
	  // One thread per fine site
	  Coordinate coor_f(_ndimension);
	  Coordinate coor_l(_ndimension);
	  Coordinate coor_b(_ndimension);
	  Coordinate coor_c(_ndimension);

	  Lexicographic::CoorFromIndex(coor_f,sf,fine_rdimensions);
	  Lexicographic::CoorFromIndex(coor_l,lane,fsimd);
	  for(int d=0;d<_ndimension;d++) coor_f[d] += fine_rdimensions[d]*coor_l[d];

	  for(int d=0;d<_ndimension;d++) coor_b[d] = coor_f[d]%block_l[d];
	  for(int d=0;d<_ndimension;d++) coor_c[d] = coor_f[d]/block_l[d];

	  int sc;
	  int sb;
	  Lexicographic::IndexFromCoor(coor_c,sc,coarse_l);
	  Lexicographic::IndexFromCoor(coor_b,sb,block_l);

	  // BLAS_F[coarse_vol][nvec][block_vol][words]
	  int64_t site = (sc*nvec + v)*bv
	               + sb;

	  scalar_object * ptr = (scalar_object *)&blasData_p[site*lwords];

	  scalar_object data = *ptr;

	  insertLane(lane,fineData[sf],data);
#ifdef GRID_SIMT
	}
#else
	}
#endif
      });
    }
  }
  template<class vobj>
  void ImportCoarseGridVectors(std::vector <Lattice<vobj> > &vecs, deviceVector<scalar> &blas)
  {
    GRID_TRACE("ImportCoarseGridVectors");
    int nvec = vecs.size();
    typedef typename vobj::scalar_object coarse_scalar_object;

    //    std::cout << " BlockProjector importing "<<nvec<< " coarse grid vectors" <<std::endl;

    GRID_ASSERT(vecs[0].Grid()==coarse_grid);

    int _ndimension = coarse_grid->_ndimension;

    uint64_t sz = blas.size();

    Coordinate coarse_rdimensions = coarse_grid->_rdimensions;
    Coordinate coarse_l = coarse_ldimensions;
    Coordinate csimd    = coarse_simd;

    for(int v=0;v<vecs.size();v++){

      //      std::cout << " BlockProjector importing coarse vector"<<v<<" "<<norm2(vecs[v])<<std::endl;
      autoView( coarseData   , vecs[v], AcceleratorRead);

      auto blasData_p  = &blas[0];
      auto coarseData_p  = &coarseData[0];

      int64_t osites = coarse_grid->oSites();

      // loop over fine sites
      const int Nsimd = vobj::Nsimd();
      uint64_t cwords=sizeof(typename vobj::scalar_object)/sizeof(scalar);
      GRID_ASSERT(cwords==nbasis);
      
      accelerator_for(sc,osites,Nsimd,{
#ifdef GRID_SIMT
        {
	  int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
	  for(int lane=0;lane<Nsimd;lane++) {
#endif
           // C_br per site
	    Coordinate coor_c(_ndimension);
	    Coordinate coor_l(_ndimension);
	    Lexicographic::CoorFromIndex(coor_c,sc,coarse_rdimensions);
	    Lexicographic::CoorFromIndex(coor_l,lane,csimd);
	    for(int d=0;d<_ndimension;d++) coor_c[d] += coarse_rdimensions[d]*coor_l[d];

	    int scl;
	    Lexicographic::IndexFromCoor(coor_c,scl,coarse_l);

	    int64_t blas_site = (scl*nvec + v)*cwords;

	    coarse_scalar_object data = extractLane(lane,coarseData[sc]);

	    coarse_scalar_object * ptr = (coarse_scalar_object *)&blasData_p[blas_site];

	    *ptr = data;
#ifdef GRID_SIMT
	}
#else
	}
#endif
      });
      //      std::cout << " import coarsee Blas norm "<<blasNorm2(blas)<<std::endl;
    }
  }
  template<class vobj>
  void ExportCoarseGridVectors(std::vector <Lattice<vobj> > &vecs, deviceVector<scalar> &blas)
  {
    GRID_TRACE("ExportCoarseGridVectors");
    int nvec = vecs.size();
    typedef typename vobj::scalar_object coarse_scalar_object;
    //    std::cout << GridLogMessage<<" BlockProjector exporting "<<nvec<< " coarse grid vectors" <<std::endl;

    GRID_ASSERT(vecs[0].Grid()==coarse_grid);

    int _ndimension = coarse_grid->_ndimension;
    
    uint64_t sz = blas.size();

    Coordinate coarse_rdimensions = coarse_grid->_rdimensions;
    Coordinate coarse_l = coarse_ldimensions;
    Coordinate csimd    = coarse_simd;

    //    std::cout << " export coarsee Blas norm "<<blasNorm2(blas)<<std::endl;
    for(int v=0;v<vecs.size();v++){

      //  std::cout << " BlockProjector exporting coarse vector"<<v<<std::endl;
      autoView( coarseData   , vecs[v], AcceleratorWrite);

      auto blasData_p  = &blas[0];
      auto coarseData_p  = &coarseData[0];

      int64_t osites = coarse_grid->oSites();

      // loop over fine sites
      const int Nsimd = vobj::Nsimd();
      uint64_t cwords=sizeof(typename vobj::scalar_object)/sizeof(scalar);
      GRID_ASSERT(cwords==nbasis);
      
      accelerator_for(sc,osites,Nsimd,{
	  // Wrap in a macro "FOR_ALL_LANES(lane,{ ... });
#ifdef GRID_SIMT
        {
	  int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
	  for(int lane=0;lane<Nsimd;lane++) {
#endif
	    Coordinate coor_c(_ndimension);
	    Coordinate coor_l(_ndimension);
	    Lexicographic::CoorFromIndex(coor_c,sc,coarse_rdimensions);
	    Lexicographic::CoorFromIndex(coor_l,lane,csimd);
	    for(int d=0;d<_ndimension;d++) coor_c[d] += coarse_rdimensions[d]*coor_l[d];

	    int scl;
	    Lexicographic::IndexFromCoor(coor_c,scl,coarse_l);

	    int64_t blas_site = (scl*nvec + v)*cwords;
	    coarse_scalar_object * ptr = (coarse_scalar_object *)&blasData_p[blas_site];
	    coarse_scalar_object data = *ptr;
	    insertLane(lane,coarseData[sc],data);
#ifdef GRID_SIMT
	}
#else
	}
#endif
      });
    }
  }
  void ImportBasis(std::vector < Field > &vecs)
  {
    //    std::cout << " BlockProjector Import basis size "<<vecs.size()<<std::endl;
    ImportFineGridVectors(vecs,BLAS_V);
  }

  template<class cobj>
  void blockProject(std::vector<Field> &fine,std::vector< Lattice<cobj> > & coarse)
  {
    GRID_TRACE("BlockProject");
    int nrhs=fine.size();
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    //    std::cout << "blockProject nbasis " <<nbasis<<" " << _nbasis<<std::endl;
    GRID_ASSERT(nbasis==_nbasis);
    
    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    /////////////////////////////////////////////
    // Copy in the multi-rhs sources to same data layout
    /////////////////////////////////////////////
    //    std::cout << "BlockProject import fine"<<std::endl;
    ImportFineGridVectors(fine,BLAS_F);
    
    deviceVector<scalar *> Vd(coarse_vol);
    deviceVector<scalar *> Fd(coarse_vol);
    deviceVector<scalar *> Cd(coarse_vol);

    //    std::cout << "BlockProject pointers"<<std::endl;
    // ONE bulk transfer per table.  acceleratorPut is a *synchronous* 8-byte
    // hipMemcpy, so the elementwise form emitted 3*coarse_vol of them per call:
    // a traced run showed 272k hipMemcpy calls costing 4.5 s of API time to move
    // 0.7 s worth of bytes.  Same fix as BatchedBlas.h's staging rewrite.
    if ( coarse_vol ) {
      std::vector<scalar *> hVd(coarse_vol), hFd(coarse_vol), hCd(coarse_vol);
      for(int c=0;c<coarse_vol;c++){
	hVd[c] = & BLAS_V[c*nbasis*block_vol*words];
	hFd[c] = & BLAS_F[c*nrhs*block_vol*words];
	hCd[c] = & BLAS_C[c*nrhs*nbasis];
      }
      acceleratorCopyToDevice(&hVd[0],&Vd[0],coarse_vol*sizeof(scalar *));
      acceleratorCopyToDevice(&hFd[0],&Fd[0],coarse_vol*sizeof(scalar *));
      acceleratorCopyToDevice(&hCd[0],&Cd[0],coarse_vol*sizeof(scalar *));
    }

    GridBLAS BLAS;

    //    std::cout << "BlockProject BLAS"<<std::endl;
    int64_t vw = block_vol * words;
    /////////////////////////////////////////
    // C_br = V^dag R
    /////////////////////////////////////////
    BLAS.gemmBatched(GridBLAS_OP_C,GridBLAS_OP_N, 
    		     nbasis,nrhs,vw,
		     scalar(1.0),
		     Vd,
		     Fd,
		     scalar(0.0),  // wipe out C
		     Cd);
    BLAS.synchronise();
    //    std::cout << "BlockProject done"<<std::endl;
    ExportCoarseGridVectors(coarse, BLAS_C);
    //    std::cout << "BlockProject done"<<std::endl;

  }

  template<class cobj>
  void blockPromote(std::vector<Field> &fine,std::vector<Lattice<cobj> > & coarse)
  {
    GRID_TRACE("BlockPromote");
    int nrhs=fine.size();
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(nbasis==_nbasis);
    
    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    ImportCoarseGridVectors(coarse, BLAS_C);

    GridBLAS BLAS;

    deviceVector<scalar *> Vd(coarse_vol);
    deviceVector<scalar *> Fd(coarse_vol);
    deviceVector<scalar *> Cd(coarse_vol);

    // ONE bulk transfer per table.  acceleratorPut is a *synchronous* 8-byte
    // hipMemcpy, so the elementwise form emitted 3*coarse_vol of them per call:
    // a traced run showed 272k hipMemcpy calls costing 4.5 s of API time to move
    // 0.7 s worth of bytes.  Same fix as BatchedBlas.h's staging rewrite.
    if ( coarse_vol ) {
      std::vector<scalar *> hVd(coarse_vol), hFd(coarse_vol), hCd(coarse_vol);
      for(int c=0;c<coarse_vol;c++){
	hVd[c] = & BLAS_V[c*nbasis*block_vol*words];
	hFd[c] = & BLAS_F[c*nrhs*block_vol*words];
	hCd[c] = & BLAS_C[c*nrhs*nbasis];
      }
      acceleratorCopyToDevice(&hVd[0],&Vd[0],coarse_vol*sizeof(scalar *));
      acceleratorCopyToDevice(&hFd[0],&Fd[0],coarse_vol*sizeof(scalar *));
      acceleratorCopyToDevice(&hCd[0],&Cd[0],coarse_vol*sizeof(scalar *));
    }

    /////////////////////////////////////////
    // Block promote:
    // F_xr = Vxb Cbr (x coarse_vol)
    /////////////////////////////////////////

    int64_t vw = block_vol * words;
    BLAS.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N, 
    		     vw,nrhs,nbasis,
		     scalar(1.0),
		     Vd,
		     Cd,
		     scalar(0.0),  // wipe out C
		     Fd);
    BLAS.synchronise();
    //    std::cout << " blas call done"<<std::endl;
    
    ExportFineGridVectors(fine, BLAS_F);
    //    std::cout << " exported "<<std::endl;
  }

  ////////////////////////////////////////////////////////////////////////////
  // multiRHS ordered interfaces.  The GEMM is identical; only the import and
  // export differ, so those are the whole of the layout question.
  ////////////////////////////////////////////////////////////////////////////
  template<class cobj>
  void blockProject(Field &fine_mrhs,Lattice<cobj> &coarse_mrhs)
  {
    GRID_TRACE("BlockProjectMrhs");
    int nrhs = fine_mrhs.Grid()->_rdimensions[0];
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(nbasis==_nbasis);
    GRID_ASSERT(coarse_mrhs.Grid()->_rdimensions[0]==nrhs);

    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    ImportFineGridMrhsVectors(fine_mrhs,BLAS_F);
    ProjectBLAS(nrhs);
    ExportCoarseGridMrhsVectors(coarse_mrhs,BLAS_C);
  }

  template<class cobj>
  void blockPromote(Field &fine_mrhs,Lattice<cobj> &coarse_mrhs)
  {
    GRID_TRACE("BlockPromoteMrhs");
    int nrhs = fine_mrhs.Grid()->_rdimensions[0];
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(nbasis==_nbasis);
    GRID_ASSERT(coarse_mrhs.Grid()->_rdimensions[0]==nrhs);

    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    ImportCoarseGridMrhsVectors(coarse_mrhs,BLAS_C);
    PromoteBLAS(nrhs);
    ExportFineGridMrhsVectors(fine_mrhs,BLAS_F);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Mixed orderings. A single RHS fine operator produces a vector of fine
  // fields with no packing; the coarse side is still wanted in mrhs order.
  ////////////////////////////////////////////////////////////////////////////
  template<class cobj>
  void blockProject(std::vector<Field> &fine,Lattice<cobj> &coarse_mrhs)
  {
    GRID_TRACE("BlockProjectMixed");
    int nrhs = fine.size();
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(nbasis==_nbasis);
    GRID_ASSERT(coarse_mrhs.Grid()->_rdimensions[0]==nrhs);

    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    ImportFineGridVectors(fine,BLAS_F);
    ProjectBLAS(nrhs);
    ExportCoarseGridMrhsVectors(coarse_mrhs,BLAS_C);
  }

  template<class cobj>
  void blockProject(Field &fine_mrhs,std::vector< Lattice<cobj> > &coarse)
  {
    GRID_TRACE("BlockProjectMixed");
    int nrhs = fine_mrhs.Grid()->_rdimensions[0];
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(nbasis==_nbasis);
    GRID_ASSERT(coarse.size()==nrhs);

    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    ImportFineGridMrhsVectors(fine_mrhs,BLAS_F);
    ProjectBLAS(nrhs);
    ExportCoarseGridVectors(coarse,BLAS_C);
  }

  template<class cobj>
  void blockPromote(std::vector<Field> &fine,Lattice<cobj> &coarse_mrhs)
  {
    GRID_TRACE("BlockPromoteMixed");
    int nrhs = fine.size();
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(nbasis==_nbasis);
    GRID_ASSERT(coarse_mrhs.Grid()->_rdimensions[0]==nrhs);

    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    ImportCoarseGridMrhsVectors(coarse_mrhs,BLAS_C);
    PromoteBLAS(nrhs);
    ExportFineGridVectors(fine,BLAS_F);
  }

  template<class cobj>
  void blockPromote(Field &fine_mrhs,std::vector< Lattice<cobj> > &coarse)
  {
    GRID_TRACE("BlockPromoteMixed");
    int nrhs = fine_mrhs.Grid()->_rdimensions[0];
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    GRID_ASSERT(nbasis==_nbasis);
    GRID_ASSERT(coarse.size()==nrhs);

    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    ImportCoarseGridVectors(coarse,BLAS_C);
    PromoteBLAS(nrhs);
    ExportFineGridMrhsVectors(fine_mrhs,BLAS_F);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Pointer tables and the batched GEMM, shared by both orderings
  ////////////////////////////////////////////////////////////////////////////
  void BLASPointers(int nrhs,
		    deviceVector<scalar *> &Vd,
		    deviceVector<scalar *> &Fd,
		    deviceVector<scalar *> &Cd)
  {
    // ONE bulk transfer per table.  acceleratorPut is a *synchronous* 8-byte
    // hipMemcpy, so the elementwise form emitted 3*coarse_vol of them per call:
    // a traced run showed 272k hipMemcpy calls costing 4.5 s of API time to move
    // 0.7 s worth of bytes.  Same fix as BatchedBlas.h's staging rewrite.
    if ( coarse_vol ) {
      std::vector<scalar *> hVd(coarse_vol), hFd(coarse_vol), hCd(coarse_vol);
      for(int c=0;c<coarse_vol;c++){
	hVd[c] = & BLAS_V[c*nbasis*block_vol*words];
	hFd[c] = & BLAS_F[c*nrhs*block_vol*words];
	hCd[c] = & BLAS_C[c*nrhs*nbasis];
      }
      acceleratorCopyToDevice(&hVd[0],&Vd[0],coarse_vol*sizeof(scalar *));
      acceleratorCopyToDevice(&hFd[0],&Fd[0],coarse_vol*sizeof(scalar *));
      acceleratorCopyToDevice(&hCd[0],&Cd[0],coarse_vol*sizeof(scalar *));
    }
  }

  // C_br = V^dag F
  void ProjectBLAS(int nrhs)
  {
    GRID_TRACE("ProjectBLAS");
    deviceVector<scalar *> Vd(coarse_vol);
    deviceVector<scalar *> Fd(coarse_vol);
    deviceVector<scalar *> Cd(coarse_vol);
    BLASPointers(nrhs,Vd,Fd,Cd);

    GridBLAS BLAS;
    int64_t vw = block_vol * words;
    BLAS.gemmBatched(GridBLAS_OP_C,GridBLAS_OP_N,
		     nbasis,nrhs,vw,
		     scalar(1.0),
		     Vd,
		     Fd,
		     scalar(0.0),  // wipe out C
		     Cd);
    BLAS.synchronise();
  }

  // F_xr = Vxb Cbr
  void PromoteBLAS(int nrhs)
  {
    GRID_TRACE("PromoteBLAS");
    deviceVector<scalar *> Vd(coarse_vol);
    deviceVector<scalar *> Fd(coarse_vol);
    deviceVector<scalar *> Cd(coarse_vol);
    BLASPointers(nrhs,Vd,Fd,Cd);

    GridBLAS BLAS;
    int64_t vw = block_vol * words;
    BLAS.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,
		     vw,nrhs,nbasis,
		     scalar(1.0),
		     Vd,
		     Cd,
		     scalar(0.0),  // wipe out F
		     Fd);
    BLAS.synchronise();
  }
};

NAMESPACE_END(Grid);
