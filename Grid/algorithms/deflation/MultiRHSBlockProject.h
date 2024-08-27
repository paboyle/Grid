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

    BLAS_V.resize (fine_vol * words * nbasis );
  }
  void ImportFineGridVectors(std::vector <Field > &vecs, deviceVector<scalar> &blas)
  {
    int nvec = vecs.size();
    typedef typename Field::vector_object vobj;
    //    std::cout << GridLogMessage <<" BlockProjector importing "<<nvec<< " fine grid vectors" <<std::endl;

    assert(vecs[0].Grid()==fine_grid);

    subdivides(coarse_grid,fine_grid); // require they map

    int _ndimension = coarse_grid->_ndimension;
    assert(block_vol == fine_grid->oSites() / coarse_grid->oSites());
    
    Coordinate  block_r      (_ndimension);
    for(int d=0 ; d<_ndimension;d++){
      block_r[d] = fine_grid->_rdimensions[d] / coarse_grid->_rdimensions[d];
    }

    uint64_t sz = blas.size();

    acceleratorMemSet(&blas[0],0,blas.size()*sizeof(scalar));

    Coordinate fine_rdimensions = fine_grid->_rdimensions;
    Coordinate coarse_rdimensions = coarse_grid->_rdimensions;
    int64_t bv= block_vol;
    for(int v=0;v<vecs.size();v++){

      //      std::cout << " BlockProjector importing vector"<<v<<" "<<norm2(vecs[v])<<std::endl;
      autoView( fineData   , vecs[v], AcceleratorRead);

      auto blasData_p  = &blas[0];
      auto fineData_p  = &fineData[0];

      int64_t osites = fine_grid->oSites();

      // loop over fine sites
      const int Nsimd = vobj::Nsimd();
      //      std::cout << "sz "<<sz<<std::endl;
      //      std::cout << "prod "<<Nsimd * coarse_grid->oSites() * block_vol * nvec * words<<std::endl;
      assert(sz == Nsimd * coarse_grid->oSites() * block_vol * nvec * words);
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
	  Coordinate coor_b(_ndimension);
	  Coordinate coor_c(_ndimension);

	  // Fine site to fine coor
	  Lexicographic::CoorFromIndex(coor_f,sf,fine_rdimensions);

	  for(int d=0;d<_ndimension;d++) coor_b[d] = coor_f[d]%block_r[d];
	  for(int d=0;d<_ndimension;d++) coor_c[d] = coor_f[d]/block_r[d];
	  
	  int sc;// coarse site
	  int sb;// block site
	  Lexicographic::IndexFromCoor(coor_c,sc,coarse_rdimensions);
	  Lexicographic::IndexFromCoor(coor_b,sb,block_r);

          scalar_object data = extractLane(lane,fineData[sf]);

	  // BLAS layout address calculation
	  // words * block_vol * nbasis x coarse_vol
	  // coarse oSite x block vole x lanes
	  int64_t site = (lane*osites + sc*bv)*nvec
   	               + v*bv
	               + sb;

	  //	  assert(site*lwords<sz);

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
  void ExportFineGridVectors(std::vector <Field> &vecs, deviceVector<scalar> &blas)
  {
    typedef typename Field::vector_object vobj;

    int nvec = vecs.size();

    assert(vecs[0].Grid()==fine_grid);

    subdivides(coarse_grid,fine_grid); // require they map

    int _ndimension = coarse_grid->_ndimension;
    assert(block_vol == fine_grid->oSites() / coarse_grid->oSites());
    
    Coordinate  block_r      (_ndimension);
    for(int d=0 ; d<_ndimension;d++){
      block_r[d] = fine_grid->_rdimensions[d] / coarse_grid->_rdimensions[d];
    }
    Coordinate fine_rdimensions = fine_grid->_rdimensions;
    Coordinate coarse_rdimensions = coarse_grid->_rdimensions;

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
	  Coordinate coor_b(_ndimension);
	  Coordinate coor_c(_ndimension);

	  Lexicographic::CoorFromIndex(coor_f,sf,fine_rdimensions);

	  for(int d=0;d<_ndimension;d++) coor_b[d] = coor_f[d]%block_r[d];
	  for(int d=0;d<_ndimension;d++) coor_c[d] = coor_f[d]/block_r[d];
	  
	  int sc;
	  int sb;
	  Lexicographic::IndexFromCoor(coor_c,sc,coarse_rdimensions);
	  Lexicographic::IndexFromCoor(coor_b,sb,block_r);

	  // BLAS layout address calculation
	  // words * block_vol * nbasis x coarse_vol 	  
	  int64_t site = (lane*osites + sc*bv)*nvec
   	               + v*bv
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
    int nvec = vecs.size();
    typedef typename vobj::scalar_object coarse_scalar_object;

    //    std::cout << " BlockProjector importing "<<nvec<< " coarse grid vectors" <<std::endl;

    assert(vecs[0].Grid()==coarse_grid);

    int _ndimension = coarse_grid->_ndimension;

    uint64_t sz = blas.size();

    Coordinate coarse_rdimensions = coarse_grid->_rdimensions;
    
    for(int v=0;v<vecs.size();v++){

      //      std::cout << " BlockProjector importing coarse vector"<<v<<" "<<norm2(vecs[v])<<std::endl;
      autoView( coarseData   , vecs[v], AcceleratorRead);

      auto blasData_p  = &blas[0];
      auto coarseData_p  = &coarseData[0];

      int64_t osites = coarse_grid->oSites();

      // loop over fine sites
      const int Nsimd = vobj::Nsimd();
      uint64_t cwords=sizeof(typename vobj::scalar_object)/sizeof(scalar);
      assert(cwords==nbasis);
      
      accelerator_for(sc,osites,Nsimd,{
#ifdef GRID_SIMT
        {
	  int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
	  for(int lane=0;lane<Nsimd;lane++) {
#endif
           // C_br per site
	    int64_t blas_site = (lane*osites + sc)*nvec*cwords + v*cwords;
	    
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
    int nvec = vecs.size();
    typedef typename vobj::scalar_object coarse_scalar_object;
    //    std::cout << GridLogMessage<<" BlockProjector exporting "<<nvec<< " coarse grid vectors" <<std::endl;

    assert(vecs[0].Grid()==coarse_grid);

    int _ndimension = coarse_grid->_ndimension;
    
    uint64_t sz = blas.size();

    Coordinate coarse_rdimensions = coarse_grid->_rdimensions;
    
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
      assert(cwords==nbasis);
      
      accelerator_for(sc,osites,Nsimd,{
	  // Wrap in a macro "FOR_ALL_LANES(lane,{ ... });
#ifdef GRID_SIMT
        {
	  int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
	  for(int lane=0;lane<Nsimd;lane++) {
#endif
	    int64_t blas_site = (lane*osites + sc)*nvec*cwords + v*cwords;
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
    int nrhs=fine.size();
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    //    std::cout << "blockProject nbasis " <<nbasis<<" " << _nbasis<<std::endl;
    assert(nbasis==_nbasis);
    
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
    for(int c=0;c<coarse_vol;c++){
      // BLAS_V[coarse_vol][nbasis][block_vol][words]
      // BLAS_F[coarse_vol][nrhs][block_vol][words]
      // BLAS_C[coarse_vol][nrhs][nbasis]
      scalar * Vh = & BLAS_V[c*nbasis*block_vol*words];
      scalar * Fh = & BLAS_F[c*nrhs*block_vol*words];
      scalar * Ch = & BLAS_C[c*nrhs*nbasis];

      acceleratorPut(Vd[c],Vh);
      acceleratorPut(Fd[c],Fh);
      acceleratorPut(Cd[c],Ch);
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
    int nrhs=fine.size();
    int _nbasis = sizeof(typename cobj::scalar_object)/sizeof(scalar);
    assert(nbasis==_nbasis);
    
    BLAS_F.resize (fine_vol * words * nrhs );
    BLAS_C.resize (coarse_vol * nbasis * nrhs );

    ImportCoarseGridVectors(coarse, BLAS_C);

    GridBLAS BLAS;

    deviceVector<scalar *> Vd(coarse_vol);
    deviceVector<scalar *> Fd(coarse_vol);
    deviceVector<scalar *> Cd(coarse_vol);

    for(int c=0;c<coarse_vol;c++){
      // BLAS_V[coarse_vol][nbasis][block_vol][words]
      // BLAS_F[coarse_vol][nrhs][block_vol][words]
      // BLAS_C[coarse_vol][nrhs][nbasis]
      scalar * Vh = & BLAS_V[c*nbasis*block_vol*words];
      scalar * Fh = & BLAS_F[c*nrhs*block_vol*words];
      scalar * Ch = & BLAS_C[c*nrhs*nbasis];
      acceleratorPut(Vd[c],Vh);
      acceleratorPut(Fd[c],Fh);
      acceleratorPut(Cd[c],Ch);
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
};

NAMESPACE_END(Grid);
