/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: MomentumProject.h

    Copyright (C) 2025

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
   MultiMomProject

   Import vectors -> nxyz x (ncomponent x nt)
   Import complex phases -> nmom x nxy

   apply = via (possibly batched) GEMM
*/
template<class Field, class ComplexField>
class MomentumProject
{
public:

  typedef typename Field::scalar_type   scalar;
  typedef typename Field::scalar_object scalar_object;

  GridBase *grid;
  uint64_t nmom;
  uint64_t nxyz;
  uint64_t nt;
  uint64_t nbtw;
  uint64_t words;

  deviceVector<scalar> BLAS_V;      // 
  deviceVector<scalar> BLAS_M;      // 
  deviceVector<scalar> BLAS_P;      // 
  
  MomentumProject(){};
 ~MomentumProject(){ Deallocate(); };
  
  void Deallocate(void)
  {
    grid=nullptr;
    nmom=0;
    nxyz=0;
    nt=0;
    nbtw=0;
    words=0;
    BLAS_V.resize(0);
    BLAS_M.resize(0);
    BLAS_P.resize(0);
  }
  void Allocate(int _nmom,GridBase *_grid)
  {
    grid=_grid;
    Coordinate ldims = grid->LocalDimensions();

    nmom=_nmom;
    nt   = ldims[grid->Nd()-1];
    nxyz = grid->lSites()/nt;
    words = sizeof(scalar_object)/sizeof(scalar);
    nbtw = nt * words;

    BLAS_V.resize (nxyz * nt * words );
    BLAS_M.resize (nmom * nxyz       );
    BLAS_P.resize (nmom * nt * words );
  }
  void ImportMomenta(const std::vector <ComplexField> &momenta)
  {
    GRID_ASSERT(momenta.size()==nmom);
    //    might as well just make the momenta here
    typedef typename Field::vector_object vobj;

    int nd = grid->_ndimension;

    uint64_t sz = BLAS_M.size();

    GRID_ASSERT(momenta.size()==nmom)
    GRID_ASSERT(momenta[0].Grid()==grid);
    GRID_ASSERT(sz = nxyz * nmom);
    
    Coordinate rdimensions = grid->_rdimensions;
    Coordinate ldims       = grid->LocalDimensions();
    int64_t osites         = grid->oSites();
    Coordinate simd        = grid->_simd_layout;
    const int Nsimd        = vobj::Nsimd();
    uint64_t lwords        = words; // local variable for copy in to GPU
    int64_t Nxyz = nxyz;
    auto blasData_p  = &BLAS_M[0];
    for(int m=0;m<momenta.size();m++){

      autoView( Data   , momenta[m], AcceleratorRead);
      auto Data_p  = &Data[0];

      accelerator_for(xyz,nxyz,1,{
	  //////////////////////////////////////////
	  // isite -- map lane within buffer to lane within lattice
	  ////////////////////////////////////////////
	    Coordinate lcoor(nd,0);
	    Lexicographic::CoorFromIndex(lcoor,xyz,ldims);
	    
	    Coordinate icoor(nd);
	    Coordinate ocoor(nd);
	    for (int d = 0; d < nd; d++) {
	      icoor[d] = lcoor[d]/rdimensions[d];
	      ocoor[d] = lcoor[d]%rdimensions[d];
	    }
	    int64_t osite;
	    int64_t isite;
	    Lexicographic::IndexFromCoor(ocoor,osite,rdimensions);
	    Lexicographic::IndexFromCoor(icoor,isite,simd);
	    
	    // BLAS_M[nmom][slice_vol]
	    // Fortran Column major BLAS layout is M_xyz,mom
	    scalar data = extractLane(isite,Data[osite]);
	    uint64_t idx = xyz+m*Nxyz;
	    blasData_p[idx] = data;
	});
    }
  }
  void ImportVector(Field &vec)
  {
    typedef typename Field::vector_object vobj;

    int nd = grid->_ndimension;

    uint64_t sz = BLAS_V.size();

    GRID_ASSERT(sz = nxyz * words * nt);
    
    Coordinate rdimensions = grid->_rdimensions;
    Coordinate ldims= grid->LocalDimensions();
    int64_t osites = grid->oSites();
    Coordinate simd = grid->_simd_layout;
    const int Nsimd = vobj::Nsimd();
    uint64_t lwords= words; // local variable for copy in to GPU

    auto blasData_p  = &BLAS_V[0];
    autoView( Data   , vec, AcceleratorRead);
    auto Data_p  = &Data[0];

    int64_t nwords = words;// for capture
    int64_t Nt     = nt;// for capture
    
    accelerator_for(sf,osites,Nsimd,{
#ifdef GRID_SIMT
        {
	  int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
	  for(int lane=0;lane<Nsimd;lane++) {
#endif
	  //////////////////////////////////////////
	  // isite -- map lane within buffer to lane within lattice
	  ////////////////////////////////////////////
	    Coordinate lcoor(nd,0);
	    Coordinate icoor(nd);
	    Coordinate ocoor(nd);
	    
	    Lexicographic::CoorFromIndex(icoor,lane,simd);
	    Lexicographic::CoorFromIndex(ocoor,sf,rdimensions);

	  
	    int64_t l_xyz = 0;
	    for (int d = 0; d < nd; d++) {
	      lcoor[d] = rdimensions[d]*icoor[d] + ocoor[d];
	    }
	    uint64_t l_t   = lcoor[nd-1];

	    Coordinate xyz_coor = lcoor;
	    xyz_coor[nd-1] =0;
	    Lexicographic::IndexFromCoor(xyz_coor,l_xyz,ldims);

	    
	    scalar_object data = extractLane(lane,Data[sf]);
	    scalar *data_words = (scalar *) &data;
	    for(int w = 0 ; w < nwords; w++) {
	      // BLAS_V[slice_vol][nt][words]
	      // Fortran Column major BLAS layout is V_(t,w)_xyz
	      uint64_t idx = w+l_t*nwords + l_xyz * nwords * Nt;
	      blasData_p[idx] = data_words[w];
	    }
#ifdef GRID_SIMT
	}
#else
	}
#endif
	});
  }
  void ExportMomentumProjection(std::vector<typename Field::scalar_object> &projection)
  {
    projection.resize(nmom*nt);
    acceleratorCopyFromDevice(&BLAS_P[0],(scalar *)&projection[0],BLAS_P.size()*sizeof(scalar));
    // Could decide on a layout late?
  }

  // Row major layout "C" order:
  // BLAS_V[slice_vol][nt][words]
  // BLAS_M[nmom][slice_vol]
  // BLAS_P[nmom][nt][words]
  //
  // Fortran Column major BLAS layout is V_(w,t)_xyz
  // Fortran Column major BLAS layout is M_xyz,mom
  // Fortran Column major BLAS layout is P_(w,t),mom
  //
  // Projected
  //
  // P = (V * M)_(w,t),mom
  //
  void Project(Field &data,std::vector< typename Field::scalar_object > & projected_gdata)
  {
    this->ImportVector(data);

    std::vector< typename Field::scalar_object > projected_planes;

    deviceVector<scalar *> Vd(1);
    deviceVector<scalar *> Md(1);
    deviceVector<scalar *> Pd(1);

    scalar * Vh = & BLAS_V[0];
    scalar * Mh = & BLAS_M[0];
    scalar * Ph = & BLAS_P[0];

    acceleratorPut(Vd[0],Vh);
    acceleratorPut(Md[0],Mh);
    acceleratorPut(Pd[0],Ph);

    GridBLAS BLAS;

    /////////////////////////////////////////
    // P_im = VMmx . Vxi
    /////////////////////////////////////////
    BLAS.gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N, 
    		     words*nt,nmom,nxyz,
		     scalar(1.0),
		     Vd,
		     Md,
		     scalar(0.0),  // wipe out result
		     Pd);
    BLAS.synchronise();

    ExportMomentumProjection(projected_planes); // resizes

    /////////////////////////////////
    // Reduce across MPI ranks
    /////////////////////////////////
    int nd = grid->Nd();
    int gt = grid->GlobalDimensions()[nd-1];
    int lt = grid->LocalDimensions()[nd-1];
    projected_gdata.resize(gt*nmom);
    for(int t=0;t<gt*nmom;t++){ // global Nt array with zeroes for stuff not on this node
      projected_gdata[t]=Zero();
    }
    for(int t=0;t<lt;t++){
    for(int m=0;m<nmom;m++){
      int st = grid->LocalStarts()[nd-1];
      projected_gdata[t+st + gt*m] = projected_planes[t+lt*m];
    }}
    grid->GlobalSumVector((scalar *)&projected_gdata[0],gt*nmom*words);
  }
};

NAMESPACE_END(Grid);
