/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/GeneralCoarsenedMatrixMultiRHS.h

    Copyright (C) 2015

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


// Fine Object == (per site) type of fine field
// nbasis      == number of deflation vectors
template<class Fobj,class CComplex,int nbasis>
class MultiGeneralCoarsenedOperatorV2 : public SparseMatrixBase<Lattice<iVector<CComplex,nbasis > > >  {
public:
  typedef typename CComplex::scalar_object SComplex;
  typedef GeneralCoarsenedMatrix<Fobj,CComplex,nbasis> GeneralCoarseOp;
  typedef MultiGeneralCoarsenedOperatorV2<Fobj,CComplex,nbasis> MultiGeneralCoarseOp;

  typedef iVector<CComplex,nbasis >           siteVector;
  typedef iMatrix<CComplex,nbasis >           siteMatrix;
  typedef iVector<SComplex,nbasis >           calcVector;
  typedef iMatrix<SComplex,nbasis >           calcMatrix;
  typedef Lattice<iScalar<CComplex> >         CoarseComplexField;
  typedef Lattice<siteVector>                 CoarseVector;
  typedef Lattice<iMatrix<CComplex,nbasis > > CoarseMatrix;
  typedef iMatrix<CComplex,nbasis >  Cobj;
  typedef iVector<CComplex,nbasis >  Cvec;
  typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj >        FineField;
  typedef CoarseVector Field;

  // Block operations on the fine vectors carry the fine layout, which need
  // not be the coarse one
  typedef decltype(innerProduct(Fobj(),Fobj())) FineInner;
  typedef Lattice<FineInner>                    FineComplexField;
  typedef Lattice<FineInner>                    BlockComplexField;

  ////////////////////
  // Data members
  //
  // Nrhs independent: the D dimensional coarse grid, the geometry, the padded
  // cell that supplies the stencil grid, the stencil, and the matrix elements.
  //
  // Nrhs dependent: the D+1 grid, its padded cell, and the BLAS B/C buffers
  // with their pointer tables. Owned by SetNRHS().
  ////////////////////
  GridCartesian *       _CoarseGrid;        // D dimensional
  NonLocalStencilGeometry geom;
  NonLocalStencilGeometry geom_srhs;
  PaddedCell CellD;                         // D dimensional, supplies stencil grid
  GeneralLocalStencil Stencil;              // D dimensional

  int                   _Nrhs;
  GridCartesian *       _CoarseGridMulti;   // D+1 dimensional, SetNRHS
  PaddedCell *          CellMulti;          // D+1 dimensional, SetNRHS

  deviceVector<calcVector> BLAS_B;
  deviceVector<calcVector> BLAS_C;
  std::vector<deviceVector<calcMatrix> > BLAS_A;

  std::vector<deviceVector<ComplexD *> > BLAS_AP;
  std::vector<deviceVector<ComplexD *> > BLAS_BP;
  deviceVector<ComplexD *>               BLAS_CP;

  ///////////////////////
  // Interface
  ///////////////////////
  GridBase      * Grid(void)           { CheckGridSet(); return _CoarseGridMulti; };
  GridCartesian * CoarseGrid(void)     { CheckGridSet(); return _CoarseGridMulti; };
  GridCartesian * CoarseGridD(void)    { return _CoarseGrid; };        // lower dimensional grid
  int             Nrhs(void)           { CheckGridSet(); return _Nrhs; };

  void CheckGridSet(void)
  {
    if ( _CoarseGridMulti == nullptr ) {
      std::cout << GridLogError
		<< "MultiGeneralCoarsenedOperatorV2: the multiRHS grid has not been set."
		<< std::endl;
      std::cout << GridLogError
		<< "  Call SetGrid(CoarseGridMulti) with the D+1 dimensional grid your"
		<< std::endl;
      std::cout << GridLogError
		<< "  coarse vectors live on, before Grid(), Nrhs() or M()."
		<< std::endl;
      GRID_ASSERT(_CoarseGridMulti != nullptr);
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Bilingual accessors, matching GeneralCoarsenedMatrix. Note Grid() is the
  // D+1 multiRHS grid here, so a consumer wanting the space the elements live
  // on must ask for CoarseGridD().
  //////////////////////////////////////////////////////////////////////////
  NonLocalStencilGeometry & Geometry(void)      { return geom_srhs; };
  void ExtractMatrix(int p,CoarseMatrix &A)     { BLAStoGrid(A,BLAS_A[p]); };

  // I/O on the operator matrices, via the BLAS layout array. The parameter is
  // a vector over the geometry points; the body indexes A[p].
  void SetMatrix (int p,std::vector<CoarseMatrix> & A)
  {
    GRID_ASSERT(A.size()==geom_srhs.npoint);
    GridtoBLAS(A[p],BLAS_A[p]);
  }
  void GetMatrix (int p,std::vector<CoarseMatrix> & A)
  {
    GRID_ASSERT(A.size()==geom_srhs.npoint);
    BLAStoGrid(A[p],BLAS_A[p]);
  }
  void CopyMatrix (GeneralCoarseOp &_Op)
  {
    for(int p=0;p<geom.npoint;p++){
      auto Aup = _Op.Cell.Extract(_Op._A[p]);
      //Unpadded
      GridtoBLAS(Aup,BLAS_A[p]);
    }
  }
  /*
  void CheckMatrix (GeneralCoarseOp &_Op)
  {
    std::cout <<"************* Checking the little direc operator mRHS"<<std::endl;
    for(int p=0;p<geom.npoint;p++){
      //Unpadded
      auto Aup = _Op.Cell.Extract(_Op._A[p]);
      auto Ack = Aup;
      BLAStoGrid(Ack,BLAS_A[p]);
      std::cout << p<<" Ack "<<norm2(Ack)<<std::endl;
      std::cout << p<<" Aup "<<norm2(Aup)<<std::endl;
    }
    std::cout <<"************* "<<std::endl;
  }
  */
  
  ///////////////////////////////////////////////////////////////////////////
  // Constructor takes the D dimensional coarse grid. Everything built here
  // is independent of Nrhs, in particular the matrix elements, which must
  // survive a change of Nrhs untouched.
  ///////////////////////////////////////////////////////////////////////////
  MultiGeneralCoarsenedOperatorV2(NonLocalStencilGeometry &_geom,GridCartesian *CoarseGrid) :
    _CoarseGrid(CoarseGrid),
    geom_srhs(_geom),
    geom(CoarseGrid,_geom.hops,_geom.skip),
    CellD(geom.Depth(),CoarseGrid),
    Stencil(CellD.grids.back(),geom.shifts), // D dimensional padded cell stencil
    _Nrhs(-1),
    _CoarseGridMulti(nullptr),
    CellMulti(nullptr)
  {
    int32_t unpadded_sites = _CoarseGrid->lSites();

    /////////////////////////////////////////////////
    // Matrix elements and their pointer table
    /////////////////////////////////////////////////
    BLAS_A.resize(geom.npoint);
    BLAS_AP.resize(geom.npoint);
    for(int p=0;p<geom.npoint;p++){
      BLAS_A[p].resize (unpadded_sites); // no ghost zone, npoint elements
      BLAS_AP[p].resize(unpadded_sites);
    }

    // Site identity mapping for A
    for(int p=0;p<geom.npoint;p++){
      for(int ss=0;ss<unpadded_sites;ss++){
	ComplexD *ptr = (ComplexD *)&BLAS_A[p][ss];
	acceleratorPut(BLAS_AP[p][ss],ptr);
      }
    }
  }

  virtual ~MultiGeneralCoarsenedOperatorV2()
  {
    ReleaseGrid();
  }

  ///////////////////////////////////////////////////////////////////////////
  // Free everything SetGrid allocated. The D+1 grid is borrowed from the
  // caller and is never deleted here. Safe to call repeatedly and before
  // the destructor.
  ///////////////////////////////////////////////////////////////////////////
  void ReleaseGrid(void)
  {
    if ( CellMulti != nullptr ) { delete CellMulti; CellMulti = nullptr; }

    _CoarseGridMulti = nullptr;  // borrowed, not owned
    _Nrhs            = -1;

    BLAS_B.resize(0);
    BLAS_C.resize(0);
    for(int p=0;p<BLAS_BP.size();p++){
      BLAS_BP[p].resize(0);
    }
    BLAS_BP.resize(0);
    BLAS_CP.resize(0);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Everything that depends on Nrhs. Idempotent; lazy called on demand.
  //
  // The stencil lives on the D dimensional padded grid. PaddedCell only pads
  // a dimension when it is distributed, and the rhs direction never is, so
  // the D+1 padded grid is exactly Nrhs copies of the D dimensional padded
  // grid with rhs innermost. The neighbour offset therefore carries an Nrhs
  // factor, in the same way the Nsimd factor is carried.
  ///////////////////////////////////////////////////////////////////////////
  void SetGrid(GridCartesian *CoarseGridMulti)
  {
    GRID_ASSERT(CoarseGridMulti != nullptr);

    if ( CoarseGridMulti == _CoarseGridMulti ) return; // idempotent on identity

    ReleaseGrid();

    /////////////////////////////////////////////////
    // The D+1 grid is supplied and owned by the caller. Two operators over
    // the same coarse space must share one grid object or their fields will
    // not conform, so this is never manufactured internally.
    /////////////////////////////////////////////////
    int nd = _CoarseGrid->_ndimension;

    GRID_ASSERT(CoarseGridMulti->_ndimension == nd+1);
    GRID_ASSERT(CoarseGridMulti->_processors[0] == 1);   // rhs is not distributed
    for(int d=0;d<nd;d++){
      GRID_ASSERT(CoarseGridMulti->_fdimensions[d+1] == _CoarseGrid->_fdimensions[d]);
      GRID_ASSERT(CoarseGridMulti->_processors [d+1] == _CoarseGrid->_processors [d]);
      GRID_ASSERT(CoarseGridMulti->_simd_layout[d+1] == _CoarseGrid->_simd_layout[d]);
    }

    _CoarseGridMulti = CoarseGridMulti;
    _Nrhs            = CoarseGridMulti->_fdimensions[0];
    GRID_ASSERT(_Nrhs>=1);

    int nrhs = _Nrhs;

    CellMulti = new PaddedCell(geom.Depth(),_CoarseGridMulti);

    int32_t padded_sites   = CellD.grids.back()->lSites(); // D dimensional
    int32_t unpadded_sites = _CoarseGrid->lSites();        // D dimensional

    // The neighbour offset multiplication by nrhs is exact only if the D+1
    // padded volume is nrhs copies of the D dimensional one. Check it.
    GRID_ASSERT(CellMulti->grids.back()->lSites() == nrhs*padded_sites);
    GRID_ASSERT(_CoarseGridMulti->lSites()        == nrhs*unpadded_sites);

    /////////////////////////////////////////////////
    // Device data vector storage
    /////////////////////////////////////////////////
    BLAS_B.resize(nrhs *padded_sites);   // includes ghost zone
    BLAS_C.resize(nrhs *unpadded_sites); // no ghost zone
    BLAS_BP.resize(geom.npoint);
    for(int p=0;p<geom.npoint;p++){
      BLAS_BP[p].resize(unpadded_sites);
    }
    BLAS_CP.resize(unpadded_sites);

    // Site identity mapping for C
    for(int ss=0;ss<unpadded_sites;ss++){
      ComplexD *ptr = (ComplexD *)&BLAS_C[ss*nrhs];
      acceleratorPut(BLAS_CP[ss],ptr);
    }

    // Neighbour table is more complicated
    int32_t j=0; // Interior point counter (unpadded)
    for(int32_t s=0;s<padded_sites;s++){ // D volume, padded
      int ghost_zone=0;
      for(int32_t point = 0 ; point < geom.npoint; point++){
	int i=s*geom.npoint+point;
	if( Stencil._entries[i]._wrap ) { // stencil is indexed by the oSite of the D dim grid
	  ghost_zone=1; // If general stencil wrapped in any direction, wrap=1
	}
      }

      if( ghost_zone==0) {
	for(int32_t point = 0 ; point < geom.npoint; point++){
	  int i=s*geom.npoint+point;
 	  int32_t nbr = Stencil._entries[i]._offset*CComplex::Nsimd(); // oSite -> lSite, D dim
	  nbr = nbr*nrhs;                                              // D -> D+1, rhs innermost
	  GRID_ASSERT(nbr<BLAS_B.size());
	  ComplexD * ptr = (ComplexD *)&BLAS_B[nbr];
	  acceleratorPut(BLAS_BP[point][j],ptr); // neighbour indexing in ghost zone volume
	}
	j++;
      }
    }
    GRID_ASSERT(j==unpadded_sites);
  }
  template<class vobj> void GridtoBLAS(const Lattice<vobj> &from,deviceVector<typename vobj::scalar_object> &to)
  {
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  GridBase *Fg = from.Grid();
  GRID_ASSERT(!Fg->_isCheckerBoarded);
  int nd = Fg->_ndimension;

  to.resize(Fg->lSites());

  Coordinate LocalLatt = Fg->LocalDimensions();
  size_t nsite = 1;
  for(int i=0;i<nd;i++) nsite *= LocalLatt[i];

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // do the index calc on the GPU
  ////////////////////////////////////////////////////////////////////////////////////////////////
  Coordinate f_ostride = Fg->_ostride;
  Coordinate f_istride = Fg->_istride;
  Coordinate f_rdimensions = Fg->_rdimensions;

  autoView(from_v,from,AcceleratorRead);
  auto to_v = &to[0];

  const int words=sizeof(vobj)/sizeof(vector_type);
  accelerator_for(idx,nsite,1,{
      
      Coordinate from_coor, base;
      Lexicographic::CoorFromIndex(base,idx,LocalLatt);
      for(int i=0;i<nd;i++){
	from_coor[i] = base[i];
      }
      int from_oidx = 0; for(int d=0;d<nd;d++) from_oidx+=f_ostride[d]*(from_coor[d]%f_rdimensions[d]);
      int from_lane = 0; for(int d=0;d<nd;d++) from_lane+=f_istride[d]*(from_coor[d]/f_rdimensions[d]);

      const vector_type* from = (const vector_type *)&from_v[from_oidx];
      scalar_type* to = (scalar_type *)&to_v[idx];
      
      scalar_type stmp;
      for(int w=0;w<words;w++){
	stmp = getlane(from[w], from_lane);
	to[w] = stmp;
      }
    });
  }    
  template<class vobj> void BLAStoGrid(Lattice<vobj> &grid,deviceVector<typename vobj::scalar_object> &in)
  {
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  GridBase *Tg = grid.Grid();
  GRID_ASSERT(!Tg->_isCheckerBoarded);
  int nd = Tg->_ndimension;
  
  GRID_ASSERT(in.size()==Tg->lSites());

  Coordinate LocalLatt = Tg->LocalDimensions();
  size_t nsite = 1;
  for(int i=0;i<nd;i++) nsite *= LocalLatt[i];

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // do the index calc on the GPU
  ////////////////////////////////////////////////////////////////////////////////////////////////
  Coordinate t_ostride = Tg->_ostride;
  Coordinate t_istride = Tg->_istride;
  Coordinate t_rdimensions = Tg->_rdimensions;

  autoView(to_v,grid,AcceleratorWrite);
  auto from_v = &in[0];

  const int words=sizeof(vobj)/sizeof(vector_type);
  accelerator_for(idx,nsite,1,{
      
      Coordinate to_coor, base;
      Lexicographic::CoorFromIndex(base,idx,LocalLatt);
      for(int i=0;i<nd;i++){
	to_coor[i] = base[i];
      }
      int to_oidx = 0; for(int d=0;d<nd;d++) to_oidx+=t_ostride[d]*(to_coor[d]%t_rdimensions[d]);
      int to_lane = 0; for(int d=0;d<nd;d++) to_lane+=t_istride[d]*(to_coor[d]/t_rdimensions[d]);

      vector_type* to = (vector_type *)&to_v[to_oidx];
      scalar_type* from = (scalar_type *)&from_v[idx];
      
      scalar_type stmp;
      for(int w=0;w<words;w++){
	stmp=from[w];
	putlane(to[w], stmp, to_lane);
      }
    });
  }
  ///////////////////////////////////////////////////////////////////////////
  // Shared by both CoarsenOperator variants
  //
  //     conj(pha[block]) proj[k (which mom)][j (basis vec cpt)][block]
  //       =  \sum_{l in ball}  e^{i q_k . delta_l} < phi_{block,j} | MdagM | phi_{(block+delta_l),i} >
  //       =  \sum_{l in ball} e^{iqk.delta_l} A_ji^{b.b+l}
  //       = M_{kl} A_ji^{b.b+l}
  //
  //     Where q_k = delta_k . (2*M_PI/global_nb[mu])
  //     Then A{ji}^{b,b+l} = M^{-1}_{lm} ComputeProj_{m,b,i,j}
  ///////////////////////////////////////////////////////////////////////////
  void CoarsenFourierMatrix(GridBase *CoarseGrid,Eigen::MatrixXcd &invMkl)
  {
    const int npoint = geom_srhs.npoint;
    Coordinate clatt = CoarseGrid->GlobalDimensions();
    int Nd = CoarseGrid->Nd();

    Eigen::MatrixXcd Mkl = Eigen::MatrixXcd::Zero(npoint,npoint);
    ComplexD ci(0.0,1.0);
    for(int k=0;k<npoint;k++){ // Loop over momenta
      for(int l=0;l<npoint;l++){ // Loop over nbr relative
	ComplexD phase(0.0,0.0);
	for(int mu=0;mu<Nd;mu++){
	  RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	  phase=phase+TwoPiL*geom_srhs.shifts[k][mu]*geom_srhs.shifts[l][mu];
	}
	phase=exp(phase*ci);
	Mkl(k,l) = phase;
      }
    }
    invMkl = Mkl.inverse();
  }

  ///////////////////////////////////////////////////////////////////////////
  // blockOrthogonalise and blockZAXPY are block operations on the fine
  // vectors, using a coarse shaped field only as an index set. They need a
  // grid carrying the fine SIMD layout, which the coarse space no longer
  // does. Constructed local to the caller so it cannot be mistaken for the
  // coarse grid.
  ///////////////////////////////////////////////////////////////////////////
  void CoarsenBlockGridLayout(GridBase *grid,GridBase *CoarseGrid,
			      Coordinate &latt,Coordinate &simd,Coordinate &mpi)
  {
    int nd = CoarseGrid->_ndimension;
    latt.resize(nd); simd.resize(nd); mpi.resize(nd);
    for(int d=0;d<nd;d++){
      latt[d] = CoarseGrid->_fdimensions[d];
      simd[d] = grid->_simd_layout[d];
      mpi [d] = CoarseGrid->_processors[d];
    }
  }

  // D+1 coarse grid holding the batch, rhs innermost and unvectorised
  void CoarsenBatchGridLayout(GridBase *CoarseGrid,int batch,
			      Coordinate &latt,Coordinate &simd,Coordinate &mpi)
  {
    latt.resize(1,batch); simd.resize(1,1); mpi.resize(1,1);
    latt[0]=batch; simd[0]=1; mpi[0]=1;
    for(int d=0;d<CoarseGrid->_ndimension;d++){
      latt.push_back(CoarseGrid->_fdimensions[d]);
      simd.push_back(CoarseGrid->_simd_layout[d]);
      mpi .push_back(CoarseGrid->_processors[d]);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // The Fourier inverse needs the phase in the coarse layout and the basis
  // phasing needs it in the fine layout; each is built from its own
  // coordinates rather than transferred.
  ///////////////////////////////////////////////////////////////////////////
  void CoarsenPhases(GridBase *grid,GridBase *CoarseGrid,GridCartesian *BlockGrid,
		     std::vector<CoarseComplexField> &pha,
		     std::vector<FineComplexField> &phaF)
  {
    const int npoint = geom_srhs.npoint;
    Coordinate clatt = CoarseGrid->GlobalDimensions();
    int Nd = CoarseGrid->Nd();
    ComplexD ci(0.0,1.0);

    typedef typename CComplex::scalar_type SComplex;
    FineComplexField one(grid); one=SComplex(1.0);
    FineComplexField zz(grid);  zz = Zero();
    BlockComplexField pha_blk (BlockGrid);
    BlockComplexField blk_coor(BlockGrid);

    for(int p=0;p<npoint;p++){ // Loop over momenta in npoint
      CoarseComplexField coor(CoarseGrid);
      pha[p] =Zero();
      pha_blk=Zero();
      for(int mu=0;mu<Nd;mu++){
	RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	LatticeCoordinate(coor,mu);
	pha[p]  = pha[p]  + (TwoPiL * geom_srhs.shifts[p][mu]) * coor;
	LatticeCoordinate(blk_coor,mu);
	pha_blk = pha_blk + (TwoPiL * geom_srhs.shifts[p][mu]) * blk_coor;
      }
      pha[p] =exp(pha[p] *ci);
      pha_blk=exp(pha_blk*ci);

      blockZAXPY(phaF[p],pha_blk,one,zz);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Remove the bulk phase from the batch of coarse projections and
  // accumulate the Fourier inverse into A. Both variants reach here with
  // TmpProj in the same batch coarse order, so this is shared verbatim.
  ///////////////////////////////////////////////////////////////////////////
  void CoarsenAccumulate(int p,int i0,int nbv,int batch,
			 Eigen::MatrixXcd &invMkl,
			 std::vector<CoarseComplexField> &pha,
			 CoarseComplexField &phaB,
			 CoarseVector &TmpProj,
			 std::vector<CoarseMatrix> &_A,
			 GridBase *CoarseGrid)
  {
    typedef typename CComplex::scalar_type SComplex;
    const int npoint = geom_srhs.npoint;

    for(int b=0;b<batch;b++) InsertSliceFast(pha[p],phaB,b,0);
    TmpProj = conjugate(phaB)*TmpProj;

    int osites=CoarseGrid->oSites();
    for(int k=0;k<npoint;k++){
      SComplex sc(invMkl(p,k).real(),invMkl(p,k).imag());
      CComplex coef(sc);
      autoView( A_v  , _A[k], AcceleratorWrite);
      autoView( TP_v , TmpProj, AcceleratorRead);
      accelerator_for(sss, osites, 1, {
	  for(int b=0;b<nbv;b++){
	    for(int j=0;j<nbasis;j++){
	      A_v[sss](i0+b,j) = A_v[sss](i0+b,j) + coef*TP_v[b+batch*sss](j);
	    }
	  }
      });
    }
  }

  void CoarsenReport(RealD tphase,RealD tphaseBZ,RealD tslice,
		     RealD tmat,RealD tproj,RealD tinv)
  {
    std::cout << GridLogMessage<<"CoarsenOperator phase  "<<tphase<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator phaseBZ "<<tphaseBZ<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator slice  "<<tslice <<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator mat    "<<tmat <<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator proj   "<<tproj<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator inv    "<<tinv<<" us"<<std::endl;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Coarsen a NATIVELY multiRHS fine operator.
  //
  // linop acts on the D+1 dimensional fine grid FineGridMulti, with the batch
  // of phased basis vectors carried in the rhs direction. A single RHS
  // operator can be promoted with MrhsPromotedOperator, but that pays an
  // ExtractSlice/InsertSlice pair per rhs; prefer the single RHS variant
  // below in that case.
  ///////////////////////////////////////////////////////////////////////////
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       GridCartesian *FineGridMulti,
		       std::vector<FineField> &Subspace,
		       GridBase *CoarseGrid)
  {
    RealD tproj=0.0, tmat=0.0, tphase=0.0, tphaseBZ=0.0, tslice=0.0, tinv=0.0;

    std::cout << GridLogMessage<< "GeneralCoarsenMatrixMrhs (multiRHS fine operator)"<< std::endl;

    GRID_ASSERT(Subspace.size()==nbasis);
    GridBase *grid = Subspace[0].Grid();

    GRID_ASSERT(FineGridMulti->_ndimension == grid->_ndimension+1);
    GRID_ASSERT(FineGridMulti->_processors[0] == 1);
    for(int d=0;d<grid->_ndimension;d++){
      GRID_ASSERT(FineGridMulti->_fdimensions[d+1] == grid->_fdimensions[d]);
      GRID_ASSERT(FineGridMulti->_processors [d+1] == grid->_processors [d]);
    }
    int batch = FineGridMulti->_fdimensions[0];

    Coordinate blatt,bsimd,bmpi;
    CoarsenBlockGridLayout(grid,CoarseGrid,blatt,bsimd,bmpi);
    GridCartesian BlockGrid(blatt,bsimd,bmpi);

    BlockComplexField InnerProd(&BlockGrid);
    blockOrthogonalise(InnerProd,Subspace);

    MultiRHSBlockProject<Lattice<Fobj> >    Projector;
    Projector.Allocate(nbasis,grid,CoarseGrid);
    Projector.ImportBasis(Subspace);

    const int npoint = geom_srhs.npoint;

    Eigen::MatrixXcd invMkl;
    CoarsenFourierMatrix(CoarseGrid,invMkl);

    FineField phaV(grid);
    std::vector<FineComplexField>   phaF(npoint,grid);
    std::vector<CoarseComplexField> pha (npoint,CoarseGrid);

    tphase=-usecond();
    CoarsenPhases(grid,CoarseGrid,&BlockGrid,pha,phaF);
    tphase+=usecond();

    std::vector<CoarseMatrix> _A;
    _A.resize(npoint,CoarseGrid);
    for(int k=0;k<npoint;k++) _A[k] = Zero();

    Coordinate cmlatt,cmsimd,cmmpi;
    CoarsenBatchGridLayout(CoarseGrid,batch,cmlatt,cmsimd,cmmpi);
    GridCartesian CoarseBatchGrid(cmlatt,cmsimd,cmmpi);

    CoarseVector       TmpProj(&CoarseBatchGrid);
    CoarseComplexField phaB(&CoarseBatchGrid);

    FineField hi_in (FineGridMulti);
    FineField hi_out(FineGridMulti);
    FineField zzF(grid); zzF = Zero();

    for(int i0=0;i0<nbasis;i0+=batch){ // Loop over batches of basis vectors

      int nbv = MIN(batch,nbasis-i0);
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i0<<"/"<<nbasis<< std::endl;

      for(int p=0;p<npoint;p++){ // Loop over momenta

	// One phase, applied to the whole batch. Tail slices are zeroed so
	// the operator never sees undefined data.
	for(int b=0;b<nbv;b++){
	  tphaseBZ-=usecond();
	  phaV = phaF[p]*Subspace[i0+b];
	  tphaseBZ+=usecond();
	  tslice-=usecond();
	  InsertSliceFast(phaV,hi_in,b,0);
	  tslice+=usecond();
	}
	tslice-=usecond();
	for(int b=nbv;b<batch;b++){
	  InsertSliceFast(zzF,hi_in,b,0);
	}
	tslice+=usecond();

	tmat-=usecond();
	linop.Op(hi_in,hi_out);
	tmat+=usecond();

	tproj-=usecond();
	Projector.blockProject(hi_out,TmpProj);
	tproj+=usecond();

	tinv-=usecond();
	CoarsenAccumulate(p,i0,nbv,batch,invMkl,pha,phaB,TmpProj,_A,CoarseGrid);
	tinv+=usecond();
      }
    }

    for(int p=0;p<npoint;p++){
      GridtoBLAS(_A[p],BLAS_A[p]);
    }
    CoarsenReport(tphase,tphaseBZ,tslice,tmat,tproj,tinv);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Coarsen a SINGLE RHS fine operator.
  //
  // No multiRHS packing: the operator is applied once per phased basis
  // vector and the batch is assembled on the coarse side by the mixed
  // blockProject, which takes a vector of fine fields and writes the batch
  // coarse field the accumulate expects. Only nbv applications per momentum,
  // so a batch that does not divide nbasis wastes nothing, and the live fine
  // storage is batch fields rather than two D+1 fields of extent batch.
  ///////////////////////////////////////////////////////////////////////////
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       std::vector<FineField> &Subspace,
		       GridBase *CoarseGrid,
		       int batch)
  {
    RealD tproj=0.0, tmat=0.0, tphase=0.0, tphaseBZ=0.0, tslice=0.0, tinv=0.0;

    std::cout << GridLogMessage<< "GeneralCoarsenMatrixMrhs (single RHS fine operator)"<< std::endl;

    GRID_ASSERT(Subspace.size()==nbasis);
    GRID_ASSERT(batch>=1);
    GridBase *grid = Subspace[0].Grid();

    Coordinate blatt,bsimd,bmpi;
    CoarsenBlockGridLayout(grid,CoarseGrid,blatt,bsimd,bmpi);
    GridCartesian BlockGrid(blatt,bsimd,bmpi);

    BlockComplexField InnerProd(&BlockGrid);
    blockOrthogonalise(InnerProd,Subspace);

    MultiRHSBlockProject<Lattice<Fobj> >    Projector;
    Projector.Allocate(nbasis,grid,CoarseGrid);
    Projector.ImportBasis(Subspace);

    const int npoint = geom_srhs.npoint;

    Eigen::MatrixXcd invMkl;
    CoarsenFourierMatrix(CoarseGrid,invMkl);

    FineField phaV(grid);
    std::vector<FineComplexField>   phaF(npoint,grid);
    std::vector<CoarseComplexField> pha (npoint,CoarseGrid);

    tphase=-usecond();
    CoarsenPhases(grid,CoarseGrid,&BlockGrid,pha,phaF);
    tphase+=usecond();

    std::vector<CoarseMatrix> _A;
    _A.resize(npoint,CoarseGrid);
    for(int k=0;k<npoint;k++) _A[k] = Zero();

    Coordinate cmlatt,cmsimd,cmmpi;
    CoarsenBatchGridLayout(CoarseGrid,batch,cmlatt,cmsimd,cmmpi);
    GridCartesian CoarseBatchGrid(cmlatt,cmsimd,cmmpi);

    CoarseVector       TmpProj(&CoarseBatchGrid);
    CoarseComplexField phaB(&CoarseBatchGrid);

    std::vector<FineField> MphaV(batch,grid);

    for(int i0=0;i0<nbasis;i0+=batch){ // Loop over batches of basis vectors

      int nbv = MIN(batch,nbasis-i0);
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i0<<"/"<<nbasis<< std::endl;

      for(int p=0;p<npoint;p++){ // Loop over momenta

	for(int b=0;b<nbv;b++){
	  tphaseBZ-=usecond();
	  phaV = phaF[p]*Subspace[i0+b];
	  tphaseBZ+=usecond();
	  tmat-=usecond();
	  linop.Op(phaV,MphaV[b]);
	  tmat+=usecond();
	}
	// The accumulate reads only the first nbv slices, but the projector
	// sees the whole vector, so the tail must not be undefined.
	for(int b=nbv;b<batch;b++) MphaV[b] = Zero();

	tproj-=usecond();
	Projector.blockProject(MphaV,TmpProj);
	tproj+=usecond();

	tinv-=usecond();
	CoarsenAccumulate(p,i0,nbv,batch,invMkl,pha,phaB,TmpProj,_A,CoarseGrid);
	tinv+=usecond();
      }
    }

    for(int p=0;p<npoint;p++){
      GridtoBLAS(_A[p],BLAS_A[p]);
    }
    CoarsenReport(tphase,tphaseBZ,tslice,tmat,tproj,tinv);
  }
  void Mdag(const CoarseVector &in, CoarseVector &out)
  {
    this->M(in,out);
  }
  void M (const CoarseVector &in, CoarseVector &out)
  {
    //    std::cout << GridLogMessage << "New Mrhs coarse"<<std::endl;
    conformable(CoarseGrid(),in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();

    RealD t_tot;
    RealD t_exch;
    RealD t_GtoB;
    RealD t_BtoG;
    RealD t_mult;

    CheckGridSet();
    if ( in.Grid() != _CoarseGridMulti ) {
      std::cout << GridLogError
		<< "MultiGeneralCoarsenedOperatorV2::M called with a field on a"
		<< std::endl;
      std::cout << GridLogError
		<< "  different grid object from the one given to SetGrid(). Two"
		<< std::endl;
      std::cout << GridLogError
		<< "  grids of identical shape do not conform; share one object."
		<< std::endl;
      GRID_ASSERT(in.Grid() == _CoarseGridMulti);
    }

    GRID_TRACE("CoarseV2Mult");
    t_tot=-usecond();
    CoarseVector tin=in;
    t_exch=-usecond();
    // lambda scope so the roctx range covers exactly the exchange; the
    // PaddedCellFwd/BwdMPI markers inside it then nest properly.
    CoarseVector pin = [&](){ GRID_TRACE("CoarseV2Exchange");
                              return CellMulti->ExchangePeriodic(tin); }(); //padded input
    t_exch+=usecond();

    CoarseVector pout(pin.Grid());

    int npoint = geom.npoint;
    typedef calcMatrix* Aview;
    typedef LatticeView<Cvec> Vview;
      
    const int Nsimd = CComplex::Nsimd();

    int64_t nrhs  =pin.Grid()->GlobalDimensions()[0];
    GRID_ASSERT(nrhs>=1);

    RealD flops,bytes;
    int64_t osites=in.Grid()->oSites(); // unpadded
    int64_t unpadded_vol = CoarseGrid()->lSites()/nrhs;
    
    flops = 1.0* npoint * nbasis * nbasis * 8.0 * osites * CComplex::Nsimd();
    bytes = 1.0*osites*sizeof(siteMatrix)*npoint/pin.Grid()->GlobalDimensions()[0]
          + 2.0*osites*sizeof(siteVector)*npoint;
    

    t_GtoB=-usecond();
    { GRID_TRACE("CoarseV2GridToBLAS");
      GridtoBLAS(pin,BLAS_B);
    }
    t_GtoB+=usecond();

    GridBLAS BLAS;

    t_mult=-usecond();
    { GRID_TRACE("CoarseV2StencilGEMM");
    for(int p=0;p<geom.npoint;p++){
      RealD c = 1.0;
      if (p==0) c = 0.0;
      ComplexD beta(c);

      BLAS.gemmBatched(nbasis,nrhs,nbasis,
		       ComplexD(1.0),
		       BLAS_AP[p], 
		       BLAS_BP[p], 
		       ComplexD(c), 
		       BLAS_CP);
    }
    BLAS.synchronise();
    }
    t_mult+=usecond();

    t_BtoG=-usecond();
    { GRID_TRACE("CoarseV2BLASToGrid");
      BLAStoGrid(out,BLAS_C);
    }
    t_BtoG+=usecond();
    t_tot+=usecond();
    /*
    std::cout << GridLogMessage << "New Mrhs coarse DONE "<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult exch "<<t_exch<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult mult "<<t_mult<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult GtoB  "<<t_GtoB<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult BtoG  "<<t_BtoG<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult tot  "<<t_tot<<" us"<<std::endl;
    */
    //    std::cout << GridLogMessage<<std::endl;
    //    std::cout << GridLogMessage<<"Coarse Kernel flops "<< flops<<std::endl;
    //    std::cout << GridLogMessage<<"Coarse Kernel flop/s "<< flops/t_mult<<" mflop/s"<<std::endl;
    //    std::cout << GridLogMessage<<"Coarse Kernel bytes/s "<< bytes/t_mult/1000<<" GB/s"<<std::endl;
    //    std::cout << GridLogMessage<<"Coarse overall flops/s "<< flops/t_tot<<" mflop/s"<<std::endl;
    //    std::cout << GridLogMessage<<"Coarse total bytes   "<< bytes/1e6<<" MB"<<std::endl;
  };
  virtual  void Mdiag    (const Field &in, Field &out){ GRID_ASSERT(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};
};
  
NAMESPACE_END(Grid);
