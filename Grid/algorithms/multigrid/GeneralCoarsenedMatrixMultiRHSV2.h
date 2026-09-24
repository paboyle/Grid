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
  // The BLAS scalar of this level follows the coefficient precision:
  // ComplexF for an fp32 coarse space, ComplexD for fp64.
  typedef typename GridTypeMapper<CComplex>::scalar_type CoarseBLASScalar;
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

  std::vector<deviceVector<CoarseBLASScalar *> > BLAS_AP;
  std::vector<deviceVector<CoarseBLASScalar *> > BLAS_BP;
  deviceVector<CoarseBLASScalar *>               BLAS_CP;

  ///////////////////////////////////////////////////////////////////////////
  // Stencil legs carried in the BATCH dimension, LegGroup at a time.
  //
  // One call per stencil point batches over the local coarse volume alone,
  // which at the production point is 1024 sites.  That is too few for the
  // single-precision kernel: it then runs at half the bandwidth the double
  // one reaches and takes the same time, so an fp32 coarse space buys
  // nothing in the operator.  The same shape at four times the batch runs
  // 1.7x faster in fp32 than fp64, so the fix is more batch, not fewer
  // bytes.  Grouping G legs multiplies the batch by G and divides the launch
  // count by G, at the price of G partial outputs summed at the end.
  //
  // G must divide npoint; 3 divides both the 33 point (next-to-nearest) and
  // 81 point (full box) stencils.  G=1 is the ungrouped form.  The only
  // numerical difference is the ORDER of the sum over stencil points.
  ///////////////////////////////////////////////////////////////////////////
  // The group need not divide npoint: the last call carries the remainder,
  // with its own smaller tables.  The batch a call presents is LegGroup times
  // the local coarse volume, so the group that saturates the kernel depends
  // on the decomposition, and the default is chosen from it rather than
  // fixed.  TargetBatch is where the measured single-precision curve flattens
  // (530 GB/s at 1024, 695 at 9216, 728 at 33792 for 60x12x60 on MI250X).
  //
  // 9216 is nine times the production local coarse volume, so the default
  // lands on nine legs per call: 81 = 9x9 for the full-box stencil and
  // 33 = 9+9+9+6 for the next-to-nearest one, which is the grouping the
  // earlier HDCG coarse operator used.
  static const int TargetBatch = 9216;

  int LegGroup = 1;
  std::vector<deviceVector<CoarseBLASScalar *> > BLAS_APg; // per group
  std::vector<deviceVector<CoarseBLASScalar *> > BLAS_BPg;
  deviceVector<CoarseBLASScalar *>               BLAS_CPg; // LegGroup*sites
  deviceVector<CoarseBLASScalar *>               BLAS_CPr; // remainder*sites
  deviceVector<calcVector>                       BLAS_Cg;  // LegGroup partials

  int LegGroups(void)    const { return geom.npoint/LegGroup; }          // full
  int LegRemainder(void) const { return geom.npoint%LegGroup; }

  ///////////////////////////////////////////////////////////////////////////
  // Matrix pointer table for the grouped call.  Group g holds legs
  // g*LegGroup .. and the batch index runs site fastest within a leg.  The
  // last group is short when LegGroup does not divide npoint.
  ///////////////////////////////////////////////////////////////////////////
  void BuildGroupedA(void)
  {
    int32_t sites = _CoarseGrid->lSites();
    int ng = LegGroups() + (LegRemainder()?1:0);
    BLAS_APg.resize(ng);
    for(int g=0;g<ng;g++){
      int legs = (g<LegGroups()) ? LegGroup : LegRemainder();
      BLAS_APg[g].resize(legs*sites);
      for(int l=0;l<legs;l++){
	int p = g*LegGroup+l;
	for(int32_t ss=0;ss<sites;ss++){
	  CoarseBLASScalar *ptr = (CoarseBLASScalar *)&BLAS_A[p][ss];
	  acceleratorPut(BLAS_APg[g][l*sites+ss],ptr);
	}
      }
    }
  }

  void SetLegGroup(int G)
  {
    GRID_ASSERT( G>=1 );
    GRID_ASSERT( G<=geom.npoint );   // every partial must get an initialising call
    LegGroup = G;
    BuildGroupedA();
    if ( _CoarseGridMulti ) SetGrid(_CoarseGridMulti);   // rebuild the rest
  }

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
  // Resident device memory this operator holds (matrix elements and the
  // Nrhs-dependent BLAS buffers).  Not evictable.
  uint64_t DeviceBytes(void)
  {
    uint64_t b = (uint64_t)(BLAS_B.capacity()+BLAS_C.capacity()+BLAS_Cg.capacity())*sizeof(calcVector);
    for(int p=0;p<(int)BLAS_A.size();p++) b += (uint64_t)BLAS_A[p].capacity()*sizeof(calcMatrix);
    return b;
  }

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
	CoarseBLASScalar *ptr = (CoarseBLASScalar *)&BLAS_A[p][ss];
	acceleratorPut(BLAS_AP[p][ss],ptr);
      }
    }

    // Enough legs to bring the batch up to where the kernel saturates, and
    // never more legs than the stencil has.
    LegGroup = (TargetBatch + unpadded_sites - 1)/unpadded_sites;
    if ( LegGroup < 1 )            LegGroup = 1;
    if ( LegGroup > geom.npoint )  LegGroup = geom.npoint;
    BuildGroupedA();
    std::cout << GridLogMessage << "MultiGeneralCoarsenedOperatorV2: stencil "
	      << geom.npoint << " points, local coarse volume " << unpadded_sites
	      << ", legs per GEMM " << LegGroup
	      << " (batch " << LegGroup*unpadded_sites << ", "
	      << LegGroups() << " full call(s)"
	      << (LegRemainder() ? " + a remainder of "+std::to_string(LegRemainder()) : "")
	      << ")" << std::endl;
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
    BLAS_Cg.resize(0);
    BLAS_CPg.resize(0);
    BLAS_CPr.resize(0);
    for(int g=0;g<(int)BLAS_BPg.size();g++){ BLAS_BPg[g].resize(0); }
    BLAS_BPg.resize(0);
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
      CoarseBLASScalar *ptr = (CoarseBLASScalar *)&BLAS_C[ss*nrhs];
      acceleratorPut(BLAS_CP[ss],ptr);
    }

    // Grouped output: LegGroup partial results, each the full C volume, and
    // their pointer table.  The batch index runs site fastest within a leg,
    // matching the grouped A and B tables.  The remainder call writes into
    // the first few partials, so it needs a truncated copy of the table.
    int nfull = LegGroups();
    int rem   = LegRemainder();
    int ng    = nfull + (rem?1:0);
    BLAS_Cg.resize (LegGroup*nrhs*unpadded_sites);
    BLAS_CPg.resize(LegGroup*unpadded_sites);
    BLAS_CPr.resize(rem*unpadded_sites);
    for(int l=0;l<LegGroup;l++){
      for(int ss=0;ss<unpadded_sites;ss++){
	CoarseBLASScalar *ptr = (CoarseBLASScalar *)&BLAS_Cg[(l*unpadded_sites+ss)*nrhs];
	acceleratorPut(BLAS_CPg[l*unpadded_sites+ss],ptr);
	if ( l<rem ) acceleratorPut(BLAS_CPr[l*unpadded_sites+ss],ptr);
      }
    }
    BLAS_BPg.resize(ng);
    for(int g=0;g<ng;g++){
      int legs = (g<nfull) ? LegGroup : rem;
      BLAS_BPg[g].resize(legs*unpadded_sites);
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
	  CoarseBLASScalar * ptr = (CoarseBLASScalar *)&BLAS_B[nbr];
	  acceleratorPut(BLAS_BP[point][j],ptr); // neighbour indexing in ghost zone volume
	  acceleratorPut(BLAS_BPg[point/LegGroup][(point%LegGroup)*unpadded_sites+j],ptr);
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
  ///////////////////////////////////////////////////////////////////////////
  // Probe momenta.  The stencil DISPLACEMENTS are the +-1 box; the momenta
  // used to separate them are free, and are chosen here to span the
  // Brillouin zone: k_mu = K_mu * shift_mu with K_mu ~ L_mu/3, so a phase
  // is ~2pi/3 per unit displacement rather than 2pi/L_mu.
  //
  // This is what conditions the extraction.  With K_mu = 1 every entry of
  // the phase matrix tends to 1 as the coarse lattice grows, the matrix
  // tends to rank one, and its inverse amplifies any error in the measured
  // projections: at 24.24.16.32 the condition number is 2.7e5, so fp32
  // projections give coarse matrix elements wrong by ~1e-3.  Spread momenta
  // bring it to ~20, independent of the lattice size.
  //
  // K_mu is stepped away from L_mu/2, where +k and -k alias into the same
  // momentum and the matrix is singular.
  ///////////////////////////////////////////////////////////////////////////
  void CoarsenMomenta(GridBase *CoarseGrid,Coordinate &K)
  {
    Coordinate clatt = CoarseGrid->GlobalDimensions();
    int Nd = CoarseGrid->Nd();
    K.resize(Nd);
    for(int mu=0;mu<Nd;mu++){
      int L = clatt[mu];
      int k = (L>=3) ? (int)std::lround(L/3.0) : 1;
      if ( k < 1 ) k = 1;
      if ( (L>2) && ((2*k)%L == 0) ) k = k-1;     // +k and -k must differ
      if ( k < 1 ) k = 1;
      K[mu] = k;
    }
  }

  void CoarsenFourierMatrix(GridBase *CoarseGrid,Eigen::MatrixXcd &invMkl)
  {
    const int npoint = geom_srhs.npoint;
    Coordinate clatt = CoarseGrid->GlobalDimensions();
    int Nd = CoarseGrid->Nd();
    Coordinate K;  CoarsenMomenta(CoarseGrid,K);

    Eigen::MatrixXcd Mkl = Eigen::MatrixXcd::Zero(npoint,npoint);
    ComplexD ci(0.0,1.0);
    for(int k=0;k<npoint;k++){ // Loop over momenta
      for(int l=0;l<npoint;l++){ // Loop over nbr relative
	ComplexD phase(0.0,0.0);
	for(int mu=0;mu<Nd;mu++){
	  RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	  phase=phase+TwoPiL*K[mu]*geom_srhs.shifts[k][mu]*geom_srhs.shifts[l][mu];
	}
	phase=exp(phase*ci);
	Mkl(k,l) = phase;
      }
    }
    invMkl = Mkl.inverse();

    // The extraction's error amplification, printed because it is the thing
    // that decides whether an fp32 coarse sector is usable.
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(Mkl);
    RealD cond = svd.singularValues()(0)/svd.singularValues()(npoint-1);
    std::cout << GridLogMessage << "CoarsenOperator: probe momenta "<<K
	      <<" on coarse lattice "<<clatt
	      <<" ; Fourier matrix condition number "<<cond<<std::endl;
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
    Coordinate K;  CoarsenMomenta(CoarseGrid,K);
    ComplexD ci(0.0,1.0);

    // The fine-side scratch carries the FINE level's precision, which need
    // not be the coarse one (fp64 fine, fp32 coarse at L1).
    typedef typename GridTypeMapper<FineInner>::scalar_type FineScalar;
    FineComplexField one(grid); one=FineScalar(1.0);
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
	pha[p]  = pha[p]  + (TwoPiL * K[mu] * geom_srhs.shifts[p][mu]) * coor;
	LatticeCoordinate(blk_coor,mu);
	pha_blk = pha_blk + (TwoPiL * K[mu] * geom_srhs.shifts[p][mu]) * blk_coor;
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
  // Owning the transfer operator: the caller passes one, this imports the
  // orthonormalised basis into it, and the caller keeps it for the solve.
  // The overload below makes its own and is kept for the pre-2026 drivers;
  // it costs a second full basis store, which is why nothing new should use
  // it (documentation/MultiGridMemoryAudit.md).
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       GridCartesian *FineGridMulti,
		       std::vector<FineField> &Subspace,
		       GridBase *CoarseGrid)
  {
    MultiRHSBlockProject<Lattice<Fobj> > Projector;
    CoarsenOperator(linop,FineGridMulti,Subspace,CoarseGrid,Projector);
  }
  template<class Projector_t>
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       GridCartesian *FineGridMulti,
		       std::vector<FineField> &Subspace,
		       GridBase *CoarseGrid,
		       Projector_t &Projector)
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

    // the caller owns it; import the basis we have just orthonormalised
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
    MultiRHSBlockProject<Lattice<Fobj> > Projector;
    CoarsenOperator(linop,Subspace,CoarseGrid,batch,Projector);
  }
  template<class Projector_t>
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       std::vector<FineField> &Subspace,
		       GridBase *CoarseGrid,
		       int batch,
		       Projector_t &Projector)
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

    // the caller owns it; import the basis we have just orthonormalised
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
    t_exch=-usecond();
    // The exchange takes its input by const reference, so it reads the
    // caller's field directly.  lambda scope so the roctx range covers
    // exactly the exchange; the PaddedCellFwd/BwdMPI markers inside it then
    // nest properly.
    CoarseVector pin = [&](){ GRID_TRACE("CoarseV2Exchange");
                              return CellMulti->ExchangePeriodic(in); }(); //padded input
    t_exch+=usecond();

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
    // The scalar type selects the GEMM: Cgemm for an fp32 coarse space,
    // Zgemm for fp64.
    if ( LegGroup == 1 ) {
      for(int p=0;p<geom.npoint;p++){
	RealD c = (p==0) ? 0.0 : 1.0;
	BLAS.gemmBatched(nbasis,nrhs,nbasis,
			 CoarseBLASScalar(1.0),
			 BLAS_AP[p],
			 BLAS_BP[p],
			 CoarseBLASScalar(c),
			 BLAS_CP);
      }
    } else {
      // LegGroup legs per call: the batch is LegGroup times the local volume
      // and there are LegGroup partial outputs, accumulated across the
      // groups and summed below.  A final short call carries the remainder
      // when the group does not divide the stencil; it writes into the first
      // few partials, which the full calls have already initialised.
      int nfull = LegGroups();
      int rem   = LegRemainder();
      for(int g=0;g<nfull;g++){
	RealD c = (g==0) ? 0.0 : 1.0;
	BLAS.gemmBatched(nbasis,nrhs,nbasis,
			 CoarseBLASScalar(1.0),
			 BLAS_APg[g],
			 BLAS_BPg[g],
			 CoarseBLASScalar(c),
			 BLAS_CPg);
      }
      if ( rem ) {
	RealD c = (nfull==0) ? 0.0 : 1.0;
	BLAS.gemmBatched(nbasis,nrhs,nbasis,
			 CoarseBLASScalar(1.0),
			 BLAS_APg[nfull],
			 BLAS_BPg[nfull],
			 CoarseBLASScalar(c),
			 BLAS_CPr);
      }
    }
    BLAS.synchronise();
    }
    t_mult+=usecond();

    if ( LegGroup > 1 ) { GRID_TRACE("CoarseV2LegSum");
      // Sum the LegGroup partial results into the single output buffer.
      // Flat over scalars: a calcVector is nbasis of them and the partials
      // are contiguous, one whole C volume after another.
      int64_t nscalar = (int64_t)BLAS_C.size()*nbasis;   // one whole C volume
      const int G = LegGroup;
      CoarseBLASScalar *dst = (CoarseBLASScalar *)&BLAS_C[0];
      CoarseBLASScalar *src = (CoarseBLASScalar *)&BLAS_Cg[0];
      accelerator_for(i,nscalar,1,{
	  CoarseBLASScalar sum = src[i];
	  for(int l=1;l<G;l++) sum = sum + src[(int64_t)l*nscalar + i];
	  dst[i] = sum;
	});
    }

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
