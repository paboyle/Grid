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
class MultiGeneralCoarsenedMatrix : public SparseMatrixBase<Lattice<iVector<CComplex,nbasis > > >  {
public:
  typedef typename CComplex::scalar_object SComplex;
  typedef GeneralCoarsenedMatrix<Fobj,CComplex,nbasis> GeneralCoarseOp;
  typedef MultiGeneralCoarsenedMatrix<Fobj,CComplex,nbasis> MultiGeneralCoarseOp;

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
  typedef Lattice<CComplex >    FineComplexField;
  typedef CoarseVector Field;

  ////////////////////
  // Data members
  ////////////////////
  GridCartesian *       _CoarseGridMulti; 
  NonLocalStencilGeometry geom;
  NonLocalStencilGeometry geom_srhs;
  PaddedCell Cell;
  GeneralLocalStencil Stencil;

  deviceVector<calcVector> BLAS_B;
  deviceVector<calcVector> BLAS_C;
  std::vector<deviceVector<calcMatrix> > BLAS_A;

  std::vector<deviceVector<ComplexD *> > BLAS_AP;
  std::vector<deviceVector<ComplexD *> > BLAS_BP;
  deviceVector<ComplexD *>               BLAS_CP;

  ///////////////////////
  // Interface
  ///////////////////////
  GridBase      * Grid(void)           { return _CoarseGridMulti; };   // this is all the linalg routines need to know
  GridCartesian * CoarseGrid(void)     { return _CoarseGridMulti; };   // this is all the linalg routines need to know

  // Can be used to do I/O on the operator matrices externally
  void SetMatrix (int p,CoarseMatrix & A)
  {
    assert(A.size()==geom_srhs.npoint);
    GridtoBLAS(A[p],BLAS_A[p]);
  }
  void GetMatrix (int p,CoarseMatrix & A)
  {
    assert(A.size()==geom_srhs.npoint);
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
  
  MultiGeneralCoarsenedMatrix(NonLocalStencilGeometry &_geom,GridCartesian *CoarseGridMulti) :
    _CoarseGridMulti(CoarseGridMulti),
    geom_srhs(_geom),
    geom(_CoarseGridMulti,_geom.hops,_geom.skip+1),
    Cell(geom.Depth(),_CoarseGridMulti),
    Stencil(Cell.grids.back(),geom.shifts) // padded cell stencil
  {
    int32_t padded_sites   = Cell.grids.back()->lSites();
    int32_t unpadded_sites = CoarseGridMulti->lSites();
    
    int32_t nrhs  = CoarseGridMulti->FullDimensions()[0];  // # RHS
    int32_t orhs  = nrhs/CComplex::Nsimd();

    padded_sites   = padded_sites/nrhs;
    unpadded_sites = unpadded_sites/nrhs;
    
    /////////////////////////////////////////////////
    // Device data vector storage
    /////////////////////////////////////////////////
    BLAS_A.resize(geom.npoint);
    for(int p=0;p<geom.npoint;p++){
      BLAS_A[p].resize (unpadded_sites); // no ghost zone, npoint elements
    }
    
    BLAS_B.resize(nrhs *padded_sites);   // includes ghost zone
    BLAS_C.resize(nrhs *unpadded_sites); // no ghost zone
    BLAS_AP.resize(geom.npoint);
    BLAS_BP.resize(geom.npoint);
    for(int p=0;p<geom.npoint;p++){
      BLAS_AP[p].resize(unpadded_sites);
      BLAS_BP[p].resize(unpadded_sites);
    }
    BLAS_CP.resize(unpadded_sites);

    /////////////////////////////////////////////////
    // Pointers to data
    /////////////////////////////////////////////////

    // Site identity mapping for A
    for(int p=0;p<geom.npoint;p++){
      for(int ss=0;ss<unpadded_sites;ss++){
	ComplexD *ptr = (ComplexD *)&BLAS_A[p][ss];
	acceleratorPut(BLAS_AP[p][ss],ptr);
      }
    }
    // Site identity mapping for C
    for(int ss=0;ss<unpadded_sites;ss++){
      ComplexD *ptr = (ComplexD *)&BLAS_C[ss*nrhs];
      acceleratorPut(BLAS_CP[ss],ptr);
    }

    // Neighbour table is more complicated
    int32_t j=0; // Interior point counter (unpadded)
    for(int32_t s=0;s<padded_sites;s++){ // 4 volume, padded
      int ghost_zone=0;
      for(int32_t point = 0 ; point < geom.npoint; point++){
	int i=s*orhs*geom.npoint+point;
	if( Stencil._entries[i]._wrap ) { // stencil is indexed by the oSite of the CoarseGridMulti, hence orhs factor
	  ghost_zone=1; // If general stencil wrapped in any direction, wrap=1
	}
      }

      if( ghost_zone==0) {
	for(int32_t point = 0 ; point < geom.npoint; point++){
	  int i=s*orhs*geom.npoint+point;
 	  int32_t nbr = Stencil._entries[i]._offset*CComplex::Nsimd(); // oSite -> lSite
	  assert(nbr<BLAS_B.size());
	  ComplexD * ptr = (ComplexD *)&BLAS_B[nbr];
	  acceleratorPut(BLAS_BP[point][j],ptr); // neighbour indexing in ghost zone volume
	}
	j++;
      }
    }
    assert(j==unpadded_sites);
  }
  template<class vobj> void GridtoBLAS(const Lattice<vobj> &from,deviceVector<typename vobj::scalar_object> &to)
  {
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  GridBase *Fg = from.Grid();
  assert(!Fg->_isCheckerBoarded);
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
  assert(!Tg->_isCheckerBoarded);
  int nd = Tg->_ndimension;
  
  assert(in.size()==Tg->lSites());

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
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & Subspace,
		       GridBase *CoarseGrid)
  {
#if 0
    std::cout << GridLogMessage<< "GeneralCoarsenMatrixMrhs "<< std::endl;

    GridBase *grid = Subspace.FineGrid;

    /////////////////////////////////////////////////////////////
    // Orthogonalise the subblocks over the basis
    /////////////////////////////////////////////////////////////
    CoarseScalar InnerProd(CoarseGrid); 
    blockOrthogonalise(InnerProd,Subspace.subspace);

    const int npoint = geom_srhs.npoint;

    Coordinate clatt = CoarseGrid->GlobalDimensions();
    int Nd = CoarseGrid->Nd();
      /*
       *     Here, k,l index which possible momentum/shift within the N-points connected by MdagM.
       *     Matrix index i is mapped to this shift via 
       *               geom.shifts[i]
       *
       *     conj(pha[block]) proj[k (which mom)][j (basis vec cpt)][block] 
       *       =  \sum_{l in ball}  e^{i q_k . delta_l} < phi_{block,j} | MdagM | phi_{(block+delta_l),i} > 
       *       =  \sum_{l in ball} e^{iqk.delta_l} A_ji^{b.b+l}
       *       = M_{kl} A_ji^{b.b+l}
       *
       *     Must assemble and invert matrix M_k,l = e^[i q_k . delta_l]
       *  
       *     Where q_k = delta_k . (2*M_PI/global_nb[mu])
       *
       *     Then A{ji}^{b,b+l} = M^{-1}_{lm} ComputeProj_{m,b,i,j}
       */
    Eigen::MatrixXcd Mkl    = Eigen::MatrixXcd::Zero(npoint,npoint);
    Eigen::MatrixXcd invMkl = Eigen::MatrixXcd::Zero(npoint,npoint);
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

    ///////////////////////////////////////////////////////////////////////
    // Now compute the matrix elements of linop between the orthonormal
    // set of vectors.
    ///////////////////////////////////////////////////////////////////////
    FineField phaV(grid); // Phased block basis vector
    FineField MphaV(grid);// Matrix applied
    std::vector<FineComplexField> phaF(npoint,grid);
    std::vector<CoarseComplexField> pha(npoint,CoarseGrid);
    
    CoarseVector coarseInner(CoarseGrid);
    
    typedef typename CComplex::scalar_type SComplex;
    FineComplexField one(grid); one=SComplex(1.0);
    FineComplexField zz(grid); zz = Zero();
    for(int p=0;p<npoint;p++){ // Loop over momenta in npoint
      /////////////////////////////////////////////////////
      // Stick a phase on every block
      /////////////////////////////////////////////////////
      CoarseComplexField coor(CoarseGrid);
      pha[p]=Zero();
      for(int mu=0;mu<Nd;mu++){
	LatticeCoordinate(coor,mu);
	RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	pha[p] = pha[p] + (TwoPiL * geom_srhs.shifts[p][mu]) * coor;
      }
      pha[p]  =exp(pha[p]*ci);	

      blockZAXPY(phaF[p],pha[p],one,zz);
    }

    // Could save on temporary storage here
    std::vector<CoarseMatrix> _A;
    _A.resize(geom_srhs.npoint,CoarseGrid);

    std::vector<CoarseVector> ComputeProj(npoint,CoarseGrid);
    CoarseVector          FT(CoarseGrid);
    for(int i=0;i<nbasis;i++){// Loop over basis vectors
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i<<"/"<<nbasis<< std::endl;
      for(int p=0;p<npoint;p++){ // Loop over momenta in npoint

	phaV = phaF[p]*Subspace.subspace[i];

	/////////////////////////////////////////////////////////////////////
	// Multiple phased subspace vector by matrix and project to subspace
	// Remove local bulk phase to leave relative phases
	/////////////////////////////////////////////////////////////////////
	linop.Op(phaV,MphaV);

	// Fixme, could use batched block projector here
	blockProject(coarseInner,MphaV,Subspace.subspace);

	coarseInner = conjugate(pha[p]) * coarseInner;

	ComputeProj[p] = coarseInner;
      }

      // Could do this with a block promote or similar BLAS call via the MultiRHSBlockProjector with a const matrix.
      for(int k=0;k<npoint;k++){

	FT = Zero();
	for(int l=0;l<npoint;l++){
	  FT= FT+ invMkl(l,k)*ComputeProj[l];
	}
      
	int osites=CoarseGrid->oSites();
	autoView( A_v  , _A[k], AcceleratorWrite);
	autoView( FT_v  , FT, AcceleratorRead);
	accelerator_for(sss, osites, 1, {
	    for(int j=0;j<nbasis;j++){
	      A_v[sss](i,j) = FT_v[sss](j);
	    }
        });
      }
    }

    // Only needed if nonhermitian
    //    if ( ! hermitian ) {
    //      std::cout << GridLogMessage<<"PopulateAdag  "<<std::endl;
    //      PopulateAdag();
    //    }
    // Need to write something to populate Adag from A

    for(int p=0;p<geom_srhs.npoint;p++){
      GridtoBLAS(_A[p],BLAS_A[p]);
    }
    /*
Grid : Message : 11698.730546 s : CoarsenOperator eigen  1334 us
Grid : Message : 11698.730563 s : CoarsenOperator phase  34729 us
Grid : Message : 11698.730565 s : CoarsenOperator phaseBZ 2423814 us
Grid : Message : 11698.730566 s : CoarsenOperator mat    127890998 us
Grid : Message : 11698.730567 s : CoarsenOperator proj   515840840 us
Grid : Message : 11698.730568 s : CoarsenOperator inv    103948313 us
Takes 600s to compute matrix elements, DOMINATED by the block project.
Easy to speed up with the batched block project.
Store npoint vectors, get npoint x Nbasis block projection, and 81 fold faster.

// Block project below taks to 240s
Grid : Message : 328.193418 s : CoarsenOperator phase      38338 us
Grid : Message : 328.193434 s : CoarsenOperator phaseBZ  1711226 us
Grid : Message : 328.193436 s : CoarsenOperator mat    122213270 us
//Grid : Message : 328.193438 s : CoarsenOperator proj   1181154 us <-- this is mistimed
//Grid : Message : 11698.730568 s : CoarsenOperator inv  103948313 us <-- Cut this ~10x if lucky by loop fusion
     */
#else
    RealD tproj=0.0;
    RealD tmat=0.0;
    RealD tphase=0.0;
    RealD tphaseBZ=0.0;
    RealD tinv=0.0;

    std::cout << GridLogMessage<< "GeneralCoarsenMatrixMrhs "<< std::endl;

    GridBase *grid = Subspace.FineGrid;

    /////////////////////////////////////////////////////////////
    // Orthogonalise the subblocks over the basis
    /////////////////////////////////////////////////////////////
    CoarseScalar InnerProd(CoarseGrid); 
    blockOrthogonalise(InnerProd,Subspace.subspace);


    MultiRHSBlockProject<Lattice<Fobj> >    Projector;
    Projector.Allocate(nbasis,grid,CoarseGrid);
    Projector.ImportBasis(Subspace.subspace);
    
    const int npoint = geom_srhs.npoint;

    Coordinate clatt = CoarseGrid->GlobalDimensions();
    int Nd = CoarseGrid->Nd();
      /*
       *     Here, k,l index which possible momentum/shift within the N-points connected by MdagM.
       *     Matrix index i is mapped to this shift via 
       *               geom.shifts[i]
       *
       *     conj(pha[block]) proj[k (which mom)][j (basis vec cpt)][block] 
       *       =  \sum_{l in ball}  e^{i q_k . delta_l} < phi_{block,j} | MdagM | phi_{(block+delta_l),i} > 
       *       =  \sum_{l in ball} e^{iqk.delta_l} A_ji^{b.b+l}
       *       = M_{kl} A_ji^{b.b+l}
       *
       *     Must assemble and invert matrix M_k,l = e^[i q_k . delta_l]
       *  
       *     Where q_k = delta_k . (2*M_PI/global_nb[mu])
       *
       *     Then A{ji}^{b,b+l} = M^{-1}_{lm} ComputeProj_{m,b,i,j}
       */
    Eigen::MatrixXcd Mkl    = Eigen::MatrixXcd::Zero(npoint,npoint);
    Eigen::MatrixXcd invMkl = Eigen::MatrixXcd::Zero(npoint,npoint);
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

    ///////////////////////////////////////////////////////////////////////
    // Now compute the matrix elements of linop between the orthonormal
    // set of vectors.
    ///////////////////////////////////////////////////////////////////////
    FineField phaV(grid); // Phased block basis vector
    FineField MphaV(grid);// Matrix applied
    std::vector<FineComplexField> phaF(npoint,grid);
    std::vector<CoarseComplexField> pha(npoint,CoarseGrid);
    
    CoarseVector coarseInner(CoarseGrid);
    
    tphase=-usecond();
    typedef typename CComplex::scalar_type SComplex;
    FineComplexField one(grid); one=SComplex(1.0);
    FineComplexField zz(grid); zz = Zero();
    for(int p=0;p<npoint;p++){ // Loop over momenta in npoint
      /////////////////////////////////////////////////////
      // Stick a phase on every block
      /////////////////////////////////////////////////////
      CoarseComplexField coor(CoarseGrid);
      pha[p]=Zero();
      for(int mu=0;mu<Nd;mu++){
	LatticeCoordinate(coor,mu);
	RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	pha[p] = pha[p] + (TwoPiL * geom_srhs.shifts[p][mu]) * coor;
      }
      pha[p]  =exp(pha[p]*ci);	

      blockZAXPY(phaF[p],pha[p],one,zz);
    }
    tphase+=usecond();

    // Could save on temporary storage here
    std::vector<CoarseMatrix> _A;
    _A.resize(geom_srhs.npoint,CoarseGrid);

    // Count use small chunks than npoint == 81 and save memory
    int batch = 9;
    std::vector<FineField>    _MphaV(batch,grid);
    std::vector<CoarseVector> TmpProj(batch,CoarseGrid);

    std::vector<CoarseVector> ComputeProj(npoint,CoarseGrid);
    CoarseVector          FT(CoarseGrid);
    for(int i=0;i<nbasis;i++){// Loop over basis vectors
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i<<"/"<<nbasis<< std::endl;

      //      std::cout << GridLogMessage << " phasing the fine vector "<<std::endl;
      // Fixme : do this in batches
      for(int p=0;p<npoint;p+=batch){ // Loop over momenta in npoint

	for(int b=0;b<MIN(batch,npoint-p);b++){
	  tphaseBZ-=usecond();
	  phaV = phaF[p+b]*Subspace.subspace[i];
	  tphaseBZ+=usecond();

	  /////////////////////////////////////////////////////////////////////
	  // Multiple phased subspace vector by matrix and project to subspace
	  // Remove local bulk phase to leave relative phases
	  /////////////////////////////////////////////////////////////////////
	  // Memory footprint was an issue
	  tmat-=usecond();
	  linop.Op(phaV,MphaV);
	  _MphaV[b] = MphaV;
	  tmat+=usecond();
	}      

	//	std::cout << GridLogMessage << " Calling block project "<<std::endl;
	tproj-=usecond();
	Projector.blockProject(_MphaV,TmpProj);
	tproj+=usecond();
	
	//	std::cout << GridLogMessage << " conj phasing the coarse vectors "<<std::endl;
	for(int b=0;b<MIN(batch,npoint-p);b++){
	  ComputeProj[p+b] = conjugate(pha[p+b])*TmpProj[b];
	}
      }

      // Could do this with a block promote or similar BLAS call via the MultiRHSBlockProjector with a const matrix.
      
      // std::cout << GridLogMessage << " Starting FT inv "<<std::endl;
      tinv-=usecond();
      for(int k=0;k<npoint;k++){
	FT = Zero();
	// 81 kernel calls as many ComputeProj vectors
	// Could fuse with a vector of views, but ugly
	// Could unroll the expression and run fewer kernels -- much more attractive
	// Could also do non blocking.
#if 0	
	for(int l=0;l<npoint;l++){
	  FT= FT+ invMkl(l,k)*ComputeProj[l];
	}
#else
	const int radix = 9;
	int ll;
	for(ll=0;ll+radix-1<npoint;ll+=radix){
	  // When ll = npoint-radix, ll+radix-1 = npoint-1, and we do it all.
	  FT = FT 
	    + invMkl(ll+0,k)*ComputeProj[ll+0]
	    + invMkl(ll+1,k)*ComputeProj[ll+1]
	    + invMkl(ll+2,k)*ComputeProj[ll+2]
	    + invMkl(ll+3,k)*ComputeProj[ll+3]
	    + invMkl(ll+4,k)*ComputeProj[ll+4]
	    + invMkl(ll+5,k)*ComputeProj[ll+5]
	    + invMkl(ll+6,k)*ComputeProj[ll+6]
	    + invMkl(ll+7,k)*ComputeProj[ll+7]
	    + invMkl(ll+8,k)*ComputeProj[ll+8];
	}
	for(int l=ll;l<npoint;l++){
	  FT= FT+ invMkl(l,k)*ComputeProj[l];
	}
#endif
      
	// 1 kernel call -- must be cheaper
	int osites=CoarseGrid->oSites();
	autoView( A_v  , _A[k], AcceleratorWrite);
	autoView( FT_v  , FT, AcceleratorRead);
	accelerator_for(sss, osites, 1, {
	    for(int j=0;j<nbasis;j++){
	      A_v[sss](i,j) = FT_v[sss](j);
	    }
        });
      }
      tinv+=usecond();
    }

    // Only needed if nonhermitian
    //    if ( ! hermitian ) {
    //      std::cout << GridLogMessage<<"PopulateAdag  "<<std::endl;
    //      PopulateAdag();
    //    }
    // Need to write something to populate Adag from A
    //    std::cout << GridLogMessage << " Calling GridtoBLAS "<<std::endl;
    for(int p=0;p<geom_srhs.npoint;p++){
      GridtoBLAS(_A[p],BLAS_A[p]);
    }
    std::cout << GridLogMessage<<"CoarsenOperator phase  "<<tphase<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator phaseBZ "<<tphaseBZ<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator mat    "<<tmat <<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator proj   "<<tproj<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator inv    "<<tinv<<" us"<<std::endl;
#endif
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

    t_tot=-usecond();
    CoarseVector tin=in;
    t_exch=-usecond();
    CoarseVector pin = Cell.ExchangePeriodic(tin); //padded input
    t_exch+=usecond();

    CoarseVector pout(pin.Grid());

    int npoint = geom.npoint;
    typedef calcMatrix* Aview;
    typedef LatticeView<Cvec> Vview;
      
    const int Nsimd = CComplex::Nsimd();

    int64_t nrhs  =pin.Grid()->GlobalDimensions()[0];
    assert(nrhs>=1);

    RealD flops,bytes;
    int64_t osites=in.Grid()->oSites(); // unpadded
    int64_t unpadded_vol = CoarseGrid()->lSites()/nrhs;
    
    flops = 1.0* npoint * nbasis * nbasis * 8.0 * osites * CComplex::Nsimd();
    bytes = 1.0*osites*sizeof(siteMatrix)*npoint/pin.Grid()->GlobalDimensions()[0]
          + 2.0*osites*sizeof(siteVector)*npoint;
    

    t_GtoB=-usecond();
    GridtoBLAS(pin,BLAS_B);
    t_GtoB+=usecond();

    GridBLAS BLAS;

    t_mult=-usecond();
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
    t_mult+=usecond();

    t_BtoG=-usecond();
    BLAStoGrid(out,BLAS_C);
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
  virtual  void Mdiag    (const Field &in, Field &out){ assert(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};
};
  
NAMESPACE_END(Grid);
