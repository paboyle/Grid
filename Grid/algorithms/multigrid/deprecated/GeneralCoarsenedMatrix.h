/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/GeneralCoarsenedMatrix.h

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

#include <Grid/qcd/QCD.h> // needed for Dagger(Yes|No), Inverse(Yes|No)

#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/deflation/MultiRHSBlockProject.h>

NAMESPACE_BEGIN(Grid);

// Fine Object == (per site) type of fine field
// nbasis      == number of deflation vectors
template<class Fobj,class CComplex,int nbasis>
class GeneralCoarsenedMatrix : public SparseMatrixBase<Lattice<iVector<CComplex,nbasis > > >  {
public:

  typedef GeneralCoarsenedMatrix<Fobj,CComplex,nbasis> GeneralCoarseOp;
  typedef iVector<CComplex,nbasis >           siteVector;
  typedef iMatrix<CComplex,nbasis >           siteMatrix;
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
  int hermitian;
  GridBase      *       _FineGrid; 
  GridCartesian *       _CoarseGrid; 
  NonLocalStencilGeometry &geom;
  PaddedCell Cell;
  GeneralLocalStencil Stencil;
  
  std::vector<CoarseMatrix> _A;
  std::vector<CoarseMatrix> _Adag;
  std::vector<CoarseVector> MultTemporaries;

  int64_t MultCalls;
  double  MultFlopsAccum;
  double  MultUsecAccum;

  ///////////////////////
  // Interface
  ///////////////////////
  //////////////////////////////////////////////////////////////////////////
  // Bilingual accessors: everything a consumer needs to read the operator
  // without knowing which of the three coarse classes it holds. The D
  // dimensional grid the elements live on, the geometry they are indexed by,
  // and one unpadded point at a time (a whole npoint vector is too much
  // memory at production nbasis).
  //////////////////////////////////////////////////////////////////////////
  GridCartesian * CoarseGridD(void)             { return _CoarseGrid; };
  NonLocalStencilGeometry & Geometry(void)      { return geom; };
  void ExtractMatrix(int p,CoarseMatrix &A)     { A = Cell.Extract(_A[p]); };

  GridBase      * Grid(void)           { return _CoarseGrid; };   // this is all the linalg routines need to know
  GridBase      * FineGrid(void)       { return _FineGrid; };   // this is all the linalg routines need to know
  GridCartesian * CoarseGrid(void)     { return _CoarseGrid; };   // this is all the linalg routines need to know

  /*  void ShiftMatrix(RealD shift)
  {
    int Nd=_FineGrid->Nd(); 
    Coordinate zero_shift(Nd,0);
    for(int p=0;p<geom.npoint;p++){
      if ( zero_shift==geom.shifts[p] ) {
	_A[p] = _A[p]+shift;
	//	_Adag[p] = _Adag[p]+shift;
      }
    }    
  }
  void ProjectNearestNeighbour(RealD shift, GeneralCoarseOp &CopyMe)
  {
    int nfound=0;
    std::cout << GridLogMessage <<"GeneralCoarsenedMatrix::ProjectNearestNeighbour "<< CopyMe._A[0].Grid()<<std::endl;
    for(int p=0;p<geom.npoint;p++){
      for(int pp=0;pp<CopyMe.geom.npoint;pp++){
 	// Search for the same relative shift
	// Avoids brutal handling of Grid pointers
	if ( CopyMe.geom.shifts[pp]==geom.shifts[p] ) {
	  _A[p] = CopyMe.Cell.Extract(CopyMe._A[pp]);
	  //	  _Adag[p] = CopyMe.Cell.Extract(CopyMe._Adag[pp]);
	  nfound++;
	}
      }
    }
    GRID_ASSERT(nfound==geom.npoint);
    ExchangeCoarseLinks();
  }
  */
  
  GeneralCoarsenedMatrix(NonLocalStencilGeometry &_geom,GridBase *FineGrid, GridCartesian * CoarseGrid,int _herm=1)
    : geom(_geom),
      _FineGrid(FineGrid),
      _CoarseGrid(CoarseGrid),
      hermitian(_herm),
      Cell(_geom.Depth(),_CoarseGrid),
      Stencil(Cell.grids.back(),geom.shifts),
      MultCalls(0), MultFlopsAccum(0.0), MultUsecAccum(0.0)
  {
    {
      int npoint = _geom.npoint;
    }
    _A.resize(geom.npoint,CoarseGrid);
    if ( !hermitian ) _Adag.resize(geom.npoint,CoarseGrid);
  }
  void M (const CoarseVector &in, CoarseVector &out)
  {
    Mult(_A,in,out);
  }
  void Mdag (const CoarseVector &in, CoarseVector &out)
  {
    if(hermitian)
      Mult(_A,in,out);
    else
      Mult(_Adag,in,out);
  }
  void Mult (std::vector<CoarseMatrix> &A,const CoarseVector &in, CoarseVector &out)
  {
    RealD tviews=0;    RealD ttot=0;    RealD tmult=0;   RealD texch=0;    RealD text=0; RealD ttemps=0; RealD tcopy=0;
    RealD tmult2=0;

    ttot=-usecond();
    conformable(CoarseGrid(),in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();
    CoarseVector tin=in;

    texch-=usecond();
    CoarseVector pin = Cell.ExchangePeriodic(tin);
    texch+=usecond();

    CoarseVector pout(pin.Grid());

    int npoint = geom.npoint;
    typedef LatticeView<Cobj> Aview;
    typedef LatticeView<Cvec> Vview;
      
    const int Nsimd = CComplex::Nsimd();
    
    int64_t osites=pin.Grid()->oSites();

    RealD flops = 1.0* npoint * nbasis * nbasis * 8.0 * osites * CComplex::Nsimd();
    RealD bytes = 1.0*osites*sizeof(siteMatrix)*npoint
                + 2.0*osites*sizeof(siteVector)*npoint;
      
    {
      tviews-=usecond();
      autoView( in_v , pin, AcceleratorRead);
      autoView( out_v , pout, AcceleratorWriteDiscard);
      autoView( Stencil_v  , Stencil, AcceleratorRead);
      tviews+=usecond();

      // Static and prereserve to keep UVM region live and not resized across multiple calls
      ttemps-=usecond();
      MultTemporaries.resize(npoint,pin.Grid());       
      ttemps+=usecond();
      std::vector<Aview> AcceleratorViewContainer_h;
      std::vector<Vview> AcceleratorVecViewContainer_h; 

      tviews-=usecond();
      for(int p=0;p<npoint;p++) {
	AcceleratorViewContainer_h.push_back(      A[p].View(AcceleratorRead));
	AcceleratorVecViewContainer_h.push_back(MultTemporaries[p].View(AcceleratorWrite));
      }
      tviews+=usecond();

      static deviceVector<Aview> AcceleratorViewContainer; AcceleratorViewContainer.resize(npoint);
      static deviceVector<Vview> AcceleratorVecViewContainer; AcceleratorVecViewContainer.resize(npoint); 
      
      auto Aview_p = &AcceleratorViewContainer[0];
      auto Vview_p = &AcceleratorVecViewContainer[0];
      tcopy-=usecond();
      acceleratorCopyToDevice(&AcceleratorViewContainer_h[0],&AcceleratorViewContainer[0],npoint *sizeof(Aview));
      acceleratorCopyToDevice(&AcceleratorVecViewContainer_h[0],&AcceleratorVecViewContainer[0],npoint *sizeof(Vview));
      tcopy+=usecond();

      tmult-=usecond();
      accelerator_for(spb, osites*nbasis*npoint, Nsimd, {
	  typedef decltype(coalescedRead(in_v[0](0))) calcComplex;
	  int32_t ss   = spb/(nbasis*npoint);
	  int32_t bp   = spb%(nbasis*npoint);
	  int32_t point= bp/nbasis;
	  int32_t b    = bp%nbasis;
	  auto SE  = Stencil_v.GetEntry(point,ss);
	  auto nbr = coalescedReadGeneralPermute(in_v[SE->_offset],SE->_permute,Nd);
	  auto res = coalescedRead(Aview_p[point][ss](0,b))*nbr(0);
	  for(int bb=1;bb<nbasis;bb++) {
	    res = res + coalescedRead(Aview_p[point][ss](bb,b))*nbr(bb);
	  }
	  coalescedWrite(Vview_p[point][ss](b),res);
      });
      tmult2-=usecond();
      accelerator_for(sb, osites*nbasis, Nsimd, {
	  int ss = sb/nbasis;
	  int b  = sb%nbasis;
	  auto res = coalescedRead(Vview_p[0][ss](b));
	  for(int point=1;point<npoint;point++){
	    res = res + coalescedRead(Vview_p[point][ss](b));
	  }
	  coalescedWrite(out_v[ss](b),res);
      });
      tmult2+=usecond();
      tmult+=usecond();
      for(int p=0;p<npoint;p++) {
	AcceleratorViewContainer_h[p].ViewClose();
	AcceleratorVecViewContainer_h[p].ViewClose();
      }
    }

    text-=usecond();
    out = Cell.Extract(pout);
    text+=usecond();
    ttot+=usecond();
    
    MultCalls++;
    MultFlopsAccum += flops;
    MultUsecAccum  += ttot;
    std::cout << GridLogPerformance
	      << "Coarse Mult call " << MultCalls
	      << "  tot " << ttot << " us"
	      << "  kernel " << tmult << " us"
	      << "  kernel " << flops/tmult*1e-3 << " GFlop/s"
	      << "  overall " << MultFlopsAccum/MultUsecAccum*1e-3 << " GFlop/s (cumul)"
	      << "  bw " << bytes/tmult*1e-3 << " GB/s"
	      << std::endl;

  };
  
  void PopulateAdag(void)
  {
#if 0
    // Serial global peek/poke reference implementation
    for(int64_t bidx=0;bidx<CoarseGrid()->gSites() ;bidx++){
      Coordinate bcoor;
      CoarseGrid()->GlobalIndexToGlobalCoor(bidx,bcoor);

      for(int p=0;p<geom.npoint;p++){
	Coordinate scoor = bcoor;
	for(int mu=0;mu<bcoor.size();mu++){
	  int L = CoarseGrid()->GlobalDimensions()[mu];
	  scoor[mu] = (bcoor[mu] - geom.shifts[p][mu] + L) % L; // Modulo arithmetic
	}
	// Flip to poke/peekLocalSite and not too bad
	auto link = peekSite(_A[p],scoor);
	int pp = geom.Reverse(p);
	pokeSite(adj(link),_Adag[pp],bcoor);
      }
    }
#else
    // Parallel: _Adag[pp](x) = adj( _A[p](x + s_pp) ),  pp = Reverse(p), s_pp = -s_p.
    // The neighbour fetch reuses the same padded-cell + stencil machinery as Mult,
    // reading one matrix element per coalesced access so no whole site matrix
    // (230KB at nbasis=60) ever lands on a GPU thread stack (HIP limit 128KB).
    // Halo sites compute garbage neighbours; Cell.Extract discards them.
    // Must run on the unpadded _A, i.e. before ExchangeCoarseLinks.
    const int Nsimd = CComplex::Nsimd();
    for(int p=0;p<geom.npoint;p++){
      int pp = geom.Reverse(p);
      CoarseMatrix Apad = Cell.ExchangePeriodic(_A[p]);
      CoarseMatrix Dpad(Apad.Grid());
      int64_t osites = Apad.Grid()->oSites();
      {
	autoView( Apad_v   , Apad, AcceleratorRead);
	autoView( Dpad_v   , Dpad, AcceleratorWriteDiscard);
	autoView( Stencil_v, Stencil, AcceleratorRead);
	accelerator_for(sj, osites*nbasis, Nsimd, {
	    int32_t ss = sj/nbasis;
	    int32_t j  = sj%nbasis;
	    auto SE = Stencil_v.GetEntry(pp,ss);
	    for(int i=0;i<nbasis;i++){
	      auto z = coalescedReadGeneralPermute(Apad_v[SE->_offset](i,j),SE->_permute,Nd);
	      coalescedWrite(Dpad_v[ss](j,i),conjugate(z));
	    }
	});
      }
      _Adag[pp] = Cell.Extract(Dpad);
    }
#endif
  }
  /////////////////////////////////////////////////////////////
  // 
  // A) Only reduced flops option is to use a padded cell of depth 4
  // and apply MpcDagMpc in the padded cell.
  //
  // Makes for ONE application of MpcDagMpc per vector instead of 30 or 80.
  // With the effective cell size around (B+8)^4 perhaps 12^4/4^4 ratio
  // Cost is 81x more, same as stencil size.
  //
  // But: can eliminate comms and do as local dirichlet.
  //
  // Local exchange gauge field once.
  // Apply to all vectors, local only computation.
  // Must exchange ghost subcells in reverse process of PaddedCell to take inner products
  //
  // B) Can reduce cost: pad by 1, apply Deo      (4^4+6^4+8^4+8^4 )/ (4x 4^4)
  //                     pad by 2, apply Doe
  //                     pad by 3, apply Deo
  //                     then break out 8x directions; cost is ~10x MpcDagMpc per vector
  //
  // => almost factor of 10 in setup cost, excluding data rearrangement
  //
  // Intermediates -- ignore the corner terms, leave approximate and force Hermitian
  // Intermediates -- pad by 2 and apply 1+8+24 = 33 times.
  /////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    // BFM HDCG style approach: Solve a system of equations to get Aij
    //////////////////////////////////////////////////////////
    /*
     *     Here, k,l index which possible shift within the 3^Nd "ball" connected by MdagM.
     *
     *     conj(phases[block]) proj[k][ block*Nvec+j ] =  \sum_ball  e^{i q_k . delta} < phi_{block,j} | MdagM | phi_{(block+delta),i} > 
     *                                                 =  \sum_ball e^{iqk.delta} A_ji
     *
     *     Must invert matrix M_k,l = e^[i q_k . delta_l]
     *
     *     Where q_k = delta_k . (2*M_PI/global_nb[mu])
     */
#if 0
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & Subspace)
  {
    std::cout << GridLogMessage<< "GeneralCoarsenMatrix "<< std::endl;
    GridBase *grid = FineGrid();

    RealD tproj=0.0;
    RealD teigen=0.0;
    RealD tmat=0.0;
    RealD tphase=0.0;
    RealD tinv=0.0;

    /////////////////////////////////////////////////////////////
    // Orthogonalise the subblocks over the basis
    /////////////////////////////////////////////////////////////
    CoarseScalar InnerProd(CoarseGrid()); 
    blockOrthogonalise(InnerProd,Subspace.subspace);

    const int npoint = geom.npoint;
      
    Coordinate clatt = CoarseGrid()->GlobalDimensions();
    int Nd = CoarseGrid()->Nd();

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
    teigen-=usecond();
    Eigen::MatrixXcd Mkl    = Eigen::MatrixXcd::Zero(npoint,npoint);
    Eigen::MatrixXcd invMkl = Eigen::MatrixXcd::Zero(npoint,npoint);
    ComplexD ci(0.0,1.0);
    for(int k=0;k<npoint;k++){ // Loop over momenta

      for(int l=0;l<npoint;l++){ // Loop over nbr relative
	ComplexD phase(0.0,0.0);
	for(int mu=0;mu<Nd;mu++){
	  RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	  phase=phase+TwoPiL*geom.shifts[k][mu]*geom.shifts[l][mu];
	}
	phase=exp(phase*ci);
	Mkl(k,l) = phase;
      }
    }
    invMkl = Mkl.inverse();
    teigen+=usecond();

    ///////////////////////////////////////////////////////////////////////
    // Now compute the matrix elements of linop between the orthonormal
    // set of vectors.
    ///////////////////////////////////////////////////////////////////////
    FineField phaV(grid); // Phased block basis vector
    FineField MphaV(grid);// Matrix applied
    CoarseVector coarseInner(CoarseGrid());

    std::vector<CoarseVector> ComputeProj(npoint,CoarseGrid());
    std::vector<CoarseVector>          FT(npoint,CoarseGrid());
    for(int i=0;i<nbasis;i++){// Loop over basis vectors
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i<<"/"<<nbasis<< std::endl;
      for(int p=0;p<npoint;p++){ // Loop over momenta in npoint
	/////////////////////////////////////////////////////
	// Stick a phase on every block
	/////////////////////////////////////////////////////
	tphase-=usecond();
	CoarseComplexField coor(CoarseGrid());
	CoarseComplexField pha(CoarseGrid());	pha=Zero();
	for(int mu=0;mu<Nd;mu++){
	  LatticeCoordinate(coor,mu);
	  RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	  pha = pha + (TwoPiL * geom.shifts[p][mu]) * coor;
	}
	pha  =exp(pha*ci);
	phaV=Zero();
	blockZAXPY(phaV,pha,Subspace.subspace[i],phaV);
	tphase+=usecond();

	/////////////////////////////////////////////////////////////////////
	// Multiple phased subspace vector by matrix and project to subspace
	// Remove local bulk phase to leave relative phases
	/////////////////////////////////////////////////////////////////////
	tmat-=usecond();
	linop.Op(phaV,MphaV);
	tmat+=usecond();

	tproj-=usecond();
	blockProject(coarseInner,MphaV,Subspace.subspace);
	coarseInner = conjugate(pha) * coarseInner;

	ComputeProj[p] = coarseInner;
	tproj+=usecond();

      }

      tinv-=usecond();
      for(int k=0;k<npoint;k++){
	FT[k] = Zero();
	for(int l=0;l<npoint;l++){
	  FT[k]= FT[k]+ invMkl(l,k)*ComputeProj[l];
	}
      
	int osites=CoarseGrid()->oSites();
	autoView( A_v  , _A[k], AcceleratorWrite);
	autoView( FT_v  , FT[k], AcceleratorRead);
	accelerator_for(sss, osites, nbasis, {
#ifdef GRID_SIMT
	    int j = acceleratorSIMTlane(nbasis);
	    A_v[sss](i,j) = FT_v[sss](j);
#else
	    // CPU build: acceleratorSIMTlane()==0 -- an un-looped SIMT tensor
	    // index writes ONLY j=0 and silently drops the other nbasis-1
	    // columns (caught by Test_schur_dense_coarse import certificate,
	    // 2026-08-14).  Loop explicitly.
	    for(int j=0;j<nbasis;j++){
	      A_v[sss](i,j) = FT_v[sss](j);
	    }
#endif
        });
      }
      tinv+=usecond();
    }

    // Only needed if nonhermitian
    if ( ! hermitian ) {
      std::cout << GridLogMessage<<"PopulateAdag  "<<std::endl;
      PopulateAdag();
    }

    // Need to write something to populate Adag from A
    ExchangeCoarseLinks();
    std::cout << GridLogMessage<<"CoarsenOperator eigen  "<<teigen<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator phase  "<<tphase<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator mat    "<<tmat <<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator proj   "<<tproj<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator inv    "<<tinv<<" us"<<std::endl;
  }
#else
  //////////////////////////////////////////////////////////////////////
  // Galerkin projection of matrix
  //////////////////////////////////////////////////////////////////////
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & Subspace)
  {
    CoarsenOperator(linop,Subspace,Subspace);
  }
  //////////////////////////////////////////////////////////////////////
  // Petrov - Galerkin projection of matrix
  //////////////////////////////////////////////////////////////////////
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & U,
		       Aggregation<Fobj,CComplex,nbasis> & V)
  {
    std::cout << GridLogMessage<< "GeneralCoarsenMatrix "<< std::endl;
    GridBase *grid = FineGrid();

    RealD tproj=0.0;
    RealD teigen=0.0;
    RealD tmat=0.0;
    RealD tphase=0.0;
    RealD tphaseBZ=0.0;
    RealD tinv=0.0;

    /////////////////////////////////////////////////////////////
    // Orthogonalise the subblocks over the basis
    /////////////////////////////////////////////////////////////
    CoarseScalar InnerProd(CoarseGrid()); 
    blockOrthogonalise(InnerProd,V.subspace);
    blockOrthogonalise(InnerProd,U.subspace);

    const int npoint = geom.npoint;
      
    Coordinate clatt = CoarseGrid()->GlobalDimensions();
    int Nd = CoarseGrid()->Nd();

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
    teigen-=usecond();
    Eigen::MatrixXcd Mkl    = Eigen::MatrixXcd::Zero(npoint,npoint);
    Eigen::MatrixXcd invMkl = Eigen::MatrixXcd::Zero(npoint,npoint);
    ComplexD ci(0.0,1.0);
    for(int k=0;k<npoint;k++){ // Loop over momenta

      for(int l=0;l<npoint;l++){ // Loop over nbr relative
	ComplexD phase(0.0,0.0);
	for(int mu=0;mu<Nd;mu++){
	  RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	  phase=phase+TwoPiL*geom.shifts[k][mu]*geom.shifts[l][mu];
	}
	phase=exp(phase*ci);
	Mkl(k,l) = phase;
      }
    }
    invMkl = Mkl.inverse();
    teigen+=usecond();

    ///////////////////////////////////////////////////////////////////////
    // Now compute the matrix elements of linop between the orthonormal
    // set of vectors.
    ///////////////////////////////////////////////////////////////////////
    std::vector<FineComplexField> phaF(npoint,grid);
    std::vector<CoarseComplexField> pha(npoint,CoarseGrid());

    typedef typename CComplex::scalar_type SComplex;
    FineComplexField one(grid); one=SComplex(1.0);
    FineComplexField zz(grid); zz = Zero();
    tphase=-usecond();
    for(int p=0;p<npoint;p++){ // Loop over momenta in npoint
      /////////////////////////////////////////////////////
      // Stick a phase on every block
      /////////////////////////////////////////////////////
      CoarseComplexField coor(CoarseGrid());
      pha[p]=Zero();
      for(int mu=0;mu<Nd;mu++){
	LatticeCoordinate(coor,mu);
	RealD TwoPiL =  M_PI * 2.0/ clatt[mu];
	pha[p] = pha[p] + (TwoPiL * geom.shifts[p][mu]) * coor;
      }
      pha[p]  =exp(pha[p]*ci);

      blockZAXPY(phaF[p],pha[p],one,zz);

    }
    tphase+=usecond();

    // Import basis into BLAS layout once; blockProject then reads it once per
    // basis vector rather than once per (i,p) as in scalar blockProject.
    // Process all npoint in a single batch.
    MultiRHSBlockProject<FineField> Projector;
    Projector.Allocate(nbasis, grid, CoarseGrid());
    Projector.ImportBasis(U.subspace);

    std::vector<FineField>    phaV_batch(npoint, grid);
    std::vector<FineField>    MphaV_batch(npoint, grid);
    std::vector<CoarseVector> proj_batch(npoint, CoarseGrid());
    std::vector<CoarseVector> ComputeProj(npoint, CoarseGrid());
    std::vector<CoarseVector>          FT(npoint, CoarseGrid());

    // Pre-allocate BLAS_F and BLAS_C to avoid repeated hipMalloc/hipFree of
    // ~5.6 GB per blockProject call, which hangs on ROCm for large allocations.
    Projector.BLAS_F.resize(Projector.fine_vol * Projector.words * npoint);
    Projector.BLAS_C.resize(Projector.coarse_vol * nbasis * npoint);

    for(int i=0;i<nbasis;i++){// Loop over basis vectors
      accelerator_barrier(); // ensure prior iteration's async writes are retired
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i<<"/"<<nbasis<< std::endl;

      tphaseBZ-=usecond();
      for(int p=0;p<npoint;p++)
	phaV_batch[p] = phaF[p] * V.subspace[i];
      tphaseBZ+=usecond();
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i<<" phaseBZ done"<< std::endl;

      tmat-=usecond();
      for(int p=0;p<npoint;p++)
	linop.Op(phaV_batch[p], MphaV_batch[p]);
      tmat+=usecond();
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i<<" mat done"<< std::endl;

      // One batched GEMM reads BLAS_V once for all npoint vectors.
      tproj-=usecond();
      Projector.blockProject(MphaV_batch, proj_batch);
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i<<" blockProject done"<< std::endl;
      for(int p=0;p<npoint;p++)
	ComputeProj[p] = conjugate(pha[p]) * proj_batch[p];
      tproj+=usecond();
      std::cout << GridLogMessage<< "CoarsenMatrixColoured vec "<<i<<" proj done"<< std::endl;

      tinv-=usecond();
      for(int k=0;k<npoint;k++){
	FT[k] = Zero();
	for(int l=0;l<npoint;l++){
	  FT[k]= FT[k]+ invMkl(l,k)*ComputeProj[l];
	}

	int osites=CoarseGrid()->oSites();
	autoView( A_v  , _A[k], AcceleratorWrite);
	autoView( FT_v  , FT[k], AcceleratorRead);
	accelerator_for(sss, osites, nbasis, {
#ifdef GRID_SIMT
	    int j = acceleratorSIMTlane(nbasis);
	    A_v[sss](i,j) = FT_v[sss](j);
#else
	    // CPU build: acceleratorSIMTlane()==0 -- an un-looped SIMT tensor
	    // index writes ONLY j=0 and silently drops the other nbasis-1
	    // columns (caught by Test_schur_dense_coarse import certificate,
	    // 2026-08-14).  Loop explicitly.
	    for(int j=0;j<nbasis;j++){
	      A_v[sss](i,j) = FT_v[sss](j);
	    }
#endif
        });
      }
      tinv+=usecond();
    }

    // Only needed if nonhermitian
    if ( ! hermitian ) {
      std::cout << GridLogMessage<<"PopulateAdag  "<<std::endl;
      PopulateAdag();
    }

    //    for(int p=0;p<geom.npoint;p++){
    //      std::cout << " _A["<<p<<"] "<<norm2(_A[p])<<std::endl;
    //    }

    // Need to write something to populate Adag from A
    ExchangeCoarseLinks();
    std::cout << GridLogMessage<<"CoarsenOperator eigen  "<<teigen<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator phase  "<<tphase<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator phaseBZ "<<tphaseBZ<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator mat    "<<tmat <<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator proj   "<<tproj<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator inv    "<<tinv<<" us"<<std::endl;
  }
#endif  
  void ExchangeCoarseLinks(void){
    for(int p=0;p<geom.npoint;p++){
      _A[p] = Cell.ExchangePeriodic(_A[p]);
      if ( !hermitian ) _Adag[p]= Cell.ExchangePeriodic(_Adag[p]);
    }
  }
  virtual  void Mdiag    (const Field &in, Field &out){ GRID_ASSERT(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};
};


  
NAMESPACE_END(Grid);
