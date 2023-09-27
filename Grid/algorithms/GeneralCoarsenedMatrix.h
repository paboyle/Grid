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

NAMESPACE_BEGIN(Grid);

// Fixme need coalesced read gpermute
template<class vobj> void gpermute(vobj & inout,int perm){
  vobj tmp=inout;
  if (perm & 0x1 ) { permute(inout,tmp,0); tmp=inout;}
  if (perm & 0x2 ) { permute(inout,tmp,1); tmp=inout;}
  if (perm & 0x4 ) { permute(inout,tmp,2); tmp=inout;}
  if (perm & 0x8 ) { permute(inout,tmp,3); tmp=inout;}
}

/////////////////////////////////////////////////////////////////
// Reuse Aggregation class from CoarsenedMatrix for now
// Might think about *smoothed* Aggregation
// Equivalent of Geometry class in cartesian case
/////////////////////////////////////////////////////////////////
class NonLocalStencilGeometry {
public:
  int depth;
  int hops;
  int npoint;
  std::vector<Coordinate> shifts;
  Coordinate stencil_size;
  Coordinate stencil_lo;
  Coordinate stencil_hi;
  GridCartesian *grid;
  GridCartesian *Grid() {return grid;};
  int Depth(void){return 1;};   // Ghost zone depth
  int Hops(void){return hops;}; // # of hops=> level of corner fill in in stencil

  virtual int DimSkip(void) =0;

  virtual ~NonLocalStencilGeometry() {};

  int  Reverse(int point)
  {
    int Nd = Grid()->Nd();
    Coordinate shft = shifts[point];
    Coordinate rev(Nd);
    for(int mu=0;mu<Nd;mu++) rev[mu]= -shft[mu];
    for(int p=0;p<npoint;p++){
      if(rev==shifts[p]){
	return p;
      }
    }
    assert(0);
    return -1;
  }
  void BuildShifts(void)
  {
    this->shifts.resize(0);
    int Nd = this->grid->Nd();

    int dd = this->DimSkip();
    for(int s0=this->stencil_lo[dd+0];s0<=this->stencil_hi[dd+0];s0++){
    for(int s1=this->stencil_lo[dd+1];s1<=this->stencil_hi[dd+1];s1++){
    for(int s2=this->stencil_lo[dd+2];s2<=this->stencil_hi[dd+2];s2++){
    for(int s3=this->stencil_lo[dd+3];s3<=this->stencil_hi[dd+3];s3++){
      Coordinate sft(Nd,0);
      sft[dd+0] = s0;
      sft[dd+1] = s1;
      sft[dd+2] = s2;
      sft[dd+3] = s3;
      int nhops = abs(s0)+abs(s1)+abs(s2)+abs(s3);
      if(nhops<=this->hops) this->shifts.push_back(sft);
    }}}}
    this->npoint = this->shifts.size();
    std::cout << GridLogMessage << "NonLocalStencilGeometry has "<< this->npoint << " terms in stencil "<<std::endl;
  }
  
  NonLocalStencilGeometry(GridCartesian *_coarse_grid,int _hops) : grid(_coarse_grid), hops(_hops)
  {
    Coordinate latt = grid->GlobalDimensions();
    stencil_size.resize(grid->Nd());
    stencil_lo.resize(grid->Nd());
    stencil_hi.resize(grid->Nd());
    for(int d=0;d<grid->Nd();d++){
     if ( latt[d] == 1 ) {
      stencil_lo[d] = 0;
      stencil_hi[d] = 0;
      stencil_size[d]= 1;
     } else if ( latt[d] == 2 ) {
      stencil_lo[d] = -1;
      stencil_hi[d] = 0;
      stencil_size[d]= 2;
     } else if ( latt[d] > 2 ) {
       stencil_lo[d] = -1;
       stencil_hi[d] =  1;
       stencil_size[d]= 3;
     }
    }
  };

};

// Need to worry about red-black now
class NonLocalStencilGeometry4D : public NonLocalStencilGeometry {
public:
  virtual int DimSkip(void) { return 0;};
  NonLocalStencilGeometry4D(GridCartesian *Coarse,int _hops) : NonLocalStencilGeometry(Coarse,_hops) { };
  virtual ~NonLocalStencilGeometry4D() {};
};
class NonLocalStencilGeometry5D : public NonLocalStencilGeometry {
public:
  virtual int DimSkip(void) { return 1; }; 
  NonLocalStencilGeometry5D(GridCartesian *Coarse,int _hops) : NonLocalStencilGeometry(Coarse,_hops)  { };
  virtual ~NonLocalStencilGeometry5D() {};
};
/*
 * Bunch of different options classes
 */
class NextToNextToNextToNearestStencilGeometry4D : public NonLocalStencilGeometry4D {
public:
  NextToNextToNextToNearestStencilGeometry4D(GridCartesian *Coarse) :  NonLocalStencilGeometry4D(Coarse,4)
  {
    this->BuildShifts();
  };
};
class NextToNextToNextToNearestStencilGeometry5D : public  NonLocalStencilGeometry5D {
public:
  NextToNextToNextToNearestStencilGeometry5D(GridCartesian *Coarse) :  NonLocalStencilGeometry5D(Coarse,4)
  {
    this->BuildShifts();
  };
};
class NextToNearestStencilGeometry4D : public  NonLocalStencilGeometry4D {
public:
  NextToNearestStencilGeometry4D(GridCartesian *Coarse) :  NonLocalStencilGeometry4D(Coarse,2)
  {
    this->BuildShifts();
  };
};
class NextToNearestStencilGeometry5D : public  NonLocalStencilGeometry5D {
public:
  NextToNearestStencilGeometry5D(GridCartesian *Coarse) :  NonLocalStencilGeometry5D(Coarse,2)
  {
    this->BuildShifts();
  };
};
class NearestStencilGeometry4D : public  NonLocalStencilGeometry4D {
public:
  NearestStencilGeometry4D(GridCartesian *Coarse) :  NonLocalStencilGeometry4D(Coarse,1)
  {
    this->BuildShifts();
  };
};
class NearestStencilGeometry5D : public  NonLocalStencilGeometry5D {
public:
  NearestStencilGeometry5D(GridCartesian *Coarse) :  NonLocalStencilGeometry5D(Coarse,1)
  {
    this->BuildShifts();
  };
};

// Fine Object == (per site) type of fine field
// nbasis      == number of deflation vectors
template<class Fobj,class CComplex,int nbasis>
class GeneralCoarsenedMatrix : public SparseMatrixBase<Lattice<iVector<CComplex,nbasis > > >  {
public:
    
  typedef iVector<CComplex,nbasis >           siteVector;
  typedef Lattice<iScalar<CComplex> >         CoarseComplexField;
  typedef Lattice<siteVector>                 CoarseVector;
  typedef Lattice<iMatrix<CComplex,nbasis > > CoarseMatrix;
  typedef iMatrix<CComplex,nbasis >  Cobj;
  typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj >        FineField;
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

  ///////////////////////
  // Interface
  ///////////////////////
  GridBase      * Grid(void)           { return _FineGrid; };   // this is all the linalg routines need to know
  GridBase      * FineGrid(void)       { return _FineGrid; };   // this is all the linalg routines need to know
  GridCartesian * CoarseGrid(void)     { return _CoarseGrid; };   // this is all the linalg routines need to know


  void ProjectNearestNeighbour(RealD shift)
  {
    int Nd = geom.grid->Nd();
    int point;
    std::cout << "ProjectNearestNeighbour "<<std::endl;
    for(int p=0;p<geom.npoint;p++){
      int nhops = 0;
      for(int s=0;s<Nd;s++){
	nhops+=abs(geom.shifts[p][s]);
      }
      if(nhops>1) {
	std::cout << "setting geom "<<p<<" shift "<<geom.shifts[p]<<" to zero "<<std::endl;
	_A[p]=Zero();
	_Adag[p]=Zero();
      }
      if(nhops==0) {
	std::cout << " Adding IR shift "<<shift<<" to "<<geom.shifts[p]<<std::endl;
	_A[p]=_A[p]+shift;
	_Adag[p]=_Adag[p]+shift;
      }
    }
  }
  
  GeneralCoarsenedMatrix(NonLocalStencilGeometry &_geom,GridBase *FineGrid, GridCartesian * CoarseGrid)
    : geom(_geom),
      _FineGrid(FineGrid),
      _CoarseGrid(CoarseGrid),
      hermitian(1),
      Cell(_geom.Depth(),_CoarseGrid),
      Stencil(Cell.grids.back(),geom.shifts)
  {
    {
      int npoint = _geom.npoint;
      StencilEntry *SE;
      autoView( Stencil_v  , Stencil, AcceleratorRead);
      int osites=Stencil.Grid()->oSites();
      for(int ss=0;ss<osites;ss++){
	for(int point=0;point<npoint;point++){
	  auto SE = Stencil_v.GetEntry(point,ss);
	  int o = SE->_offset;
	  assert( o< osites);
	}
      }    
    }

    _A.resize(geom.npoint,CoarseGrid);
    _Adag.resize(geom.npoint,CoarseGrid);
  }
  void M (const CoarseVector &in, CoarseVector &out)
  {
    Mult(_A,in,out);
  }
  void Mdag (const CoarseVector &in, CoarseVector &out)
  {
    Mult(_Adag,in,out);
  }
  void Mult (std::vector<CoarseMatrix> &A,const CoarseVector &in, CoarseVector &out)
  {
    RealD ttot=0;
    RealD tmult=0;
    RealD texch=0;
    RealD text=0;
    ttot=-usecond();
    conformable(CoarseGrid(),in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();
    CoarseVector tin=in;

    texch-=usecond();
    CoarseVector pin  = Cell.Exchange(tin);
    texch+=usecond();

    CoarseVector pout(pin.Grid());
    
    autoView( in_v , pin, AcceleratorRead);
    autoView( out_v , pout, AcceleratorWrite);
    autoView( Stencil_v  , Stencil, AcceleratorRead);
    int npoint = geom.npoint;
    typedef LatticeView<Cobj> Aview;
      
    Vector<Aview> AcceleratorViewContainer;
  
    for(int p=0;p<npoint;p++) AcceleratorViewContainer.push_back(A[p].View(AcceleratorRead));
    
    Aview *Aview_p = & AcceleratorViewContainer[0];

    const int Nsimd = CComplex::Nsimd();
    typedef siteVector calcVector;
    typedef CComplex   calcComplex;
    
    int osites=pin.Grid()->oSites();
    int gsites=pin.Grid()->gSites();

    RealD flops = 1.0* npoint * nbasis * nbasis * 8 * gsites;
      
    for(int point=0;point<npoint;point++){
      conformable(A[point],pin);
    }
    
    tmult-=usecond();
    accelerator_for(sss, osites*nbasis, 1, {
      int ss = sss/nbasis;
      int b  = sss%nbasis;
      assert(ss<osites);
      calcComplex res;
      res = Zero();
      calcVector nbr;
      int ptype;
      StencilEntry *SE;
      
      for(int point=0;point<npoint;point++){

	auto SE = Stencil_v.GetEntry(point,ss);
	
	int o = SE->_offset;

	assert( o < osites);
	// gpermute etc..
	nbr = in_v[o];
	gpermute(nbr,SE->_permute);

	for(int bb=0;bb<nbasis;bb++) {
	  res = res + Aview_p[point][ss](b,bb)*nbr(bb);
	}
      }
      out_v[ss](b)=res;
    });
    tmult+=usecond();

    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer[p].ViewClose();
    text-=usecond();
    out = Cell.Extract(pout);
    text+=usecond();
    ttot+=usecond();
    std::cout << GridLogMessage<<"Coarse Mult exch "<<texch<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult mult "<<tmult<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult ext  "<<text<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Mult tot  "<<ttot<<" us"<<std::endl;
    std::cout << GridLogMessage<<"Coarse Kernel flops/s "<< flops/tmult<<" mflop/s"<<std::endl;
    std::cout << GridLogMessage<<"Coarse flops/s "<< flops/ttot<<" mflop/s"<<std::endl;
  };

  void PopulateAdag(void)
  {
    for(int bidx=0;bidx<CoarseGrid()->gSites() ;bidx++){
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
  }
  void CoarsenOperator(LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & Subspace)
  {
    RealD tproj=0.0;
    RealD tpick=0.0;
    RealD tmat=0.0;
    RealD tpeek=0.0;
    std::cout << GridLogMessage<< "CoarsenMatrix "<< std::endl;
    GridBase *grid = FineGrid();

    ////////////////////////////////////////////////
    // Orthogonalise the subblocks over the basis
    ////////////////////////////////////////////////
    CoarseScalar InnerProd(CoarseGrid()); 
    blockOrthogonalise(InnerProd,Subspace.subspace);

    ////////////////////////////////////////////////
    // Now compute the matrix elements of linop between this orthonormal
    // set of vectors.
    ////////////////////////////////////////////////
    FineField bV(grid);
    FineField MbV(grid);
    FineField tmp(grid);
    CoarseVector coarseInner(CoarseGrid());
    
    // Very inefficient loop of order coarse volume.
    // First pass hack
    // Could replace with a coloring scheme our phase scheme
    // as in BFM
    for(int bidx=0;bidx<CoarseGrid()->gSites() ;bidx++){
      Coordinate bcoor;
      CoarseGrid()->GlobalIndexToGlobalCoor(bidx,bcoor);

      for(int b=0;b<nbasis;b++){
	tpick-=usecond();
	blockPick(CoarseGrid(),Subspace.subspace[b],bV,bcoor);
	tpick+=usecond();
	tmat-=usecond();
	linop.Op(bV,MbV);
	tmat+=usecond();
	tproj-=usecond();
	blockProject(coarseInner,MbV,Subspace.subspace);
	tproj+=usecond();

	tpeek-=usecond();
	for(int p=0;p<geom.npoint;p++){
	  Coordinate scoor = bcoor;
	  for(int mu=0;mu<bcoor.size();mu++){
	    int L = CoarseGrid()->GlobalDimensions()[mu];
	    scoor[mu] = (bcoor[mu] - geom.shifts[p][mu] + L) % L; // Modulo arithmetic
	  }
	  // Flip to peekLocalSite
	  // Flip to pokeLocalSite
	  auto ip    = peekSite(coarseInner,scoor);
	  auto Ab    = peekSite(_A[p],scoor);
	  int pp = geom.Reverse(p);
	  auto Adagb = peekSite(_Adag[pp],bcoor);
	  for(int bb=0;bb<nbasis;bb++){
	    Ab(bb,b)    = ip(bb);
	    Adagb(b,bb) = conjugate(ip(bb));
	  }
	  pokeSite(Ab,_A[p],scoor);
	  pokeSite(Adagb,_Adag[pp],bcoor);
	}
	tpeek+=usecond();
      }
    }
    for(int p=0;p<geom.npoint;p++){
      Coordinate coor({0,0,0,0,0});
      auto sval = peekSite(_A[p],coor);
    }
    for(int p=0;p<geom.npoint;p++){
      _A[p] = Cell.Exchange(_A[p]);
      _Adag[p]= Cell.Exchange(_Adag[p]);
    }
    std::cout << GridLogMessage<<"CoarsenOperator pick "<<tpick<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator mat  "<<tmat <<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator projection "<<tproj<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator peek/poke  "<<tpeek<<" us"<<std::endl;
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
  void CoarsenOperatorColoured(LinearOperatorBase<Lattice<Fobj> > &linop,
			       Aggregation<Fobj,CComplex,nbasis> & Subspace)
  {
    std::cout << GridLogMessage<< "CoarsenMatrixColoured "<< std::endl;
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
    ComplexD ci(0.0,1.0);
    Eigen::MatrixXcd Mkl    = Eigen::MatrixXcd::Zero(npoint,npoint);
    Eigen::MatrixXcd invMkl = Eigen::MatrixXcd::Zero(npoint,npoint);
    for(int k=0;k<npoint;k++){ // Loop over momenta

      for(int l=0;l<npoint;l++){ // Loop over nbr relative
	std::complex<double> phase(0.0,0.0);
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
	accelerator_for(sss, osites, 1, {
	    for(int j=0;j<nbasis;j++){
	      A_v[sss](j,i) = FT_v[sss](j);
	    }
        });
      }
      tinv+=usecond();
    }

    for(int p=0;p<geom.npoint;p++){
      Coordinate coor({0,0,0,0,0});
      auto sval = peekSite(_A[p],coor);
    }

    PopulateAdag();

    // Need to write something to populate Adag from A
    for(int p=0;p<geom.npoint;p++){
      _A[p] = Cell.Exchange(_A[p]);
      _Adag[p]= Cell.Exchange(_Adag[p]);
    }
    std::cout << GridLogMessage<<"CoarsenOperator eigen  "<<teigen<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator phase  "<<tphase<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator mat    "<<tmat <<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator proj   "<<tproj<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator inv    "<<tinv<<" us"<<std::endl;
  }

  virtual  void Mdiag    (const Field &in, Field &out){ assert(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};
};


NAMESPACE_END(Grid);
