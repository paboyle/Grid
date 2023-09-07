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
  int npoint;
  std::vector<Coordinate> shifts;
  virtual void BuildShifts(void) { assert(0); } ;
  int Depth(void){return depth;};
  NonLocalStencilGeometry(int _depth) : depth(_depth)
  {
  };
  virtual ~NonLocalStencilGeometry() {};
};
// Need to worry about red-black now
class NextToNearestStencilGeometry4D : public NonLocalStencilGeometry {
public:
  NextToNearestStencilGeometry4D(void) : NonLocalStencilGeometry(1)
  {
      this->BuildShifts();
  };
  virtual ~NextToNearestStencilGeometry4D() {};
  virtual void BuildShifts(void)
  {
    this->shifts.resize(0);
    // Like HDCG:  81 point stencil including self connection
    this->shifts.push_back(Coordinate({0,0,0,0}));
    // +-x, +-y, +-z, +-t  : 8
    for(int s=-1;s<=1;s+=2){
      this->shifts.push_back(Coordinate({s,0,0,0}));
      this->shifts.push_back(Coordinate({0,s,0,0}));
      this->shifts.push_back(Coordinate({0,0,s,0}));
      this->shifts.push_back(Coordinate({0,0,0,s}));
    }
    // +-x+-y, +-x+-z, +-x+-t, +-y+-z, +-y+-t, +-z+-t : 24
    for(int s1=-1;s1<=1;s1+=2){
    for(int s2=-1;s2<=1;s2+=2){
      this->shifts.push_back(Coordinate({s1,s2,0,0}));
      this->shifts.push_back(Coordinate({s1,0,s2,0}));
      this->shifts.push_back(Coordinate({s1,0,0,s2}));
      this->shifts.push_back(Coordinate({0,s1,s2,0}));
      this->shifts.push_back(Coordinate({0,s1,0,s2}));
      this->shifts.push_back(Coordinate({0,0,s1,s2}));
    }}
    this->npoint = this->shifts.size();
  }
};
// Need to worry about red-black now
class NextToNextToNextToNearestStencilGeometry4D : public NonLocalStencilGeometry {
public:
  NextToNextToNextToNearestStencilGeometry4D(void) : NonLocalStencilGeometry(1)
  {
    this->BuildShifts();
  };
  virtual ~NextToNextToNextToNearestStencilGeometry4D() {}
  virtual void BuildShifts(void)
  {
    this->shifts.resize(0);
    // Like HDCG:  81 point stencil including self connection
    this->shifts.push_back(Coordinate({0,0,0,0}));
    // +-x, +-y, +-z, +-t  : 8
    for(int s=-1;s<=1;s+=2){
      this->shifts.push_back(Coordinate({s,0,0,0}));
      this->shifts.push_back(Coordinate({0,s,0,0}));
      this->shifts.push_back(Coordinate({0,0,s,0}));
      this->shifts.push_back(Coordinate({0,0,0,s}));
    }
    // +-x+-y, +-x+-z, +-x+-t, +-y+-z, +-y+-t, +-z+-t : 24
    for(int s1=-1;s1<=1;s1+=2){
    for(int s2=-1;s2<=1;s2+=2){
      this->shifts.push_back(Coordinate({s1,s2,0,0}));
      this->shifts.push_back(Coordinate({s1,0,s2,0}));
      this->shifts.push_back(Coordinate({s1,0,0,s2}));
      this->shifts.push_back(Coordinate({0,s1,s2,0}));
      this->shifts.push_back(Coordinate({0,s1,0,s2}));
      this->shifts.push_back(Coordinate({0,0,s1,s2}));
    }}
    // +-x+-y+-z, +-x+-y+-z, +-x+-y+-z,
    for(int s1=-1;s1<=1;s1+=2){
    for(int s2=-1;s2<=1;s2+=2){
    for(int s3=-1;s3<=1;s3+=2){
      this->shifts.push_back(Coordinate({s1,s2,s3,0})); // 8x4 = 32
      this->shifts.push_back(Coordinate({s1,s2,0,s3}));
      this->shifts.push_back(Coordinate({s1,0,s2,s3}));
      this->shifts.push_back(Coordinate({0,s1,s2,s3}));
    }}}
    for(int s1=-1;s1<=1;s1+=2){
    for(int s2=-1;s2<=1;s2+=2){
    for(int s3=-1;s3<=1;s3+=2){
    for(int s4=-1;s4<=1;s4+=2){
      this->shifts.push_back(Coordinate({s1,s2,s3,s4})); // 16 
    }}}}
    this->npoint = this->shifts.size();
  }
};
class NextToNearestStencilGeometry5D : public NonLocalStencilGeometry {
public:
  NextToNearestStencilGeometry5D(void) : NonLocalStencilGeometry(1)
  {
      this->BuildShifts();
  };
  virtual ~NextToNearestStencilGeometry5D() {};
  virtual void BuildShifts(void)
  {
    this->shifts.resize(0);
    // Like HDCG:  81 point stencil including self connection
    this->shifts.push_back(Coordinate({0,0,0,0,0}));
    // +-x, +-y, +-z, +-t  : 8
    for(int s=-1;s<=1;s+=2){
      this->shifts.push_back(Coordinate({0,s,0,0,0}));
      this->shifts.push_back(Coordinate({0,0,s,0,0}));
      this->shifts.push_back(Coordinate({0,0,0,s,0}));
      this->shifts.push_back(Coordinate({0,0,0,0,s}));
    }
    // +-x+-y, +-x+-z, +-x+-t, +-y+-z, +-y+-t, +-z+-t : 24
    for(int s1=-1;s1<=1;s1+=2){
    for(int s2=-1;s2<=1;s2+=2){
      this->shifts.push_back(Coordinate({0,s1,s2,0,0}));
      this->shifts.push_back(Coordinate({0,s1,0,s2,0}));
      this->shifts.push_back(Coordinate({0,s1,0,0,s2}));
      this->shifts.push_back(Coordinate({0,0,s1,s2,0}));
      this->shifts.push_back(Coordinate({0,0,s1,0,s2}));
      this->shifts.push_back(Coordinate({0,0,0,s1,s2}));
    }}
    this->npoint = this->shifts.size();
  }
};
// Need to worry about red-black now
class NextToNextToNextToNearestStencilGeometry5D : public NonLocalStencilGeometry {
public:
  NextToNextToNextToNearestStencilGeometry5D(void) : NonLocalStencilGeometry(1)
  {
    this->BuildShifts();
  };
  virtual ~NextToNextToNextToNearestStencilGeometry5D() {}
  virtual void BuildShifts(void)
  {
    this->shifts.resize(0);
    // Like HDCG:  81 point stencil including self connection
    this->shifts.push_back(Coordinate({0,0,0,0,0}));
    // +-x, +-y, +-z, +-t  : 8
    for(int s=-1;s<=1;s+=2){
      this->shifts.push_back(Coordinate({0,s,0,0,0}));
      this->shifts.push_back(Coordinate({0,0,s,0,0}));
      this->shifts.push_back(Coordinate({0,0,0,s,0}));
      this->shifts.push_back(Coordinate({0,0,0,0,s}));
    }
    // +-x+-y, +-x+-z, +-x+-t, +-y+-z, +-y+-t, +-z+-t : 24
    for(int s1=-1;s1<=1;s1+=2){
    for(int s2=-1;s2<=1;s2+=2){
      this->shifts.push_back(Coordinate({0,s1,s2,0,0}));
      this->shifts.push_back(Coordinate({0,s1,0,s2,0}));
      this->shifts.push_back(Coordinate({0,s1,0,0,s2}));
      this->shifts.push_back(Coordinate({0,0,s1,s2,0}));
      this->shifts.push_back(Coordinate({0,0,s1,0,s2}));
      this->shifts.push_back(Coordinate({0,0,0,s1,s2}));
    }}
    // +-x+-y+-z, +-x+-y+-z, +-x+-y+-z,
    for(int s1=-1;s1<=1;s1+=2){
    for(int s2=-1;s2<=1;s2+=2){
    for(int s3=-1;s3<=1;s3+=2){
      this->shifts.push_back(Coordinate({0,s1,s2,s3,0})); // 8x4 = 32
      this->shifts.push_back(Coordinate({0,s1,s2,0,s3}));
      this->shifts.push_back(Coordinate({0,s1,0,s2,s3}));
      this->shifts.push_back(Coordinate({0,0,s1,s2,s3}));
    }}}
    for(int s1=-1;s1<=1;s1+=2){
    for(int s2=-1;s2<=1;s2+=2){
    for(int s3=-1;s3<=1;s3+=2){
    for(int s4=-1;s4<=1;s4+=2){
      this->shifts.push_back(Coordinate({0,s1,s2,s3,s4})); // 16 
    }}}}
    this->npoint = this->shifts.size();
  }
};

// Fine Object == (per site) type of fine field
// nbasis      == number of deflation vectors
template<class Fobj,class CComplex,int nbasis>
class GeneralCoarsenedMatrix : public SparseMatrixBase<Lattice<iVector<CComplex,nbasis > > >  {
public:
    
  typedef iVector<CComplex,nbasis >           siteVector;
  typedef Lattice<CComplex >                  CoarseComplexField;
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
    conformable(CoarseGrid(),in.Grid());
    conformable(in.Grid(),out.Grid());
    out.Checkerboard() = in.Checkerboard();
    CoarseVector tin=in;

    CoarseVector pin  = Cell.Exchange(tin);

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

    for(int point=0;point<npoint;point++){
      conformable(_A[point],pin);
    }
    
    // Should also exchange "A" and "Adag"
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

    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer[p].ViewClose();

    out = Cell.Extract(pout);
  };

  void Test(LinearOperatorBase<Lattice<Fobj> > &linop,
	    Aggregation<Fobj,CComplex,nbasis> & Subspace)
  {
    // Create a random
    GridBase *grid = FineGrid();
    FineField MbV(grid);
    FineField tmp(grid);
    FineField f_src(grid);
    FineField f_res(grid);
    FineField f_ref(grid);
    CoarseVector c_src(CoarseGrid());
    CoarseVector c_res(CoarseGrid());
    CoarseVector coarseInner(CoarseGrid());
    GridParallelRNG RNG(CoarseGrid());  RNG.SeedUniqueString(std::string("Coarse RNG"));
    random(RNG,c_src);
    blockPromote(c_src,f_src,Subspace.subspace);
    linop.op(f_src,f_ref);
    this->Mult (_A,c_src,c_res);
    blockPromote(c_res,f_res,Subspace.subspace);
    std::cout << " GeneralCoarsenedMatrix comparison res  "<<norm2(f_res)<<std::endl;
    std::cout << " GeneralCoarsenedMatrix comparison ref  "<<norm2(f_ref)<<std::endl;
    f_res = f_res - f_ref;
    std::cout << " GeneralCoarsenedMatrix comparison diff "<<norm2(f_res)<<std::endl;
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
    // Orthogonalise the subblocks over the basis
    CoarseScalar InnerProd(CoarseGrid()); 
    for(int b=0;b<nbasis;b++){
      std::cout << "subspace["<<b<<"] " <<norm2(Subspace.subspace[b])<<std::endl;
    }
    blockOrthogonalise(InnerProd,Subspace.subspace);

    // Now compute the matrix elements of linop between this orthonormal
    // set of vectors.
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
	  auto Adagb = peekSite(_Adag[p],bcoor);
	  for(int bb=0;bb<nbasis;bb++){
	    Ab(bb,b)    = ip(bb);
	    Adagb(b,bb) = conjugate(ip(bb));
	  }
	  pokeSite(Ab,_A[p],scoor);
	  pokeSite(Adagb,_Adag[p],bcoor);
	}
	tpeek+=usecond();
      }
    }
    for(int p=0;p<geom.npoint;p++){
      _A[p] = Cell.Exchange(_A[p]);
      _Adag[p] = Cell.Exchange(_Adag[p]);
    }
    std::cout << GridLogMessage<<"CoarsenOperator pick "<<tpick<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator mat  "<<tmat <<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator projection "<<tproj<<" us"<<std::endl;
    std::cout << GridLogMessage<<"CoarsenOperator peek/poke  "<<tpeek<<" us"<<std::endl;
  }

#if 0
  void CoarsenOperatorColoured(LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & Subspace)
  {
    std::cout << GridLogMessage<< "CoarsenMatrixColoured "<< std::endl;
    GridBase *grid = FineGrid();

    /////////////////////////////////////////////////////////////
    // Orthogonalise the subblocks over the basis
    /////////////////////////////////////////////////////////////
    CoarseScalar InnerProd(CoarseGrid()); 
    blockOrthogonalise(InnerProd,Subspace.subspace);

    const int npoint = geom.npoint
      
    // Now compute the matrix elements of linop between this orthonormal
    // set of vectors.
    FineField bV(grid);
    FineField MbV(grid);
    FineField tmp(grid);
    CoarseVector coarseInner(CoarseGrid());
    Coordinate clatt = CoarseGrid()->GlobalDimensions();
    Coordinate steclatt = CoarseGrid()->GlobalDimensions();
    for(int v=0;v<nbasis;v++){
      for(int p=0;p<npoint;p++){ // Loop over momenta
	/////////////////////////////////////////////////////
	// Stick a different phase on every block
	/////////////////////////////////////////////////////
	CoarseVector pha(CoarseGrid());
	std::vector<double> dmom = 
	for(int mu=0;mu<5;mu++){
	  dmom[mu] = imom[mu]*2*M_PI/global_nb[mu];
	}
	/*
	 *     Here, k,l index which possible momentum/shift within the N-points connected by MdagM.
	 *     Matrix index i is mapped to this shift via 
	 *               GetDelta(imat_ball_idx[i]).
	 *
	 *     conj(phases[block]) proj[k][ block*Nvec+j ] =  \sum_{l in ball}  e^{i q_k . delta_l} < phi_{block,j} | MdagM | phi_{(block+delta_l),i} > 
	 *                                                 =  \sum_{l in ball} e^{iqk.delta_l} A_ji^{b.b+l}
	 *                                                 = M_{kl} A_ji^{b.b+l}
	 *
	 *     Must assemble and invert matrix M_k,l = e^[i q_k . delta_l]
	 *  
	 *     Where q_k = delta_k . (2*M_PI/global_nb[mu])
	 *
	 *     Then A{ji}^{b,b+l} = M^{-1}_{lm} ComputeProj_{m,b,i,j}
	 */
      }
    }
    
    for(int p=0;p<geom.npoint;p++){
      _A[p] = Cell.Exchange(_A[p]);
      _Adag[p] = Cell.Exchange(_Adag[p]);
    }
  }
#endif

  virtual  void Mdiag    (const Field &in, Field &out){ assert(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};
};

NAMESPACE_END(Grid);
