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
  NextToNearestStencilGeometry4D(void) : NonLocalStencilGeometry(2)
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
  NextToNextToNextToNearestStencilGeometry4D(void) : NonLocalStencilGeometry(4)
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
  NextToNearestStencilGeometry5D(void) : NonLocalStencilGeometry(2)
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
  NextToNextToNextToNearestStencilGeometry5D(void) : NonLocalStencilGeometry(4)
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
  GridCartesian *       _FineGrid; 
  GridCartesian *       _CoarseGrid; 
  NonLocalStencilGeometry &geom;
  PaddedCell Cell;
  GeneralLocalStencil Stencil;
  
  std::vector<CoarseMatrix> _A;
  std::vector<CoarseMatrix> _Adag;

  ///////////////////////
  // Interface
  ///////////////////////
  GridCartesian * Grid(void)           { return _FineGrid; };   // this is all the linalg routines need to know
  GridCartesian * FineGrid(void)       { return _FineGrid; };   // this is all the linalg routines need to know
  GridCartesian * CoarseGrid(void)     { return _CoarseGrid; };   // this is all the linalg routines need to know

  GeneralCoarsenedMatrix(NonLocalStencilGeometry &_geom,GridCartesian *FineGrid, GridCartesian * CoarseGrid)
    : geom(_geom),
      _FineGrid(FineGrid),
      _CoarseGrid(CoarseGrid),
      hermitian(1),
      Cell(_geom.Depth(),_CoarseGrid),
      Stencil(Cell.grids.back(),geom.shifts)
  {
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
    std::cout << "Calling Exchange"<<std::endl;
    CoarseVector pin  = Cell.Exchange(tin);
    //    std::cout << "Called Exchange"<<std::endl;
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

      // FIXME -- exchange the A and the A dag
      
      for(int point=0;point<npoint;point++){

	auto SE = Stencil_v.GetEntry(point,ss);
	
	int o = SE->_offset;

	// gpermute etc..
	nbr = in_v[o];
	assert( o< osites);
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
    GridCartesian *grid = FineGrid();
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
    std::cout << GridLogMessage<< "CoarsenMatrix "<< std::endl;
    GridCartesian *grid = FineGrid();
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
      std::cout << GridLogMessage<< "CoarsenMatrix block "<< bcoor << std::endl;
      for(int b=0;b<nbasis;b++){
	blockPick(CoarseGrid(),Subspace.subspace[b],bV,bcoor);
	linop.HermOp(bV,MbV);
	blockProject(coarseInner,MbV,Subspace.subspace);
	for(int p=0;p<geom.npoint;p++){
	  Coordinate scoor = bcoor;
	  for(int mu=0;mu<bcoor.size();mu++){
	    int L = CoarseGrid()->GlobalDimensions()[mu];
	    scoor[mu] = (bcoor[mu] - geom.shifts[p][mu] + L) % L; // Modulo arithmetic
	  }
	  auto ip    = peekSite(coarseInner,scoor);
	  std::cout << "A["<<b<<"]["<<p<<"]"<<scoor<<" "<<" shift "<<geom.shifts[p]<<" "<< ip <<std::endl;
	  auto Ab    = peekSite(_A[p],scoor);
	  auto Adagb = peekSite(_Adag[p],bcoor);
	  for(int bb=0;bb<nbasis;bb++){
	    Ab(bb,b)    = ip(bb);
	    Adagb(b,bb) = conjugate(ip(bb));
	  }
	  pokeSite(Ab,_A[p],scoor);
	  pokeSite(Adagb,_Adag[p],bcoor);
	}
      }
    }
    std::cout << " Exchanging _A " <<std::endl;
    for(int p=0;p<geom.npoint;p++){
      _A[p] = Cell.Exchange(_A[p]);
      _Adag[p] = Cell.Exchange(_Adag[p]);
    }
  }
  virtual  void Mdiag    (const Field &in, Field &out){ assert(0);};
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};
};

NAMESPACE_END(Grid);
