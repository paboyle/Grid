    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/CoarsenedMatrix.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef  GRID_ALGORITHM_COARSENED_MATRIX_H
#define  GRID_ALGORITHM_COARSENED_MATRIX_H

#include <Grid.h>

namespace Grid {

  class Geometry {
    //    int dimension;
  public:
    int npoint;
    std::vector<int> directions   ;
    std::vector<int> displacements;

  Geometry(int _d)  {
  
      int base = (_d==5) ? 1:0;

      // make coarse grid stencil for 4d , not 5d
      if ( _d==5 ) _d=4;

      npoint = 2*_d+1;
      directions.resize(npoint);
      displacements.resize(npoint);
      for(int d=0;d<_d;d++){
	directions[2*d  ] = d+base;
	directions[2*d+1] = d+base;
	displacements[2*d  ] = +1;
	displacements[2*d+1] = -1;
      }
      directions   [2*_d]=0;
      displacements[2*_d]=0;
      
      //// report back
      std::cout<<GridLogMessage<<"directions    :";
      for(int d=0;d<npoint;d++) std::cout<< directions[d]<< " ";
      std::cout <<std::endl;
      std::cout<<GridLogMessage<<"displacements :";
      for(int d=0;d<npoint;d++) std::cout<< displacements[d]<< " ";
      std::cout<<std::endl;
    }
  
    /*
      // Original cleaner code
    Geometry(int _d) : dimension(_d), npoint(2*_d+1), directions(npoint), displacements(npoint) {
      for(int d=0;d<dimension;d++){
	directions[2*d  ] = d;
	directions[2*d+1] = d;
	displacements[2*d  ] = +1;
	displacements[2*d+1] = -1;
      }
      directions   [2*dimension]=0;
      displacements[2*dimension]=0;
    }
    std::vector<int> GetDelta(int point) {
      std::vector<int> delta(dimension,0);
      delta[directions[point]] = displacements[point];
      return delta;
    };
    */    

  };
  
  template<class Fobj,class CComplex,int nbasis>
  class Aggregation   {
  public:
    typedef iVector<CComplex,nbasis >             siteVector;
    typedef Lattice<siteVector>                 CoarseVector;
    typedef Lattice<iMatrix<CComplex,nbasis > > CoarseMatrix;

    typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
    typedef Lattice<Fobj >        FineField;

    GridBase *CoarseGrid;
    GridBase *FineGrid;
    std::vector<Lattice<Fobj> > subspace;

    Aggregation(GridBase *_CoarseGrid,GridBase *_FineGrid) : 
      CoarseGrid(_CoarseGrid),
      FineGrid(_FineGrid),
      subspace(nbasis,_FineGrid)
	{
	};
  
    void Orthogonalise(void){
      CoarseScalar InnerProd(CoarseGrid); 
      blockOrthogonalise(InnerProd,subspace);
    } 
    void CheckOrthogonal(void){
      CoarseVector iProj(CoarseGrid); 
      CoarseVector eProj(CoarseGrid); 
      Lattice<CComplex> pokey(CoarseGrid);

      
      for(int i=0;i<nbasis;i++){
	blockProject(iProj,subspace[i],subspace);

	eProj=zero; 
	for(int ss=0;ss<CoarseGrid->oSites();ss++){
	  eProj._odata[ss](i)=CComplex(1.0);
	}
	eProj=eProj - iProj;
	std::cout<<GridLogMessage<<"Orthog check error "<<i<<" " << norm2(eProj)<<std::endl;
      }
      std::cout<<GridLogMessage <<"CheckOrthog done"<<std::endl;
    }
    void ProjectToSubspace(CoarseVector &CoarseVec,const FineField &FineVec){
      blockProject(CoarseVec,FineVec,subspace);
    }
    void PromoteFromSubspace(const CoarseVector &CoarseVec,FineField &FineVec){
      blockPromote(CoarseVec,FineVec,subspace);
    }
    void CreateSubspaceRandom(GridParallelRNG &RNG){
      for(int i=0;i<nbasis;i++){
	random(RNG,subspace[i]);
	std::cout<<GridLogMessage<<" norm subspace["<<i<<"] "<<norm2(subspace[i])<<std::endl;
      }
      Orthogonalise();
    }

    virtual void CreateSubspaceLanczos(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,int nn=nbasis) 
    {
      // Run a Lanczos with sloppy convergence
	const int Nstop = nn;
	const int Nk = nn+20;
	const int Np = nn+20;
	const int Nm = Nk+Np;
	const int MaxIt= 10000;
	RealD resid = 1.0e-3;

	Chebyshev<FineField> Cheb(0.5,64.0,21);
	ImplicitlyRestartedLanczos<FineField> IRL(hermop,Cheb,Nstop,Nk,Nm,resid,MaxIt);
	//	IRL.lock = 1;

	FineField noise(FineGrid); gaussian(RNG,noise);
	FineField tmp(FineGrid); 
	std::vector<RealD>     eval(Nm);
	std::vector<FineField> evec(Nm,FineGrid);

	int Nconv;
	IRL.calc(eval,evec,
		 noise,
		 Nconv);

    	// pull back nn vectors
	for(int b=0;b<nn;b++){

	  subspace[b]   = evec[b];

	  std::cout << GridLogMessage <<"subspace["<<b<<"] = "<<norm2(subspace[b])<<std::endl;

	  hermop.Op(subspace[b],tmp); 
	  std::cout<<GridLogMessage << "filtered["<<b<<"] <f|MdagM|f> "<<norm2(tmp)<<std::endl;

	  noise = tmp -  sqrt(eval[b])*subspace[b] ;

	  std::cout<<GridLogMessage << " lambda_"<<b<<" = "<< eval[b] <<"  ;  [ M - Lambda ]_"<<b<<" vec_"<<b<<"  = " <<norm2(noise)<<std::endl;

	  noise = tmp +  eval[b]*subspace[b] ;

	  std::cout<<GridLogMessage << " lambda_"<<b<<" = "<< eval[b] <<"  ;  [ M - Lambda ]_"<<b<<" vec_"<<b<<"  = " <<norm2(noise)<<std::endl;

	}
	Orthogonalise();
	for(int b=0;b<nn;b++){
	  std::cout << GridLogMessage <<"subspace["<<b<<"] = "<<norm2(subspace[b])<<std::endl;
	}
    }

    virtual void CreateSubspace(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,int nn=nbasis) {

      RealD scale;

      ConjugateGradient<FineField> CG(1.0e-2,10000);
      FineField noise(FineGrid);
      FineField Mn(FineGrid);

      for(int b=0;b<nn;b++){
	
	gaussian(RNG,noise);
	scale = std::pow(norm2(noise),-0.5); 
	noise=noise*scale;

	hermop.Op(noise,Mn); std::cout<<GridLogMessage << "noise   ["<<b<<"] <n|MdagM|n> "<<norm2(Mn)<<std::endl;

	for(int i=0;i<1;i++){

	  CG(hermop,noise,subspace[b]);

	  noise = subspace[b];
	  scale = std::pow(norm2(noise),-0.5); 
	  noise=noise*scale;

	}

	hermop.Op(noise,Mn); std::cout<<GridLogMessage << "filtered["<<b<<"] <f|MdagM|f> "<<norm2(Mn)<<std::endl;
	subspace[b]   = noise;

      }

      Orthogonalise();

    }
  };
  // Fine Object == (per site) type of fine field
  // nbasis      == number of deflation vectors
  template<class Fobj,class CComplex,int nbasis>
  class CoarsenedMatrix : public SparseMatrixBase<Lattice<iVector<CComplex,nbasis > > >  {
  public:
    
    typedef iVector<CComplex,nbasis >             siteVector;
    typedef Lattice<siteVector>                 CoarseVector;
    typedef Lattice<iMatrix<CComplex,nbasis > > CoarseMatrix;

    typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
    typedef Lattice<Fobj >        FineField;

    ////////////////////
    // Data members
    ////////////////////
    Geometry         geom;
    GridBase *       _grid; 
    CartesianStencil<siteVector,siteVector> Stencil; 

    std::vector<CoarseMatrix> A;

      
    ///////////////////////
    // Interface
    ///////////////////////
    GridBase * Grid(void)         { return _grid; };   // this is all the linalg routines need to know

    RealD M (const CoarseVector &in, CoarseVector &out){

      conformable(_grid,in._grid);
      conformable(in._grid,out._grid);

      SimpleCompressor<siteVector> compressor;
      Stencil.HaloExchange(in,compressor);

PARALLEL_FOR_LOOP
      for(int ss=0;ss<Grid()->oSites();ss++){
        siteVector res = zero;
	siteVector nbr;
	int ptype;
	StencilEntry *SE;
	for(int point=0;point<geom.npoint;point++){

	  SE=Stencil.GetEntry(ptype,point,ss);
	  
	  if(SE->_is_local&&SE->_permute) { 
	    permute(nbr,in._odata[SE->_offset],ptype);
	  } else if(SE->_is_local) { 
	    nbr = in._odata[SE->_offset];
	  } else {
	    nbr = Stencil.comm_buf[SE->_offset];
	  }
	  res = res + A[point]._odata[ss]*nbr;
	}
	vstream(out._odata[ss],res);
      }
      return norm2(out);
    };

    RealD Mdag (const CoarseVector &in, CoarseVector &out){ 
      return M(in,out);
    };

    // Defer support for further coarsening for now
    void Mdiag    (const CoarseVector &in,  CoarseVector &out){};
    void Mdir     (const CoarseVector &in,  CoarseVector &out,int dir, int disp){};

    CoarsenedMatrix(GridCartesian &CoarseGrid) 	: 

      _grid(&CoarseGrid),
      geom(CoarseGrid._ndimension),
      Stencil(&CoarseGrid,geom.npoint,Even,geom.directions,geom.displacements),
      A(geom.npoint,&CoarseGrid)
    {
    };

    void CoarsenOperator(GridBase *FineGrid,LinearOperatorBase<Lattice<Fobj> > &linop,
			 Aggregation<Fobj,CComplex,nbasis> & Subspace){

      FineField iblock(FineGrid); // contributions from within this block
      FineField oblock(FineGrid); // contributions from outwith this block

      FineField     phi(FineGrid);
      FineField     tmp(FineGrid);
      FineField     zz(FineGrid); zz=zero;
      FineField    Mphi(FineGrid);

      Lattice<iScalar<vInteger> > coor(FineGrid);

      CoarseVector iProj(Grid()); 
      CoarseVector oProj(Grid()); 
      CoarseScalar InnerProd(Grid()); 

      // Orthogonalise the subblocks over the basis
      blockOrthogonalise(InnerProd,Subspace.subspace);

      // Compute the matrix elements of linop between this orthonormal
      // set of vectors.
      int self_stencil=-1;
      for(int p=0;p<geom.npoint;p++){ 
	A[p]=zero;
	if( geom.displacements[p]==0){
	  self_stencil=p;
	}
      }
      assert(self_stencil!=-1);

      for(int i=0;i<nbasis;i++){
	phi=Subspace.subspace[i];
	
	std::cout<<GridLogMessage<<"("<<i<<").."<<std::endl;

	for(int p=0;p<geom.npoint;p++){ 

	  int dir   = geom.directions[p];
	  int disp  = geom.displacements[p];

	  Integer block=(FineGrid->_rdimensions[dir])/(Grid()->_rdimensions[dir]);

	  LatticeCoordinate(coor,dir);

	  if ( disp==0 ){
	    linop.OpDiag(phi,Mphi);
	  }
	  else  {
	    linop.OpDir(phi,Mphi,dir,disp); 
	  }

	  ////////////////////////////////////////////////////////////////////////
	  // Pick out contributions coming from this cell and neighbour cell
	  ////////////////////////////////////////////////////////////////////////
	  if ( disp==0 ) {
	    iblock = Mphi;
	    oblock = zero;
	  } else if ( disp==1 ) {
	    oblock = where(mod(coor,block)==(block-1),Mphi,zz);
	    iblock = where(mod(coor,block)!=(block-1),Mphi,zz);
	  } else if ( disp==-1 ) {
	    oblock = where(mod(coor,block)==(Integer)0,Mphi,zz);
	    iblock = where(mod(coor,block)!=(Integer)0,Mphi,zz);
	  } else {
	    assert(0);
	  }

	  Subspace.ProjectToSubspace(iProj,iblock);
	  Subspace.ProjectToSubspace(oProj,oblock);
	  //	  blockProject(iProj,iblock,Subspace.subspace);
	  //	  blockProject(oProj,oblock,Subspace.subspace);
PARALLEL_FOR_LOOP
	  for(int ss=0;ss<Grid()->oSites();ss++){
	    for(int j=0;j<nbasis;j++){
	      if( disp!= 0 ) {
		A[p]._odata[ss](j,i) = oProj._odata[ss](j);
	      }
	      A[self_stencil]._odata[ss](j,i) =	A[self_stencil]._odata[ss](j,i) + iProj._odata[ss](j);
	    }
	  }
	}
      }

#if 0
      ///////////////////////////
      // test code worth preserving in if block
      ///////////////////////////
      std::cout<<GridLogMessage<< " Computed matrix elements "<< self_stencil <<std::endl;
      for(int p=0;p<geom.npoint;p++){
	std::cout<<GridLogMessage<< "A["<<p<<"]" << std::endl;
	std::cout<<GridLogMessage<< A[p] << std::endl;
      }
      std::cout<<GridLogMessage<< " picking by block0 "<< self_stencil <<std::endl;

      phi=Subspace.subspace[0];
      std::vector<int> bc(FineGrid->_ndimension,0);

      blockPick(Grid(),phi,tmp,bc);      // Pick out a block
      linop.Op(tmp,Mphi);                // Apply big dop
      blockProject(iProj,Mphi,Subspace.subspace); // project it and print it
      std::cout<<GridLogMessage<< " Computed matrix elements from block zero only "<<std::endl;
      std::cout<<GridLogMessage<< iProj <<std::endl;
      std::cout<<GridLogMessage<<"Computed Coarse Operator"<<std::endl;
#endif
      //      ForceHermitian();
      AssertHermitian();
      // ForceDiagonal();
    }
    void ForceDiagonal(void) {


      std::cout<<GridLogMessage<<"**************************************************"<<std::endl;
      std::cout<<GridLogMessage<<"****   Forcing coarse operator to be diagonal ****"<<std::endl;
      std::cout<<GridLogMessage<<"**************************************************"<<std::endl;
      for(int p=0;p<8;p++){
	A[p]=zero;
      }

      GridParallelRNG  RNG(Grid()); RNG.SeedRandomDevice();
      Lattice<iScalar<CComplex> > val(Grid()); random(RNG,val);

      Complex one(1.0);

      iMatrix<CComplex,nbasis> ident;  ident=one;

      val = val*adj(val);
      val = val + 1.0;

      A[8] = val*ident;

      //      for(int s=0;s<Grid()->oSites();s++) {
      //	A[8]._odata[s]=val._odata[s];
      //      }
    }
    void ForceHermitian(void) {
      for(int d=0;d<4;d++){
	int dd=d+1;
	A[2*d] = adj(Cshift(A[2*d+1],dd,1));
      }
      //      A[8] = 0.5*(A[8] + adj(A[8]));
    }
    void AssertHermitian(void) {
      CoarseMatrix AA    (Grid());
      CoarseMatrix AAc   (Grid());
      CoarseMatrix Diff  (Grid());
      for(int d=0;d<4;d++){
	
	int dd=d+1;
	AAc = Cshift(A[2*d+1],dd,1);
	AA  = A[2*d];
	
	Diff = AA - adj(AAc);

	std::cout<<GridLogMessage<<"Norm diff dim "<<d<<" "<< norm2(Diff)<<std::endl;
	std::cout<<GridLogMessage<<"Norm dim "<<d<<" "<< norm2(AA)<<std::endl;
	  
      }
      Diff = A[8] - adj(A[8]);
      std::cout<<GridLogMessage<<"Norm diff local "<< norm2(Diff)<<std::endl;
      std::cout<<GridLogMessage<<"Norm local "<< norm2(A[8])<<std::endl;
    }
    
  };

}
#endif
