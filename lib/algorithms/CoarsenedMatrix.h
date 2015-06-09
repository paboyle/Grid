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
    std::vector<int> opdirs;

    // FIXME -- don't like xposing the operator directions
    // as different to the geometrical dirs
    // Also don't like special casing five dim.. should pass an object in template
  Geometry(int _d)  {
  
      int base = (_d==5) ? 1:0;

      // make coarse grid stencil for 4d , not 5d
      if ( _d==5 ) _d=4;

      npoint = 2*_d+1;
      directions.resize(npoint);
      displacements.resize(npoint);
      opdirs.resize(npoint);
      for(int d=0;d<_d;d++){
	directions[2*d  ] = d+base;
	directions[2*d+1] = d+base;
	opdirs[2*d  ] = d;
	opdirs[2*d+1] = d;
	displacements[2*d  ] = +1;
	displacements[2*d+1] = -1;
      }
      directions   [2*_d]=0;
      displacements[2*_d]=0;
      opdirs       [2*_d]=0;
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
  
  // Fine Object == (per site) type of fine field
  // nbasis      == number of deflation vectors
  template<class Fobj,class CComplex,int nbasis>
  class CoarsenedMatrix : public SparseMatrixBase<Lattice<iVector<vComplex,nbasis > > >  {
  public:
    
    typedef iVector<vComplex,nbasis > siteVector;
    typedef Lattice<iVector<vComplex,nbasis > > CoarseVector;
    typedef Lattice<iMatrix<vComplex,nbasis > > CoarseMatrix;

    typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
    typedef Lattice<Fobj >        FineField;

    ////////////////////
    // Data members
    ////////////////////
    Geometry         geom;
    GridBase *       _grid; 
    CartesianStencil Stencil; 

    std::vector<CoarseMatrix> A;
    std::vector<CoarseMatrix> Aslow;

    std::vector<siteVector,alignedAllocator<siteVector> >   comm_buf;
      
    ///////////////////////
    // Interface
    ///////////////////////
    GridBase * Grid(void)         { return _grid; };   // this is all the linalg routines need to know

    RealD M    (const CoarseVector &in, CoarseVector &out){
      
      SimpleCompressor<siteVector> compressor;
      Stencil.HaloExchange(in,comm_buf,compressor);

PARALLEL_FOR_LOOP
      for(int ss=0;ss<Grid()->oSites();ss++){
        siteVector res = zero;
	siteVector tmp;
	siteVector nbr;

	int offset,local,perm;
	for(int point=0;point<geom.npoint;point++){
	  offset = Stencil._offsets [point][ss];
	  local  = Stencil._is_local[point][ss];
	  perm   = Stencil._permute[point][ss];
	  
	  if(local&&perm) { 
	    permute(nbr,in._odata[offset],perm);
	  } else if(local) { 
	    nbr = in._odata[offset];
	  } else {
	    nbr = comm_buf[offset];
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
	A(geom.npoint,&CoarseGrid),
	Aslow(geom.npoint,&CoarseGrid)
    {
      comm_buf.resize(Stencil._unified_buffer_size);
    };

    void CoarsenOperator(GridBase *FineGrid,LinearOperatorBase<Lattice<Fobj> > &linop,std::vector<Lattice<Fobj> > & subspace){

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
      blockOrthogonalise(InnerProd,subspace);
      blockProject(iProj,subspace[0],subspace);

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
	phi=subspace[i];
	for(int p=0;p<geom.npoint;p++){ 

	  int dir   = geom.directions[p];
	  int opdir = geom.opdirs[p];
	  int disp= geom.displacements[p];

	  int block=(FineGrid->_rdimensions[dir])/(Grid()->_rdimensions[dir]);

	  LatticeCoordinate(coor,dir);

	  if ( disp==0 ){
	    linop.OpDiag(phi,Mphi);
	  }
	  else  {
	    linop.OpDir(phi,Mphi,opdir,disp); 
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
	    oblock = where(mod(coor,block)==0,Mphi,zz);
	    iblock = where(mod(coor,block)!=0,Mphi,zz);
	  } else {
	    assert(0);
	  }

	  blockProject(iProj,iblock,subspace);
	  blockProject(oProj,oblock,subspace);
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
      std::cout<< " Computed matrix elements "<< self_stencil <<std::endl;
      for(int p=0;p<geom.npoint;p++){
	std::cout<< "A["<<p<<"]" << std::endl;
	std::cout<< A[p] << std::endl;
      }
      std::cout<< " picking by block0 "<< self_stencil <<std::endl;

      phi=subspace[0];
      std::vector<int> bc(FineGrid->_ndimension,0);

      blockPick(Grid(),phi,tmp,bc);      // Pick out a block
      linop.Op(tmp,Mphi);                // Apply big dop
      blockProject(iProj,Mphi,subspace); // project it and print it
      std::cout<< " Computed matrix elements from block zero only "<<std::endl;
      std::cout<< iProj <<std::endl;
      std::cout<<"Computed Coarse Operator"<<std::endl;
#endif
    }
    
  };

}
#endif
