#ifndef  GRID_ALGORITHM_COARSENED_MATRIX_H
#define  GRID_ALGORITHM_COARSENED_MATRIX_H

#include <Grid.h>

namespace Grid {

  class Geometry {
  public:
    int npoint;
    int dimension;
    std::vector<int> directions   ;
    std::vector<int> displacements;

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
    std::vector<siteVector,alignedAllocator<siteVector> >   comm_buf;
      
    ///////////////////////
    // Interface
    ///////////////////////
    GridBase * Grid(void)         { return _grid; };   // this is all the linalg routines need to know

    RealD M    (const CoarseVector &in, CoarseVector &out){
      
      SimpleCompressor<siteVector> compressor;
      Stencil.HaloExchange(in,comm_buf,compressor);

      //PARALLEL_FOR_LOOP
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
      A(geom.npoint,&CoarseGrid)
    {
      comm_buf.resize(Stencil._unified_buffer_size);
    };

    void CoarsenOperator(GridBase *FineGrid,LinearOperatorBase<Lattice<Fobj> > &linop,std::vector<Lattice<Fobj> > & subspace){

      FineField  phi(FineGrid);
      FineField Mphi(FineGrid);
      CoarseVector Proj(Grid()); 
      CoarseScalar InnerProd(Grid()); 

      // Orthogonalise the subblocks over the basis
      blockOrthogonalise(InnerProd,subspace);

      // Compute the matrix elements of linop between this orthonormal
      // set of vectors.
      for(int i=0;i<nbasis;i++){
	phi=subspace[i];
	for(int p=0;p<geom.npoint;p++){ 
	  int dir = geom.directions[p];
	  int disp= geom.displacements[p];

	  if ( disp==0 )linop.OpDiag(phi,Mphi);
	  else linop.OpDir(phi,Mphi,dir,disp); 

	  blockProject(Proj,Mphi,subspace);

	  for(int ss=0;ss<Grid()->oSites();ss++){
	    for(int j=0;j<nbasis;j++){
	      A[p]._odata[ss](j,i) = Proj._odata[ss](j);
	    }
	  }

	}	
      }
      std::cout<<"Computed Coarse Operator"<<std::endl;
    }
    
  };

}
#endif
