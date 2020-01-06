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


NAMESPACE_BEGIN(Grid);

class Geometry {
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
    std::cout<<std::endl;
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
  int checkerboard;
  int Checkerboard(void){return checkerboard;}
  Aggregation(GridBase *_CoarseGrid,GridBase *_FineGrid,int _checkerboard) : 
    CoarseGrid(_CoarseGrid),
    FineGrid(_FineGrid),
    subspace(nbasis,_FineGrid),
    checkerboard(_checkerboard)
  {
  };
  
  void Orthogonalise(void){
    CoarseScalar InnerProd(CoarseGrid); 
    std::cout << GridLogMessage <<" Block Gramm-Schmidt pass 1"<<std::endl;
    blockOrthogonalise(InnerProd,subspace);
    //    std::cout << GridLogMessage <<" Block Gramm-Schmidt pass 2"<<std::endl; // Really have to do twice? Yuck
    //    blockOrthogonalise(InnerProd,subspace);
    //      std::cout << GridLogMessage <<" Gramm-Schmidt checking orthogonality"<<std::endl;
    //      CheckOrthogonal();
  } 
  void CheckOrthogonal(void){
    CoarseVector iProj(CoarseGrid); 
    CoarseVector eProj(CoarseGrid); 
    for(int i=0;i<nbasis;i++){
      blockProject(iProj,subspace[i],subspace);
      eProj=Zero(); 
      accelerator_for(ss, CoarseGrid->oSites(),1,{
	eProj[ss](i)=CComplex(1.0);
      });
      eProj=eProj - iProj;
      std::cout<<GridLogMessage<<"Orthog check error "<<i<<" " << norm2(eProj)<<std::endl;
    }
    std::cout<<GridLogMessage <<"CheckOrthog done"<<std::endl;
  }
  void ProjectToSubspace(CoarseVector &CoarseVec,const FineField &FineVec){
    //    std::cout << GridLogMessage<< "BlockPromote"<<std::endl;
    blockProject(CoarseVec,FineVec,subspace);
    //    std::cout << GridLogMessage<< "BlockPromote"<<std::endl;
  }
  void PromoteFromSubspace(const CoarseVector &CoarseVec,FineField &FineVec){
    FineVec.Checkerboard() = subspace[0].Checkerboard();
    //    std::cout << GridLogMessage<< "BlockPromote"<<std::endl;
    blockPromote(CoarseVec,FineVec,subspace);
    //    std::cout << GridLogMessage<< "BlockPromote done"<<std::endl;
  }
  void CreateSubspaceRandom(GridParallelRNG &RNG){
    for(int i=0;i<nbasis;i++){
      random(RNG,subspace[i]);
      //      std::cout<<GridLogMessage<<" norm subspace["<<i<<"] "<<norm2(subspace[i])<<std::endl;
    }
  }

  virtual void CreateSubspace(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,int nn=nbasis) {

    RealD scale;

    ConjugateGradient<FineField> CG(1.0e-2,100,false);
    FineField noise(FineGrid);
    FineField Mn(FineGrid);

    for(int b=0;b<nn;b++){
	
      subspace[b] = Zero();
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
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // World of possibilities here. But have tried quite a lot of experiments (250+ jobs run on Summit)
  // and this is the best I found
  ////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void CreateSubspaceChebyshev(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,
				       int nn,
				       double hi,
				       double lo,
				       int order,
				       int orderstep
				       ) {

    RealD scale;

    FineField noise(FineGrid);
    FineField Mn(FineGrid);
    FineField tmp(FineGrid);

    // New normalised noise
    gaussian(RNG,noise);
    scale = std::pow(norm2(noise),-0.5); 
    noise=noise*scale;

    // Initial matrix element
    hermop.Op(noise,Mn); std::cout<<GridLogMessage << "noise <n|MdagM|n> "<<norm2(Mn)<<std::endl;

    int b =0;
    {
      // Filter
      Chebyshev<FineField> Cheb(lo,hi,order);
      Cheb(hermop,noise,Mn);
      // normalise
      scale = std::pow(norm2(Mn),-0.5); 	Mn=Mn*scale;
      subspace[b]   = Mn;
      hermop.Op(Mn,tmp); 
      std::cout<<GridLogMessage << "filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;
      b++;
    }

    // Generate a full sequence of Chebyshevs
    {
      lo=0;
      noise=Mn;

      FineField T0(FineGrid); T0 = noise;  
      FineField T1(FineGrid); 
      FineField T2(FineGrid);
      FineField y(FineGrid);
      
      FineField *Tnm = &T0;
      FineField *Tn  = &T1;
      FineField *Tnp = &T2;

      // Tn=T1 = (xscale M + mscale)in
      RealD xscale = 2.0/(hi-lo);
      RealD mscale = -(hi+lo)/(hi-lo);
      hermop.HermOp(T0,y);
      T1=y*xscale+noise*mscale;

      for(int n=2;n<=orderstep*(nn-1);n++){
	
	hermop.HermOp(*Tn,y);

	auto y_v = y.View();
	auto Tn_v = Tn->View();
	auto Tnp_v = Tnp->View();
	auto Tnm_v = Tnm->View();
	accelerator_forNB(ss, FineGrid->oSites(), Nsimd, {
	  coalescedWrite(y_v[ss],xscale*y_v(ss)+mscale*Tn_v(ss));
	  coalescedWrite(Tnp_v[ss],2.0*y_v(ss)-Tnm_v(ss));
        });

	// Possible more fine grained control is needed than a linear sweep,
	// but huge productivity gain if this is simple algorithm and not a tunable
	if ( (n%orderstep)==0 ) { 
	  Mn=*Tnp;
	  scale = std::pow(norm2(Mn),-0.5);         Mn=Mn*scale;
	  subspace[b] = Mn;
	  hermop.Op(Mn,tmp); 
	  std::cout<<GridLogMessage << "filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;
	  b++;
	}

	// Cycle pointers to avoid copies
	FineField *swizzle = Tnm;
	Tnm    =Tn;
	Tn     =Tnp;
	Tnp    =swizzle;
	  
      }
    }
    assert(b==nn);
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
  typedef iMatrix<CComplex,nbasis >  Cobj;
  typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj >        FineField;

  ////////////////////
  // Data members
  ////////////////////
  Geometry         geom;
  GridBase *       _grid; 
  int hermitian;

  CartesianStencil<siteVector,siteVector,int> Stencil; 

  std::vector<CoarseMatrix> A;

      
  ///////////////////////
  // Interface
  ///////////////////////
  GridBase * Grid(void)         { return _grid; };   // this is all the linalg routines need to know

  RealD M (const CoarseVector &in, CoarseVector &out){

    conformable(_grid,in.Grid());
    conformable(in.Grid(),out.Grid());


    //    RealD Nin = norm2(in);
    SimpleCompressor<siteVector> compressor;

    double comms_usec = -usecond();
    Stencil.HaloExchange(in,compressor);
    comms_usec += usecond();

    auto in_v = in.View();
    auto out_v = out.View();
    typedef LatticeView<Cobj> Aview;

    Vector<Aview> AcceleratorViewContainer;
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(A[p].View());
    Aview *Aview_p = & AcceleratorViewContainer[0];

    const int Nsimd = CComplex::Nsimd();

    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    GridStopWatch ArithmeticTimer;
    int osites=Grid()->oSites();
    double flops = osites*Nsimd*nbasis*nbasis*8.0*geom.npoint;
    double bytes = osites*nbasis*nbasis*geom.npoint*sizeof(CComplex);
    double usecs =-usecond();

    // assert(geom.npoint==9);

    accelerator_for(sss, Grid()->oSites()*nbasis, Nsimd, {
      int ss = sss/nbasis;
      int b  = sss%nbasis;
      calcComplex res = Zero();
      calcVector nbr;
      int ptype;
      StencilEntry *SE;

      int lane=SIMTlane(Nsimd);
      for(int point=0;point<geom.npoint;point++){

	SE=Stencil.GetEntry(ptype,point,ss);
	  
	if(SE->_is_local) { 
	  nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute,lane);
	} else {
	  nbr = coalescedRead(Stencil.CommBuf()[SE->_offset],lane);
	}
	synchronise();

	for(int bb=0;bb<nbasis;bb++) {
	  res = res + coalescedRead(Aview_p[point][ss](b,bb))*nbr(bb);
	}
      }
      coalescedWrite(out_v[ss](b),res,lane);
    });
    usecs +=usecond();

    double nrm_usec=-usecond();
    RealD Nout= norm2(out);
    nrm_usec+=usecond();

    /*
        std::cout << GridLogMessage << "\tNorm        " << nrm_usec << " us" <<std::endl;
        std::cout << GridLogMessage << "\tHalo        " << comms_usec << " us" <<std::endl;
        std::cout << GridLogMessage << "\tMatrix      " << usecs << " us" <<std::endl;
        std::cout << GridLogMessage << "\t  mflop/s   " << flops/usecs<<std::endl;
        std::cout << GridLogMessage << "\t  MB/s      " << bytes/usecs<<std::endl;
    */
    return Nout;
  };

  RealD Mdag (const CoarseVector &in, CoarseVector &out)
  {
    if(hermitian) {
      // corresponds to Petrov-Galerkin coarsening
      return M(in,out);
    } else {
      // corresponds to Galerkin coarsening
      CoarseVector tmp(Grid());
      G5C(tmp, in); 
      M(tmp, out);
      G5C(out, out);
      return norm2(out);
    }
  };
  void MdirComms(const CoarseVector &in)
  {
    SimpleCompressor<siteVector> compressor;
    Stencil.HaloExchange(in,compressor);
  }
  void MdirCalc(const CoarseVector &in, CoarseVector &out, int point)
  {
    conformable(_grid,in.Grid());
    conformable(_grid,out.Grid());

    typedef LatticeView<Cobj> Aview;
    Vector<Aview> AcceleratorViewContainer;
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(A[p].View());
    Aview *Aview_p = & AcceleratorViewContainer[0];

    auto out_v = out.View();
    auto in_v  = in.View();
    accelerator_for(ss,Grid()->oSites(),1,{
      siteVector res = Zero();
      siteVector nbr;
      int ptype;
      StencilEntry *SE;
      
      SE=Stencil.GetEntry(ptype,point,ss);
      
      if(SE->_is_local&&SE->_permute) {
	permute(nbr,in_v[SE->_offset],ptype);
      } else if(SE->_is_local) {
	nbr = in_v[SE->_offset];
      } else {
	nbr = Stencil.CommBuf()[SE->_offset];
      }

      res = res + Aview_p[point][ss]*nbr;
      
      out_v[ss]=res;
    });

  }
  void MdirAll(const CoarseVector &in,std::vector<CoarseVector> &out)
  {
    this->MdirComms(in);
    int ndir=geom.npoint-1;
    assert(out.size()==ndir);
    for(int p=0;p<ndir;p++){
      MdirCalc(in,out[p],p);
    }
  };
  void Mdir(const CoarseVector &in, CoarseVector &out, int dir, int disp){

    this->MdirComms(in);

    int ndim = in.Grid()->Nd();

    //////////////
    // 4D action like wilson
    // 0+ => 0 
    // 0- => 1
    // 1+ => 2 
    // 1- => 3
    // etc..
    //////////////
    // 5D action like DWF
    // 1+ => 0 
    // 1- => 1
    // 2+ => 2 
    // 2- => 3
    // etc..
    auto point = [dir, disp, ndim](){
      if(dir == 0 and disp == 0)
	return 8;
      else if ( ndim==4 ) { 
	return (4 * dir + 1 - disp) / 2;
      } else { 
	return (4 * (dir-1) + 1 - disp) / 2;
      }
    }();

    MdirCalc(in,out,point);

  };

  void Mdiag(const CoarseVector &in, CoarseVector &out)
  {
    int point=geom.npoint-1;
    MdirCalc(in, out, point); // No comms
  };

  
 CoarsenedMatrix(GridCartesian &CoarseGrid, int hermitian_=0) 	: 

    _grid(&CoarseGrid),
    geom(CoarseGrid._ndimension),
    hermitian(hermitian_),
    Stencil(&CoarseGrid,geom.npoint,Even,geom.directions,geom.displacements,0),
    A(geom.npoint,&CoarseGrid)
  {
  };

  void CoarsenOperator(GridBase *FineGrid,LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & Subspace){

    FineField iblock(FineGrid); // contributions from within this block
    FineField oblock(FineGrid); // contributions from outwith this block

    FineField     phi(FineGrid);
    FineField     tmp(FineGrid);
    FineField     zz(FineGrid); zz=Zero();
    FineField    Mphi(FineGrid);
    std::vector<FineField>     Mphi_p(geom.npoint,FineGrid);

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
      A[p]=Zero();
      if( geom.displacements[p]==0){
	self_stencil=p;
      }
    }
    assert(self_stencil!=-1);

    for(int i=0;i<nbasis;i++){
      phi=Subspace.subspace[i];

      std::cout << GridLogMessage<< "CoarsenMatrix vector "<<i<<" OpDir " << std::endl;
      linop.OpDirAll(phi,Mphi_p);
      std::cout << GridLogMessage<< "CoarsenMatrix vector "<<i<<" OpDir calculated" << std::endl;
      linop.OpDiag  (phi,Mphi_p[geom.npoint-1]);
      std::cout << GridLogMessage<< "CoarsenMatrix vector "<<i<<" OpDiag calculated" << std::endl;

      for(int p=0;p<geom.npoint;p++){ 

	Mphi = Mphi_p[p];

	int dir   = geom.directions[p];
	int disp  = geom.displacements[p];

	Integer block=(FineGrid->_rdimensions[dir])/(Grid()->_rdimensions[dir]);

	LatticeCoordinate(coor,dir);


	////////////////////////////////////////////////////////////////////////
	// Pick out contributions coming from this cell and neighbour cell
	////////////////////////////////////////////////////////////////////////
	if ( disp==0 ) {
	  iblock = Mphi;
	  oblock = Zero();
	} else if ( disp==1 ) {
	  oblock = where(mod(coor,block)==(block-1),Mphi,zz);
	  iblock = where(mod(coor,block)!=(block-1),Mphi,zz);
	} else if ( disp==-1 ) {
	  oblock = where(mod(coor,block)==(Integer)0,Mphi,zz);
	  iblock = where(mod(coor,block)!=(Integer)0,Mphi,zz);
	} else {
	  assert(0);
	}

	// Could do local inner products,
	// and then block pick the IP's.
	// Ideally write a routine to do two masked block sums at once
	std::cout << GridLogMessage<< "CoarsenMatrix picked "<<p<< std::endl;
	Subspace.ProjectToSubspace(iProj,iblock);
	Subspace.ProjectToSubspace(oProj,oblock);
	std::cout << GridLogMessage<< "CoarsenMatrix projected"<<p<< std::endl;

	// 4x gain possible in this loop. Profile and identify time loss.
	// i)  Assume Hermiticity, upper diagonal only (2x)
	// ii) Local inner product, then pick the local inners and sum. (2x)

	auto iProj_v = iProj.View() ;
	auto oProj_v = oProj.View() ;
	auto A_p     =  A[p].View();
	auto A_self  = A[self_stencil].View();
	accelerator_for(ss, Grid()->oSites(),1,{
	  for(int j=0;j<nbasis;j++){
	    if( disp!= 0 ) {
	      A_p[ss](j,i) = oProj_v[ss](j);
	    }
	    A_self[ss](j,i) =	A_self[ss](j,i) + iProj_v[ss](j);
	  }
	});

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
    /*
    if(hermitian) {
      std::cout << GridLogMessage << " ForceHermitian "<<std::endl;
      ForceHermitian();
    }
    for(int p=0;p<geom.npoint;p++){
      std::cout << GridLogMessage<< " dir "<< norm2(A[p]) <<std::endl;
    }
    */
      // AssertHermitian();
      // ForceDiagonal();
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

NAMESPACE_END(Grid);
#endif
