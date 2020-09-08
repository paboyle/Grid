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

template<class vobj,class CComplex>
inline void blockMaskedInnerProduct(Lattice<CComplex> &CoarseInner,
				    const Lattice<decltype(innerProduct(vobj(),vobj()))> &FineMask,
				    const Lattice<vobj> &fineX,
				    const Lattice<vobj> &fineY)
{
  typedef decltype(innerProduct(vobj(),vobj())) dotp;

  GridBase *coarse(CoarseInner.Grid());
  GridBase *fine  (fineX.Grid());

  Lattice<dotp> fine_inner(fine); fine_inner.Checkerboard() = fineX.Checkerboard();
  Lattice<dotp> fine_inner_msk(fine);

  // Multiply could be fused with innerProduct
  fine_inner = localInnerProduct(fineX,fineY);
  mult(fine_inner_msk, fine_inner,FineMask);
  blockSum(CoarseInner,fine_inner_msk);
}

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
      directions[d   ] = d+base;
      directions[d+_d] = d+base;
      displacements[d  ] = +1;
      displacements[d+_d]= -1;
    }
    directions   [2*_d]=0;
    displacements[2*_d]=0;

    std::cout <<GridLogMessage << "Geometry "<<std::endl;
    for(int p=0;p<npoint;p++){
      std::cout <<GridLogMessage << "point " <<p<<" dir "<<directions[p]<<" delta " <<displacements[p]<<std::endl;
    }
  }
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
    checkerboard(_checkerboard),
    subspace(nbasis,_FineGrid)
  {
  };
  
  void Orthogonalise(void){
    CoarseScalar InnerProd(CoarseGrid); 
    std::cout << GridLogMessage <<" Block Gramm-Schmidt pass 1"<<std::endl;
    blockOrthogonalise(InnerProd,subspace);
  } 
  void ProjectToSubspace(CoarseVector &CoarseVec,const FineField &FineVec){
    blockProject(CoarseVec,FineVec,subspace);
  }
  void PromoteFromSubspace(const CoarseVector &CoarseVec,FineField &FineVec){
    FineVec.Checkerboard() = subspace[0].Checkerboard();
    blockPromote(CoarseVec,FineVec,subspace);
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
				       int orderfilter,
				       int ordermin,
				       int orderstep,
				       double filterlo
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
      Chebyshev<FineField> Cheb(lo,hi,orderfilter);
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
      lo=filterlo;
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

      for(int n=2;n<=ordermin+orderstep*(nn-2);n++){
	
	hermop.HermOp(*Tn,y);

	autoView( y_v , y, AcceleratorWrite);
	autoView( Tn_v , (*Tn), AcceleratorWrite);
	autoView( Tnp_v , (*Tnp), AcceleratorWrite);
	autoView( Tnm_v , (*Tnm), AcceleratorWrite);
	const int Nsimd = CComplex::Nsimd();
	accelerator_forNB(ss, FineGrid->oSites(), Nsimd, {
	  coalescedWrite(y_v[ss],xscale*y_v(ss)+mscale*Tn_v(ss));
	  coalescedWrite(Tnp_v[ss],2.0*y_v(ss)-Tnm_v(ss));
        });

	// Possible more fine grained control is needed than a linear sweep,
	// but huge productivity gain if this is simple algorithm and not a tunable
	int m =1;
	if ( n>=ordermin ) m=n-ordermin;
	if ( (m%orderstep)==0 ) { 
	  Mn=*Tnp;
	  scale = std::pow(norm2(Mn),-0.5);         Mn=Mn*scale;
	  subspace[b] = Mn;
	  hermop.Op(Mn,tmp); 
	  std::cout<<GridLogMessage << n<<" filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;
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
    
  typedef iVector<CComplex,nbasis >           siteVector;
  typedef Lattice<CComplex >                  CoarseComplexField;
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

  virtual std::vector<int> Directions(void)   { return geom.directions; };
  virtual std::vector<int> Displacements(void){ return geom.displacements; };
  void M (const CoarseVector &in, CoarseVector &out)
  {
    conformable(_grid,in.Grid());
    conformable(in.Grid(),out.Grid());

    SimpleCompressor<siteVector> compressor;

    Stencil.HaloExchange(in,compressor);
    autoView( in_v , in, AcceleratorRead);
    autoView( out_v , out, AcceleratorWrite);
    typedef LatticeView<Cobj> Aview;
      
    Vector<Aview> AcceleratorViewContainer;
  
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(A[p].View(AcceleratorRead));
    Aview *Aview_p = & AcceleratorViewContainer[0];

    const int Nsimd = CComplex::Nsimd();
    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    int osites=Grid()->oSites();

    autoView(st,Stencil,AcceleratorRead);
    siteVector *CBp=Stencil.CommBuf();

    accelerator_for(sss, Grid()->oSites()*nbasis, Nsimd, {
      int ss = sss/nbasis;
      int b  = sss%nbasis;
      calcComplex res = Zero();
      calcVector nbr;
      int ptype;
      StencilEntry *SE;

      for(int point=0;point<geom.npoint;point++){

	SE=st.GetEntry(ptype,point,ss);
	  
	if(SE->_is_local) { 
	  nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
	} else {
	  nbr = coalescedRead(CBp[SE->_offset]);
	}
	acceleratorSynchronise();

	for(int bb=0;bb<nbasis;bb++) {
	  res = res + coalescedRead(Aview_p[point][ss](b,bb))*nbr(bb);
	}
      }
      coalescedWrite(out_v[ss](b),res);
    });

    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer[p].ViewClose();
  };

  void Mdag (const CoarseVector &in, CoarseVector &out)
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
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(A[p].View(AcceleratorRead));
    Aview *Aview_p = & AcceleratorViewContainer[0];

    autoView( out_v , out, AcceleratorWrite);
    autoView( in_v  , in, AcceleratorRead);

    const int Nsimd = CComplex::Nsimd();
    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    accelerator_for(sss, Grid()->oSites()*nbasis, Nsimd, {
      int ss = sss/nbasis;
      int b  = sss%nbasis;
      calcComplex res = Zero();
      calcVector nbr;
      int ptype;
      StencilEntry *SE;

      SE=Stencil.GetEntry(ptype,point,ss);
	  
      if(SE->_is_local) { 
	nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
      } else {
	nbr = coalescedRead(Stencil.CommBuf()[SE->_offset]);
      }
      acceleratorSynchronise();

      for(int bb=0;bb<nbasis;bb++) {
	res = res + coalescedRead(Aview_p[point][ss](b,bb))*nbr(bb);
      }
      coalescedWrite(out_v[ss](b),res);
    });
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer[p].ViewClose();
  }
  void MdirAll(const CoarseVector &in,std::vector<CoarseVector> &out)
  {
    this->MdirComms(in);
    int ndir=geom.npoint-1;
    if ((out.size()!=ndir)&&(out.size()!=ndir+1)) { 
      std::cout <<"MdirAll out size "<< out.size()<<std::endl;
      std::cout <<"MdirAll ndir "<< ndir<<std::endl;
      assert(0);
    }
    for(int p=0;p<ndir;p++){
      MdirCalc(in,out[p],p);
    }
  };
  void Mdir(const CoarseVector &in, CoarseVector &out, int dir, int disp)
  {
    this->MdirComms(in);

    int ndim = in.Grid()->Nd();

    int point=-1;
    for(int p=0;p<geom.npoint;p++){
      if( (dir==geom.directions[p])&&(disp==geom.displacements[p])) point=p;
    }
    assert(point!=-1);// Must find

    std::cout <<GridLogMessage << "Mdir point "<<point<<" dir "<<dir<<" disp "<<disp  <<std::endl;
    for(int p=0;p<geom.npoint;p++){
      std::cout <<GridLogMessage << "point " <<p<<" dir "<<geom.directions[p]<<" delta " <<geom.displacements[p]<<std::endl;
    }
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

  void Test(Aggregation<Fobj,CComplex,nbasis> &_Aggregates,GridBase *FineGrid,LinearOperatorBase<Lattice<Fobj> > &linop)
  {
    typedef Lattice<Fobj> FineField;
    CoarseVector Cin(_grid);
    CoarseVector Cout(_grid);
    CoarseVector CFout(_grid);

    FineField Fin(FineGrid);
    FineField Fout(FineGrid);


    std::vector<int> seeds({1,2,3,4,5});
    GridParallelRNG RNG(_grid);  RNG.SeedFixedIntegers(seeds);
    gaussian(RNG,Cin);

    _Aggregates.PromoteFromSubspace(Cin,Fin);
    _Aggregates.ProjectToSubspace(Cin,Fin);

    std::cout << GridLogMessage<< "************  "<<std::endl;
    std::cout << GridLogMessage<< " Testing M  "<<std::endl;
    std::cout << GridLogMessage<< "************  "<<std::endl;
    // Coarse operator
    this->M(Cin,Cout);
    // Fine projected operator
    _Aggregates.PromoteFromSubspace(Cin,Fin);
    linop.Op(Fin,Fout);
    _Aggregates.ProjectToSubspace(CFout,Fout);

    CFout = CFout-Cout;
    RealD diff = norm2(CFout);
    std::cout << GridLogMessage<< " diff  "<<diff<<std::endl;
    assert(diff<1.0e-5);

    std::cout << GridLogMessage<< "************  "<<std::endl;
    std::cout << GridLogMessage<< " Testing Mdag  "<<std::endl;
    std::cout << GridLogMessage<< "************  "<<std::endl;
    // Coarse operator
    Mdag(Cin,Cout);
    // Fine operator
    linop.AdjOp(Fin,Fout);
    _Aggregates.ProjectToSubspace(CFout,Fout);

    CFout = CFout-Cout;
    diff = norm2(CFout);
    std::cout << GridLogMessage<< " diff  "<<diff<<std::endl; 
    assert(diff<1.0e-5);
  }

  void CoarsenOperator(GridBase *FineGrid,LinearOperatorBase<Lattice<Fobj> > &linop,
		       Aggregation<Fobj,CComplex,nbasis> & Subspace)
  {
    typedef Lattice<typename Fobj::tensor_reduced> FineComplexField;
    typedef typename Fobj::scalar_type scalar_type;

    FineComplexField one(FineGrid); one=scalar_type(1.0,0.0);
    FineComplexField zero(FineGrid); zero=scalar_type(0.0,0.0);

    std::vector<FineComplexField> masks(geom.npoint,FineGrid);
    FineComplexField imask(FineGrid); // contributions from within this block
    FineComplexField omask(FineGrid); // contributions from outwith this block

    FineComplexField evenmask(FineGrid);
    FineComplexField oddmask(FineGrid); 

    FineField     phi(FineGrid);
    FineField     tmp(FineGrid);
    FineField     zz(FineGrid); zz=Zero();
    FineField    Mphi(FineGrid);
    FineField    Mphie(FineGrid);
    FineField    Mphio(FineGrid);
    std::vector<FineField>     Mphi_p(geom.npoint,FineGrid);

    Lattice<iScalar<vInteger> > coor (FineGrid);
    Lattice<iScalar<vInteger> > bcoor(FineGrid);
    Lattice<iScalar<vInteger> > bcb  (FineGrid); bcb = Zero();

    CoarseVector iProj(Grid()); 
    CoarseVector oProj(Grid()); 
    CoarseVector SelfProj(Grid()); 
    CoarseComplexField iZProj(Grid()); 
    CoarseComplexField oZProj(Grid()); 

    CoarseScalar InnerProd(Grid()); 

    std::cout << GridLogMessage<< "CoarsenMatrix Orthog " << std::endl;
    // Orthogonalise the subblocks over the basis
    blockOrthogonalise(InnerProd,Subspace.subspace);
    std::cout << GridLogMessage<< "CoarsenMatrix Orthog done " << std::endl;

    auto OpDirections    = linop.Directions();
    auto OpDisplacements = linop.Displacements();

    std::cout<<" Coarsening an operator with "<< OpDirections.size()<<" terms "<<std::endl;
    for(int p=0;p<OpDirections.size();p++) {
      assert(OpDirections[p]==geom.directions[p]);
      assert(OpDisplacements[p]==geom.displacements[p]);
    }

    // Compute the matrix elements of linop between this orthonormal
    // set of vectors.
    int self_stencil=-1;
    for(int p=0;p<geom.npoint;p++)
    { 
      int dir   = geom.directions[p];
      int disp  = geom.displacements[p];
      A[p]=Zero();
      if( geom.displacements[p]==0){
	self_stencil=p;
      }

      Integer block=(FineGrid->_rdimensions[dir])/(Grid()->_rdimensions[dir]);

      LatticeCoordinate(coor,dir);

      ///////////////////////////////////////////////////////
      // Work out even and odd block checkerboarding for fast diagonal term
      ///////////////////////////////////////////////////////
      if ( disp==1 ) {
	bcb   = bcb + div(coor,block);
      }
	
      if ( disp==0 ) {
	  masks[p]= Zero();
      } else if ( disp==1 ) {
	masks[p] = where(mod(coor,block)==(block-1),one,zero);
      } else if ( disp==-1 ) {
	masks[p] = where(mod(coor,block)==(Integer)0,one,zero);
      }
    }
    evenmask = where(mod(bcb,2)==(Integer)0,one,zero);
    oddmask  = one-evenmask;

    /*
    {
      phi=Subspace.subspace[0];
      linop.OpDirAll(phi,Mphi_p);
      for(int p=0;p<geom.npoint-1;p++){
	int dir=geom.directions[p];
	int disp=geom.displacements[p];
	linop.OpDir(phi,Mphi,dir,disp);
	Mphi=Mphi-Mphi_p[p];
	std::cout << GridLogMessage <<" Direction mapping check " <<norm2(Mphi)<<std::endl;
      }
    }
*/
    assert(self_stencil!=-1);
    int lhermitian=hermitian;

    for(int i=0;i<nbasis;i++){

      phi=Subspace.subspace[i];

      std::cout << GridLogMessage<< "CoarsenMatrix vector "<<i << std::endl;
      linop.OpDirAll(phi,Mphi_p);
      linop.OpDiag  (phi,Mphi_p[geom.npoint-1]);

      for(int p=0;p<geom.npoint;p++){ 

	Mphi = Mphi_p[p];

	int dir   = geom.directions[p];
	int disp  = geom.displacements[p];

	if ( (disp==-1) || (!lhermitian ) ) {

	  ////////////////////////////////////////////////////////////////////////
	  // Pick out contributions coming from this cell and neighbour cell
	  ////////////////////////////////////////////////////////////////////////
	  omask = masks[p];
	  imask = one-omask;
	
	  for(int j=0;j<nbasis;j++){
	    
	    blockMaskedInnerProduct(oZProj,omask,Subspace.subspace[j],Mphi);
	    
	    autoView( iZProj_v , iZProj, AcceleratorRead) ;
	    autoView( oZProj_v , oZProj, AcceleratorRead) ;
	    autoView( A_p     ,  A[p], AcceleratorWrite);
	    autoView( A_self  , A[self_stencil], AcceleratorWrite);

	    accelerator_for(ss, Grid()->oSites(), Fobj::Nsimd(),{ coalescedWrite(A_p[ss](j,i),oZProj_v(ss)); });
	    if ( lhermitian && (disp==-1) ) {
	      for(int pp=0;pp<geom.npoint;pp++){// Find the opposite link and set <j|A|i> = <i|A|j>*
		int dirp   = geom.directions[pp];
		int dispp  = geom.displacements[pp];
		if ( (dirp==dir) && (dispp==1) ){
		  auto sft = conjugate(Cshift(oZProj,dir,1));
		  autoView( sft_v    ,  sft  , AcceleratorWrite);
		  autoView( A_pp     ,  A[pp], AcceleratorWrite);
		  accelerator_for(ss, Grid()->oSites(), Fobj::Nsimd(),{ coalescedWrite(A_pp[ss](i,j),sft_v(ss)); });
		}
	      }
	    }
	  }
	}
      }

      std::cout << GridLogMessage<< "CoarsenMatrix Diag "<<std::endl;
      ///////////////////////////////////////////
      // Faster alternate self coupling.. use hermiticity to save 2x
      ///////////////////////////////////////////
      {
	mult(tmp,phi,evenmask);  linop.Op(tmp,Mphie);
	mult(tmp,phi,oddmask );  linop.Op(tmp,Mphio);

	{
	  autoView( tmp_      , tmp, AcceleratorWrite);
	  autoView( evenmask_ , evenmask, AcceleratorRead);
	  autoView( oddmask_  ,  oddmask, AcceleratorRead);
	  autoView( Mphie_    ,  Mphie, AcceleratorRead);
	  autoView( Mphio_    ,  Mphio, AcceleratorRead);
	  accelerator_for(ss, FineGrid->oSites(), Fobj::Nsimd(),{ 
	      coalescedWrite(tmp_[ss],evenmask_(ss)*Mphie_(ss) + oddmask_(ss)*Mphio_(ss));
	    });
	}

	blockProject(SelfProj,tmp,Subspace.subspace);

	autoView( SelfProj_ , SelfProj, AcceleratorRead);
	autoView( A_self  , A[self_stencil], AcceleratorWrite);

	accelerator_for(ss, Grid()->oSites(), Fobj::Nsimd(),{
	  for(int j=0;j<nbasis;j++){
	    coalescedWrite(A_self[ss](j,i), SelfProj_(ss)(j));
	  }
	});

      }
    }

    MemoryManager::PrintBytes();

    // Auto self test
    Test( Subspace,FineGrid,linop);

#if 0
    ///////////////////////////
    // test code worth preserving in if block
    ///////////////////////////
    std::cout<<GridLogMessage<< " Computed matrix elements "<< self_stencil <<std::endl;
    for(int p=0;p<geom.npoint;p++){
      std::cout<<GridLogMessage<< "A["<<p<<"]" << std::endl;
      std::cout<<GridLogMessage<< "\n"<<A[p] << std::endl;
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

  }

};

NAMESPACE_END(Grid);
#endif
