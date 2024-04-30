/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/Aggregates.h

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
#pragma once

NAMESPACE_BEGIN(Grid);

inline RealD AggregatePowerLaw(RealD x)
{
  //  return std::pow(x,-4);
  //  return std::pow(x,-3);
  return std::pow(x,-5);
}

template<class Fobj,class CComplex,int nbasis>
class Aggregation {
public:
  constexpr int Nbasis(void) { return nbasis; };
  
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
    //    std::cout << GridLogMessage <<" Block Gramm-Schmidt pass 1"<<std::endl;
    blockOrthogonalise(InnerProd,subspace);
  } 
  void ProjectToSubspace(CoarseVector &CoarseVec,const FineField &FineVec){
    blockProject(CoarseVec,FineVec,subspace);
  }
  void PromoteFromSubspace(const CoarseVector &CoarseVec,FineField &FineVec){
    FineVec.Checkerboard() = subspace[0].Checkerboard();
    blockPromote(CoarseVec,FineVec,subspace);
  }

  virtual void CreateSubspaceRandom(GridParallelRNG  &RNG) {
    int nn=nbasis;
    RealD scale;
    FineField noise(FineGrid);
    for(int b=0;b<nn;b++){
      subspace[b] = Zero();
      gaussian(RNG,noise);
      scale = std::pow(norm2(noise),-0.5); 
      noise=noise*scale;
      subspace[b] = noise;
    }
  }
  virtual void CreateSubspace(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,int nn=nbasis)
  {

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

    std::cout << GridLogMessage<<" Chebyshev subspace pass-1 : ord "<<orderfilter<<" ["<<lo<<","<<hi<<"]"<<std::endl;
    std::cout << GridLogMessage<<" Chebyshev subspace pass-2 : nbasis"<<nn<<" min "
	      <<ordermin<<" step "<<orderstep
	      <<" lo"<<filterlo<<std::endl;

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
	accelerator_for(ss, FineGrid->oSites(), Nsimd, {
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
  virtual void CreateSubspaceChebyshev(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,
				       int nn,
				       double hi,
				       double lo,
				       int orderfilter
				       ) {

    RealD scale;

    FineField noise(FineGrid);
    FineField Mn(FineGrid);
    FineField tmp(FineGrid);

    // New normalised noise
    std::cout << GridLogMessage<<" Chebyshev subspace pure noise : ord "<<orderfilter<<" ["<<lo<<","<<hi<<"]"<<std::endl;
    std::cout << GridLogMessage<<" Chebyshev subspace pure noise  : nbasis "<<nn<<std::endl;


    for(int b =0;b<nbasis;b++)
    {
      gaussian(RNG,noise);
      scale = std::pow(norm2(noise),-0.5); 
      noise=noise*scale;

      // Initial matrix element
      hermop.Op(noise,Mn);
      if(b==0) std::cout<<GridLogMessage << "noise <n|MdagM|n> "<<norm2(Mn)<<std::endl;

      // Filter
      Chebyshev<FineField> Cheb(lo,hi,orderfilter);
      Cheb(hermop,noise,Mn);
      scale = std::pow(norm2(Mn),-0.5); 	Mn=Mn*scale;

      // Refine
      Chebyshev<FineField> PowerLaw(lo,hi,1000,AggregatePowerLaw);
      noise = Mn;
      PowerLaw(hermop,noise,Mn);
      scale = std::pow(norm2(Mn),-0.5); 	Mn=Mn*scale;

      // normalise
      subspace[b]   = Mn;
      hermop.Op(Mn,tmp); 
      std::cout<<GridLogMessage << "filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;
    }

  }

  virtual void CreateSubspaceChebyshevPowerLaw(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,
					       int nn,
					       double hi,
					       int orderfilter
					       ) {

    RealD scale;

    FineField noise(FineGrid);
    FineField Mn(FineGrid);
    FineField tmp(FineGrid);

    // New normalised noise
    std::cout << GridLogMessage<<" Chebyshev subspace pure noise : ord "<<orderfilter<<" [0,"<<hi<<"]"<<std::endl;
    std::cout << GridLogMessage<<" Chebyshev subspace pure noise  : nbasis "<<nn<<std::endl;

    for(int b =0;b<nbasis;b++)
    {
      gaussian(RNG,noise);
      scale = std::pow(norm2(noise),-0.5); 
      noise=noise*scale;

      // Initial matrix element
      hermop.Op(noise,Mn);
      if(b==0) std::cout<<GridLogMessage << "noise <n|MdagM|n> "<<norm2(Mn)<<std::endl;
      // Filter
      Chebyshev<FineField> Cheb(0.0,hi,orderfilter,AggregatePowerLaw);
      Cheb(hermop,noise,Mn);
      // normalise
      scale = std::pow(norm2(Mn),-0.5); 	Mn=Mn*scale;
      subspace[b]   = Mn;
      hermop.Op(Mn,tmp); 
      std::cout<<GridLogMessage << "filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;
    }

  }
  virtual void CreateSubspaceChebyshevNew(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,
					  double hi
					  ) {

    RealD scale;

    FineField noise(FineGrid);
    FineField Mn(FineGrid);
    FineField tmp(FineGrid);

    // New normalised noise
    for(int b =0;b<nbasis;b++)
    {
      gaussian(RNG,noise);
      scale = std::pow(norm2(noise),-0.5); 
      noise=noise*scale;

      // Initial matrix element
      hermop.Op(noise,Mn);
      if(b==0) std::cout<<GridLogMessage << "noise <n|MdagM|n> "<<norm2(Mn)<<std::endl;
      // Filter
      //#opt2(x) =  acheb(x,3,90,300)* acheb(x,1,90,50) * acheb(x,0.5,90,200) * acheb(x,0.05,90,400) * acheb(x,0.01,90,1500)
      /*266
      Chebyshev<FineField> Cheb1(3.0,hi,300);
      Chebyshev<FineField> Cheb2(1.0,hi,50);
      Chebyshev<FineField> Cheb3(0.5,hi,300);
      Chebyshev<FineField> Cheb4(0.05,hi,500);
      Chebyshev<FineField> Cheb5(0.01,hi,2000);
      */
      /* 242 */
      /*
      Chebyshev<FineField> Cheb3(0.1,hi,300);
      Chebyshev<FineField> Cheb2(0.02,hi,1000);
      Chebyshev<FineField> Cheb1(0.003,hi,2000);
      8?
      */
      /* How many??
      */
      Chebyshev<FineField> Cheb2(0.001,hi,2500); // 169 iters on HDCG after refine
      Chebyshev<FineField> Cheb1(0.02,hi,600);

      //      Chebyshev<FineField> Cheb2(0.001,hi,1500);
      //      Chebyshev<FineField> Cheb1(0.02,hi,600);
      Cheb1(hermop,noise,Mn); scale = std::pow(norm2(Mn),-0.5); 	noise=Mn*scale;
      hermop.Op(noise,tmp); std::cout<<GridLogMessage << "Cheb1 <n|MdagM|n> "<<norm2(tmp)<<std::endl;
      Cheb2(hermop,noise,Mn); scale = std::pow(norm2(Mn),-0.5); 	noise=Mn*scale;
      hermop.Op(noise,tmp); std::cout<<GridLogMessage << "Cheb2 <n|MdagM|n> "<<norm2(tmp)<<std::endl;
      //      Cheb3(hermop,noise,Mn); scale = std::pow(norm2(Mn),-0.5); 	noise=Mn*scale;
      //      hermop.Op(noise,tmp); std::cout<<GridLogMessage << "Cheb3 <n|MdagM|n> "<<norm2(tmp)<<std::endl;
      //      Cheb4(hermop,noise,Mn); scale = std::pow(norm2(Mn),-0.5); 	noise=Mn*scale;
      //      hermop.Op(noise,tmp); std::cout<<GridLogMessage << "Cheb4 <n|MdagM|n> "<<norm2(tmp)<<std::endl;
      //      Cheb5(hermop,noise,Mn); scale = std::pow(norm2(Mn),-0.5); 	noise=Mn*scale;
      //      hermop.Op(noise,tmp); std::cout<<GridLogMessage << "Cheb5 <n|MdagM|n> "<<norm2(tmp)<<std::endl;
      subspace[b]   = noise;
      hermop.Op(subspace[b],tmp); 
      std::cout<<GridLogMessage << "filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<< " norm " << norm2(noise)<<std::endl;
    }

  }

  virtual void CreateSubspaceMultishift(GridParallelRNG  &RNG,LinearOperatorBase<FineField> &hermop,
					double Lo,double tol,int maxit)
  {

    RealD scale;

    FineField noise(FineGrid);
    FineField Mn(FineGrid);
    FineField tmp(FineGrid);

    // New normalised noise
    std::cout << GridLogMessage<<" Multishift subspace : Lo "<<Lo<<std::endl;

    // Filter
    // [ 1/6(x+Lo)  - 1/2(x+2Lo) + 1/2(x+3Lo)  -1/6(x+4Lo) = Lo^3 /[ (x+1Lo)(x+2Lo)(x+3Lo)(x+4Lo) ]
    //
    // 1/(x+Lo)  - 1/(x+2 Lo)
    double epsilon      = Lo/3;
    std::vector<RealD> alpha({1.0/6.0,-1.0/2.0,1.0/2.0,-1.0/6.0});
    std::vector<RealD> shifts({Lo,Lo+epsilon,Lo+2*epsilon,Lo+3*epsilon});
    std::vector<RealD> tols({tol,tol,tol,tol});
    std::cout << "sizes "<<alpha.size()<<" "<<shifts.size()<<" "<<tols.size()<<std::endl;

    MultiShiftFunction msf(4,0.0,95.0);
    std::cout << "msf constructed "<<std::endl;
    msf.poles=shifts;
    msf.residues=alpha;
    msf.tolerances=tols;
    msf.norm=0.0;
    msf.order=alpha.size();
    ConjugateGradientMultiShift<FineField> MSCG(maxit,msf);
    
    for(int b =0;b<nbasis;b++)
    {
      gaussian(RNG,noise);
      scale = std::pow(norm2(noise),-0.5); 
      noise=noise*scale;

      // Initial matrix element
      hermop.Op(noise,Mn);
      if(b==0) std::cout<<GridLogMessage << "noise <n|MdagM|n> "<<norm2(Mn)<<std::endl;

      MSCG(hermop,noise,Mn);
      scale = std::pow(norm2(Mn),-0.5); 	Mn=Mn*scale;
      subspace[b]   = Mn;
      hermop.Op(Mn,tmp); 
      std::cout<<GridLogMessage << "filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;

    }

  }
  virtual void RefineSubspace(LinearOperatorBase<FineField> &hermop,
			      double Lo,double tol,int maxit)
  {
    FineField tmp(FineGrid);
    for(int b =0;b<nbasis;b++)
    {
      ConjugateGradient<FineField>  CGsloppy(tol,maxit,false);
      ShiftedHermOpLinearOperator<FineField> ShiftedFineHermOp(hermop,Lo);
      tmp=Zero();
      CGsloppy(hermop,subspace[b],tmp);
      RealD scale = std::pow(norm2(tmp),-0.5); 	tmp=tmp*scale;
      subspace[b]=tmp;
      hermop.Op(subspace[b],tmp);
      std::cout<<GridLogMessage << "filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;
    }
  }
  virtual void RefineSubspaceHDCG(LinearOperatorBase<FineField> &hermop,
				  TwoLevelADEF2mrhs<FineField,CoarseVector> & theHDCG,
				  int nrhs)
  {
    std::vector<FineField> src_mrhs(nrhs,FineGrid);
    std::vector<FineField> res_mrhs(nrhs,FineGrid);
    FineField tmp(FineGrid);
    for(int b =0;b<nbasis;b+=nrhs)
    {
      tmp = subspace[b];
      RealD scale = std::pow(norm2(tmp),-0.5); 	tmp=tmp*scale;
      subspace[b] =tmp;
      hermop.Op(subspace[b],tmp);
      std::cout<<GridLogMessage << "before filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;

      for(int r=0;r<MIN(nbasis-b,nrhs);r++){
	src_mrhs[r] = subspace[b+r];
      }
      for(int r=0;r<nrhs;r++){
	res_mrhs[r] = Zero();
      }
      theHDCG(src_mrhs,res_mrhs);

      for(int r=0;r<MIN(nbasis-b,nrhs);r++){
	tmp = res_mrhs[r];
	RealD scale = std::pow(norm2(tmp),-0.5); tmp=tmp*scale;
	subspace[b+r]=tmp;
      }
      hermop.Op(subspace[b],tmp);
      std::cout<<GridLogMessage << "after filt ["<<b<<"] <n|MdagM|n> "<<norm2(tmp)<<std::endl;
    }
  }

  
  
};
NAMESPACE_END(Grid);

