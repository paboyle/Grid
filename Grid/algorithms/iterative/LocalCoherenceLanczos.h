/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/LocalCoherenceLanczos.h

    Copyright (C) 2015

Author: Christoph Lehner <clehner@bnl.gov>
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
#ifndef GRID_LOCAL_COHERENCE_IRL_H
#define GRID_LOCAL_COHERENCE_IRL_H

NAMESPACE_BEGIN(Grid); 

struct LanczosParams : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParams,
				  ChebyParams, Cheby,/*Chebyshev*/
				  int, Nstop,    /*Vecs in Lanczos must converge Nstop < Nk < Nm*/
				  int, Nk,       /*Vecs in Lanczos seek converge*/
				  int, Nm,       /*Total vecs in Lanczos include restart*/
				  RealD, resid,  /*residual*/
 				  int, MaxIt, 
				  RealD, betastp,  /* ? */
				  int, MinRes);    // Must restart
};

//This class is the input parameter class for some testing programs
struct LocalCoherenceLanczosParams : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LocalCoherenceLanczosParams,
				  bool, saveEvecs,
				  bool, doFine,
				  bool, doFineRead,
				  bool, doCoarse,
	       			  bool, doCoarseRead,
				  LanczosParams, FineParams,
				  LanczosParams, CoarseParams,
				  ChebyParams,   Smoother,
				  RealD        , coarse_relax_tol,
				  std::vector<int>, blockSize,
				  std::string, config,
				  std::vector < ComplexD  >, omega,
				  RealD, mass,
				  RealD, M5);
};

// Duplicate functionality; ProjectedFunctionHermOp could be used with the trivial function
template<class Fobj,class CComplex,int nbasis>
class ProjectedHermOp : public LinearFunction<Lattice<iVector<CComplex,nbasis > > > {
public:
  using LinearFunction<Lattice<iVector<CComplex,nbasis > > >::operator();
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<CComplex>   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj>          FineField;

  LinearOperatorBase<FineField> &_Linop;
  std::vector<FineField>        &subspace;

  ProjectedHermOp(LinearOperatorBase<FineField>& linop, std::vector<FineField> & _subspace) : 
    _Linop(linop), subspace(_subspace)
  {  
    assert(subspace.size() >0);
  };

  void operator()(const CoarseField& in, CoarseField& out) {
    GridBase *FineGrid = subspace[0].Grid();    
    int   checkerboard = subspace[0].Checkerboard();

    FineField fin (FineGrid);     fin.Checkerboard()= checkerboard;
    FineField fout(FineGrid);   fout.Checkerboard() = checkerboard;

    blockPromote(in,fin,subspace);       std::cout<<GridLogIRL<<"ProjectedHermop : Promote to fine"<<std::endl;
    _Linop.HermOp(fin,fout);                   std::cout<<GridLogIRL<<"ProjectedHermop : HermOp (fine) "<<std::endl;
    blockProject(out,fout,subspace);     std::cout<<GridLogIRL<<"ProjectedHermop : Project to coarse "<<std::endl;
  }
};

template<class Fobj,class CComplex,int nbasis>
class ProjectedFunctionHermOp : public LinearFunction<Lattice<iVector<CComplex,nbasis > > > {
public:
  using LinearFunction<Lattice<iVector<CComplex,nbasis > > >::operator();
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<CComplex>   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj>          FineField;


  OperatorFunction<FineField>   & _poly;
  LinearOperatorBase<FineField> &_Linop;
  std::vector<FineField>        &subspace;

  ProjectedFunctionHermOp(OperatorFunction<FineField> & poly,
			  LinearOperatorBase<FineField>& linop, 
			  std::vector<FineField> & _subspace) :
    _poly(poly),
    _Linop(linop),
    subspace(_subspace)
  {  };

  void operator()(const CoarseField& in, CoarseField& out) {

    GridBase *FineGrid = subspace[0].Grid();    
    int   checkerboard = subspace[0].Checkerboard();

    FineField fin (FineGrid); fin.Checkerboard() =checkerboard;
    FineField fout(FineGrid);fout.Checkerboard() =checkerboard;
    
    blockPromote(in,fin,subspace);             std::cout<<GridLogIRL<<"ProjectedFunctionHermop : Promote to fine"<<std::endl;
    _poly(_Linop,fin,fout);                    std::cout<<GridLogIRL<<"ProjectedFunctionHermop : Poly "<<std::endl;
    blockProject(out,fout,subspace);           std::cout<<GridLogIRL<<"ProjectedFunctionHermop : Project to coarse "<<std::endl;
  }
};

template<class Fobj,class CComplex,int nbasis>
class ImplicitlyRestartedLanczosSmoothedTester  : public ImplicitlyRestartedLanczosTester<Lattice<iVector<CComplex,nbasis > > >
{
public:
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<CComplex>   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj>          FineField;

  LinearFunction<CoarseField> & _Poly;
  OperatorFunction<FineField>   & _smoother;
  LinearOperatorBase<FineField> &_Linop;
  RealD                             _coarse_relax_tol;
  std::vector<FineField>        &_subspace;

  int _largestEvalIdxForReport; //The convergence of the LCL is based on the evals of the coarse grid operator, not those of the underlying fine grid operator
                                //As a result we do not know what the eval range of the fine operator is until the very end, making tuning the Cheby bounds very difficult
                                //To work around this issue, every restart we separately reconstruct the fine operator eval for the lowest and highest evec and print these
                                //out alongside the evals of the coarse operator. To do so we need to know the index of the largest eval (i.e. Nstop-1)
                                //NOTE: If largestEvalIdxForReport=-1 (default) then this is not performed
  
  ImplicitlyRestartedLanczosSmoothedTester(LinearFunction<CoarseField>   &Poly,
					   OperatorFunction<FineField>   &smoother,
					   LinearOperatorBase<FineField> &Linop,
					   std::vector<FineField>        &subspace,
					   RealD coarse_relax_tol=5.0e3,
					   int largestEvalIdxForReport=-1) 
    : _smoother(smoother), _Linop(Linop), _Poly(Poly), _subspace(subspace),
      _coarse_relax_tol(coarse_relax_tol), _largestEvalIdxForReport(largestEvalIdxForReport)
  {    };

  //evalMaxApprox: approximation of largest eval of the fine Chebyshev operator (suitably wrapped by block projection)
  int TestConvergence(int j,RealD eresid,CoarseField &B, RealD &eval,RealD evalMaxApprox)
  {
    CoarseField v(B);
    RealD eval_poly = eval;

    // Apply operator
    _Poly(B,v);

    RealD vnum = real(innerProduct(B,v)); // HermOp.
    RealD vden = norm2(B);
    RealD vv0  = norm2(v);
    eval   = vnum/vden;
    v -= eval*B;

    RealD vv = norm2(v) / ::pow(evalMaxApprox,2.0);

    std::cout.precision(13);
    std::cout<<GridLogIRL  << "[" << std::setw(3)<<j<<"] "
	     <<"eval = "<<std::setw(25)<< eval << " (" << eval_poly << ")"
	     <<" |H B[i] - eval[i]B[i]|^2 / evalMaxApprox^2 " << std::setw(25) << vv
	     <<std::endl;

    if(_largestEvalIdxForReport != -1 && (j==0 || j==_largestEvalIdxForReport)){
      std::cout<<GridLogIRL << "Estimating true eval of fine grid operator for eval idx " << j << std::endl;
      RealD tmp_eval;
      ReconstructEval(j,eresid,B,tmp_eval,1.0); //don't use evalMaxApprox of coarse operator! (cf below)
    }
    
    int conv=0;
    if( (vv<eresid*eresid) ) conv = 1;
    return conv;
  }

  //This function is called at the end of the coarse grid Lanczos. It promotes the coarse eigenvector 'B' to the fine grid,
  //applies a smoother to the result then computes the computes the *fine grid* eigenvalue (output as 'eval').

  //evalMaxApprox should be the approximation of the largest eval of the fine Hermop. However when this function is called by IRL it actually passes the largest eval of the *Chebyshev* operator (as this is the max approx used for the TestConvergence above)
  //As the largest eval of the Chebyshev is typically several orders of magnitude larger this makes the convergence test pass even when it should not.
  //We therefore ignore evalMaxApprox here and use a value of 1.0 (note this value is already used by TestCoarse)
  int ReconstructEval(int j,RealD eresid,CoarseField &B, RealD &eval,RealD evalMaxApprox)  
  {
    evalMaxApprox = 1.0; //cf above
    GridBase *FineGrid = _subspace[0].Grid();    
    int checkerboard   = _subspace[0].Checkerboard();
    FineField fB(FineGrid);fB.Checkerboard() =checkerboard;
    FineField fv(FineGrid);fv.Checkerboard() =checkerboard;

    blockPromote(B,fv,_subspace);  
    
    _smoother(_Linop,fv,fB); 

    RealD eval_poly = eval;
    _Linop.HermOp(fB,fv);

    RealD vnum = real(innerProduct(fB,fv)); // HermOp.
    RealD vden = norm2(fB);
    RealD vv0  = norm2(fv);
    eval   = vnum/vden;
    fv -= eval*fB;
    RealD vv = norm2(fv) / ::pow(evalMaxApprox,2.0);
    if ( j > nbasis ) eresid = eresid*_coarse_relax_tol;
    
    std::cout.precision(13);
    std::cout<<GridLogIRL  << "[" << std::setw(3)<<j<<"] "
	     <<"eval = "<<std::setw(25)<< eval << " (" << eval_poly << ")"
	     <<" |H B[i] - eval[i]B[i]|^2 / evalMaxApprox^2 " << std::setw(25) << vv << " target " << eresid*eresid
	     <<std::endl;
    if( (vv<eresid*eresid) ) return 1;
    return 0;
  }
};

////////////////////////////////////////////
// Make serializable Lanczos params
////////////////////////////////////////////
template<class Fobj,class CComplex,int nbasis>
class LocalCoherenceLanczos 
{
public:
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CComplex>                   CoarseScalar; // used for inner products on fine field
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<Fobj>                       FineField;

protected:
  GridBase *_CoarseGrid;
  GridBase *_FineGrid;
  int _checkerboard;
  LinearOperatorBase<FineField>                 & _FineOp;
  
  std::vector<RealD>                              &evals_fine;
  std::vector<RealD>                              &evals_coarse; 
  std::vector<FineField>                          &subspace;
  std::vector<CoarseField>                        &evec_coarse;

private:
  std::vector<RealD>                              _evals_fine;
  std::vector<RealD>                              _evals_coarse; 
  std::vector<FineField>                          _subspace;
  std::vector<CoarseField>                        _evec_coarse;

public:

  LocalCoherenceLanczos(GridBase *FineGrid,
			GridBase *CoarseGrid,
			LinearOperatorBase<FineField> &FineOp,
			int checkerboard) :
    _CoarseGrid(CoarseGrid),
    _FineGrid(FineGrid),
    _FineOp(FineOp),
    _checkerboard(checkerboard),
    evals_fine  (_evals_fine),
    evals_coarse(_evals_coarse),
    subspace    (_subspace),
    evec_coarse(_evec_coarse)
  {
    evals_fine.resize(0);
    evals_coarse.resize(0);
  };
  //////////////////////////////////////////////////////////////////////////
  // Alternate constructore, external storage for use by Hadrons module
  //////////////////////////////////////////////////////////////////////////
  LocalCoherenceLanczos(GridBase *FineGrid,
			GridBase *CoarseGrid,
			LinearOperatorBase<FineField> &FineOp,
			int checkerboard,
			std::vector<FineField>   &ext_subspace,
			std::vector<CoarseField> &ext_coarse,
			std::vector<RealD>       &ext_eval_fine,
			std::vector<RealD>       &ext_eval_coarse
			) :
    _CoarseGrid(CoarseGrid),
    _FineGrid(FineGrid),
    _FineOp(FineOp),
    _checkerboard(checkerboard),
    evals_fine  (ext_eval_fine), 
    evals_coarse(ext_eval_coarse),
    subspace    (ext_subspace),
    evec_coarse (ext_coarse)
  {
    evals_fine.resize(0);
    evals_coarse.resize(0);
  };

  //The block inner product is the inner product on the fine grid locally summed over the blocks
  //to give a Lattice<Scalar> on the coarse grid. This function orthnormalizes the fine-grid subspace
  //vectors under the block inner product. This step must be performed after computing the fine grid
  //eigenvectors and before computing the coarse grid eigenvectors.    
  void Orthogonalise(void ) {
    CoarseScalar InnerProd(_CoarseGrid);
    std::cout << GridLogMessage <<" Gramm-Schmidt pass 1"<<std::endl;
    blockOrthogonalise(InnerProd,subspace);
    std::cout << GridLogMessage <<" Gramm-Schmidt pass 2"<<std::endl;
    blockOrthogonalise(InnerProd,subspace);
  };

  template<typename T>  static RealD normalise(T& v) 
  {
    RealD nn = norm2(v);
    nn = ::sqrt(nn);
    v = v * (1.0/nn);
    return nn;
  }
  /*
  void fakeFine(void)
  {
    int Nk = nbasis;
    subspace.resize(Nk,_FineGrid);
    subspace[0]=1.0;
    subspace[0].Checkerboard()=_checkerboard;
    normalise(subspace[0]);
    PlainHermOp<FineField>    Op(_FineOp);
    for(int k=1;k<Nk;k++){
      subspace[k].Checkerboard()=_checkerboard;
      Op(subspace[k-1],subspace[k]);
      normalise(subspace[k]);
    }
  }
  */

  void testFine(RealD resid) 
  {
    assert(evals_fine.size() == nbasis);
    assert(subspace.size() == nbasis);
    PlainHermOp<FineField>    Op(_FineOp);
    ImplicitlyRestartedLanczosHermOpTester<FineField> SimpleTester(Op);
    for(int k=0;k<nbasis;k++){
      assert(SimpleTester.ReconstructEval(k,resid,subspace[k],evals_fine[k],1.0)==1);
    }
  }

  //While this method serves to check the coarse eigenvectors, it also recomputes the eigenvalues from the smoothed reconstructed eigenvectors
  //hence the smoother can be tuned after running the coarse Lanczos by using a different smoother here
  void testCoarse(RealD resid,ChebyParams cheby_smooth,RealD relax) 
  {
    assert(evals_fine.size() == nbasis);
    assert(subspace.size() == nbasis);
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // create a smoother and see if we can get a cheap convergence test and smooth inside the IRL
    //////////////////////////////////////////////////////////////////////////////////////////////////
    Chebyshev<FineField>                          ChebySmooth(cheby_smooth);
    ProjectedFunctionHermOp<Fobj,CComplex,nbasis> ChebyOp (ChebySmooth,_FineOp,subspace);
    ImplicitlyRestartedLanczosSmoothedTester<Fobj,CComplex,nbasis> ChebySmoothTester(ChebyOp,ChebySmooth,_FineOp,subspace,relax);

    for(int k=0;k<evec_coarse.size();k++){
      if ( k < nbasis ) { 
	assert(ChebySmoothTester.ReconstructEval(k,resid,evec_coarse[k],evals_coarse[k],1.0)==1);
      } else { 
	assert(ChebySmoothTester.ReconstructEval(k,resid*relax,evec_coarse[k],evals_coarse[k],1.0)==1);
      }
    }
  }

  void calcFine(ChebyParams cheby_parms,int Nstop,int Nk,int Nm,RealD resid, 
		RealD MaxIt, RealD betastp, int MinRes)
  {
    assert(nbasis<=Nm);
    Chebyshev<FineField>      Cheby(cheby_parms);
    FunctionHermOp<FineField> ChebyOp(Cheby,_FineOp);
    PlainHermOp<FineField>    Op(_FineOp);

    evals_fine.resize(Nm);
    subspace.resize(Nm,_FineGrid);

    ImplicitlyRestartedLanczos<FineField> IRL(ChebyOp,Op,Nstop,Nk,Nm,resid,MaxIt,betastp,MinRes);

    FineField src(_FineGrid); 
    typedef typename FineField::scalar_type Scalar;
    // src=1.0; 
    src=Scalar(1.0); 
    src.Checkerboard() = _checkerboard;

    int Nconv;
    IRL.calc(evals_fine,subspace,src,Nconv,false);
    
    // Shrink down to number saved
    assert(Nstop>=nbasis);
    assert(Nconv>=nbasis);
    evals_fine.resize(nbasis);
    subspace.resize(nbasis,_FineGrid);
  }


  //cheby_op: Parameters of the fine grid Chebyshev polynomial used for the Lanczos acceleration
  //cheby_smooth: Parameters of a separate Chebyshev polynomial used after the Lanczos has completed to smooth out high frequency noise in the reconstructed fine grid eigenvectors prior to computing the eigenvalue
  //relax: Reconstructed eigenvectors (post smoothing) are naturally not as precise as true eigenvectors. This factor acts as a multiplier on the stopping condition when determining whether the results satisfy the user provided stopping condition
  void calcCoarse(ChebyParams cheby_op,ChebyParams cheby_smooth,RealD relax,
		  int Nstop, int Nk, int Nm,RealD resid, 
		  RealD MaxIt, RealD betastp, int MinRes)
  {
    Chebyshev<FineField>                          Cheby(cheby_op); //Chebyshev of fine operator on fine grid
    ProjectedHermOp<Fobj,CComplex,nbasis>         Op(_FineOp,subspace); //Fine operator on coarse grid with intermediate fine grid conversion
    ProjectedFunctionHermOp<Fobj,CComplex,nbasis> ChebyOp (Cheby,_FineOp,subspace); //Chebyshev of fine operator on coarse grid with intermediate fine grid conversion
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // create a smoother and see if we can get a cheap convergence test and smooth inside the IRL
    //////////////////////////////////////////////////////////////////////////////////////////////////

    Chebyshev<FineField>                                           ChebySmooth(cheby_smooth); //lower order Chebyshev of fine operator on fine grid used to smooth regenerated eigenvectors
    ImplicitlyRestartedLanczosSmoothedTester<Fobj,CComplex,nbasis> ChebySmoothTester(ChebyOp,ChebySmooth,_FineOp,subspace,relax,Nstop-1); 

    evals_coarse.resize(Nm);
    evec_coarse.resize(Nm,_CoarseGrid);

    CoarseField src(_CoarseGrid);     src=1.0; 

    //Note the "tester" here is also responsible for generating the fine grid eigenvalues which are output into the "evals_coarse" array
    ImplicitlyRestartedLanczos<CoarseField> IRL(ChebyOp,ChebyOp,ChebySmoothTester,Nstop,Nk,Nm,resid,MaxIt,betastp,MinRes);
    int Nconv=0;
    IRL.calc(evals_coarse,evec_coarse,src,Nconv,false);
    assert(Nconv>=Nstop);
    evals_coarse.resize(Nstop);
    evec_coarse.resize (Nstop,_CoarseGrid);
    for (int i=0;i<Nstop;i++){
      std::cout << i << " Coarse eval = " << evals_coarse[i]  << std::endl;
    }
  }

  //Get the fine eigenvector 'i' by reconstruction
  void getFineEvecEval(FineField &evec, RealD &eval, const int i) const{
    blockPromote(evec_coarse[i],evec,subspace);  
    eval = evals_coarse[i];
  }
    
    
};

NAMESPACE_END(Grid);
#endif
