#ifndef GRID_PGCR_H
#define GRID_PGCR_H

namespace Grid {

    /////////////////////////////////////////////////////////////
    // Base classes for iterative processes based on operators
    // single input vec, single output vec.
    /////////////////////////////////////////////////////////////

  template<class Field> 
    class ConjugateResidual : public OperatorFunction<Field> {
  public:                                                
    RealD   Tolerance;
    Integer MaxIterations;
    int verbose;

    ConjugateResidual(RealD tol,Integer maxit) : Tolerance(tol), MaxIterations(maxit) { 
      verbose=0;
    };

    void operator() (LinearOperatorBase<Field> &Linop,const Field &src, Field &psi){

      RealD a, b, c, d;
      RealD cp, ssq,rsq;
      
      RealD rAr, rAAr, rArp;
      RealD pAp, pAAp;

      GridBase *grid = src._grid;
      psi=zero;
      Field r(grid),  p(grid), Ap(grid), Ar(grid);
      
      r=src;
      p=src;

      Linop.HermOpAndNorm(p,Ap,pAp,pAAp);
      Linop.HermOpAndNorm(r,Ar,rAr,rAAr);

      if(verbose) std::cout << "pAp, pAAp"<< pAp<<" "<<pAAp<<std::endl;
      if(verbose) std::cout << "rAr, rAAr"<< rAr<<" "<<rAAr<<std::endl;

      cp =norm2(r);
      ssq=norm2(src);
      rsq=Tolerance*Tolerance*ssq;

      std::cout<<"ConjugateResidual: iteration " <<0<<" residual "<<cp<< " target"<< rsq<<std::endl;

      for(int k=1;k<MaxIterations;k++){

	a = rAr/pAAp;

	axpy(psi,a,p,psi);

	cp = axpy_norm(r,-a,Ap,r);

	rArp=rAr;

	Linop.HermOpAndNorm(r,Ar,rAr,rAAr);

	b   =rAr/rArp;
 
	axpy(p,b,p,r);
	pAAp=axpy_norm(Ap,b,Ap,Ar);

	if(verbose) std::cout<<"ConjugateResidual: iteration " <<k<<" residual "<<cp<< " target"<< rsq<<std::endl;

	if(cp<rsq) {
	  Linop.HermOp(psi,Ap);
	  axpy(r,-1.0,src,Ap);
	  RealD true_resid = norm2(r);
	  std::cout<<"ConjugateResidual: Converged on iteration " <<k
		   << " computed residual "<<sqrt(cp/ssq)
	           << " true residual "<<true_resid
	           << " target "       <<Tolerance <<std::endl;
	  return;
	}

      }

      std::cout<<"ConjugateResidual did NOT converge"<<std::endl;
      assert(0);
    }
  };
}
#endif
