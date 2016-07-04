    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/PrecGeneralisedConjugateResidual.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_PREC_GCR_H
#define GRID_PREC_GCR_H

///////////////////////////////////////////////////////////////////////////////////////////////////////
//VPGCR Abe and Zhang, 2005.
//INTERNATIONAL JOURNAL OF NUMERICAL ANALYSIS AND MODELING
//Computing and Information Volume 2, Number 2, Pages 147-161
//NB. Likely not original reference since they are focussing on a preconditioner variant.
//    but VPGCR was nicely written up in their paper
///////////////////////////////////////////////////////////////////////////////////////////////////////
namespace Grid {

  template<class Field>
    class PrecGeneralisedConjugateResidual : public OperatorFunction<Field> {
  public:                                                
    RealD   Tolerance;
    Integer MaxIterations;
    int verbose;
    int mmax;
    int nstep;
    int steps;
    GridStopWatch PrecTimer;
    GridStopWatch MatTimer;
    GridStopWatch LinalgTimer;

    LinearFunction<Field> &Preconditioner;

   PrecGeneralisedConjugateResidual(RealD tol,Integer maxit,LinearFunction<Field> &Prec,int _mmax,int _nstep) : 
      Tolerance(tol), 
      MaxIterations(maxit),
      Preconditioner(Prec),
      mmax(_mmax),
      nstep(_nstep)
    { 
      verbose=1;
    };

    void operator() (LinearOperatorBase<Field> &Linop,const Field &src, Field &psi){

      psi=zero;
      RealD cp, ssq,rsq;
      ssq=norm2(src);
      rsq=Tolerance*Tolerance*ssq;
      
      Field r(src._grid);

        PrecTimer.Reset();
         MatTimer.Reset();
      LinalgTimer.Reset();

      GridStopWatch SolverTimer;
      SolverTimer.Start();

      steps=0;
      for(int k=0;k<MaxIterations;k++){

	cp=GCRnStep(Linop,src,psi,rsq);

	std::cout<<GridLogMessage<<"VPGCR("<<mmax<<","<<nstep<<") "<< steps <<" steps cp = "<<cp<<std::endl;

	if(cp<rsq) {

	  SolverTimer.Stop();

	  Linop.HermOp(psi,r);
	  axpy(r,-1.0,src,r);
	  RealD tr = norm2(r);
	  std::cout<<GridLogMessage<<"PrecGeneralisedConjugateResidual: Converged on iteration " <<steps
		   << " computed residual "<<sqrt(cp/ssq)
	           << " true residual "    <<sqrt(tr/ssq)
	           << " target "           <<Tolerance <<std::endl;

	  std::cout<<GridLogMessage<<"VPGCR Time elapsed: Total  "<< SolverTimer.Elapsed() <<std::endl;
	  std::cout<<GridLogMessage<<"VPGCR Time elapsed: Precon "<<   PrecTimer.Elapsed() <<std::endl;
	  std::cout<<GridLogMessage<<"VPGCR Time elapsed: Matrix "<<    MatTimer.Elapsed() <<std::endl;
	  std::cout<<GridLogMessage<<"VPGCR Time elapsed: Linalg "<< LinalgTimer.Elapsed() <<std::endl;
	  return;
	}

      }
      std::cout<<GridLogMessage<<"Variable Preconditioned GCR did not converge"<<std::endl;
      assert(0);
    }

    RealD GCRnStep(LinearOperatorBase<Field> &Linop,const Field &src, Field &psi,RealD rsq){

      RealD cp;
      RealD a, b, c, d;
      RealD zAz, zAAz;
      RealD rAq, rq;

      GridBase *grid = src._grid;

      Field r(grid);
      Field z(grid);
      Field tmp(grid);
      Field ttmp(grid);
      Field Az(grid);

      ////////////////////////////////
      // history for flexible orthog
      ////////////////////////////////
      std::vector<Field> q(mmax,grid);
      std::vector<Field> p(mmax,grid);
      std::vector<RealD> qq(mmax);
      
      //////////////////////////////////
      // initial guess x0 is taken as nonzero.
      // r0=src-A x0 = src
      //////////////////////////////////
      MatTimer.Start();
      Linop.HermOpAndNorm(psi,Az,zAz,zAAz); 
      MatTimer.Stop();
      r=src-Az;
      
      /////////////////////
      // p = Prec(r)
      /////////////////////
      PrecTimer.Start();
      Preconditioner(r,z);
      PrecTimer.Stop();

      MatTimer.Start();
      Linop.HermOp(z,tmp); 
      MatTimer.Stop();

      ttmp=tmp;
      tmp=tmp-r;

      /*
      std::cout<<GridLogMessage<<r<<std::endl;
      std::cout<<GridLogMessage<<z<<std::endl;
      std::cout<<GridLogMessage<<ttmp<<std::endl;
      std::cout<<GridLogMessage<<tmp<<std::endl;
      */

      MatTimer.Start();
      Linop.HermOpAndNorm(z,Az,zAz,zAAz); 
      MatTimer.Stop();

      //p[0],q[0],qq[0] 
      p[0]= z;
      q[0]= Az;
      qq[0]= zAAz;

      cp =norm2(r);

      for(int k=0;k<nstep;k++){

	steps++;

	int kp     = k+1;
	int peri_k = k %mmax;
	int peri_kp= kp%mmax;

	rq= real(innerProduct(r,q[peri_k])); // what if rAr not real?
	a = rq/qq[peri_k];

	axpy(psi,a,p[peri_k],psi);         

	cp = axpy_norm(r,-a,q[peri_k],r);  

	if((k==nstep-1)||(cp<rsq)){
	  return cp;
	}

	std::cout<<GridLogMessage<< " VPGCR_step["<<steps<<"]  resid " <<sqrt(cp/rsq)<<std::endl; 

	PrecTimer.Start();
	Preconditioner(r,z);// solve Az = r
	PrecTimer.Stop();

	MatTimer.Start();
	Linop.HermOpAndNorm(z,Az,zAz,zAAz);
	Linop.HermOp(z,tmp);
	MatTimer.Stop();
        tmp=tmp-r;
	std::cout<<GridLogMessage<< " Preconditioner resid " <<sqrt(norm2(tmp)/norm2(r))<<std::endl; 

	q[peri_kp]=Az;
	p[peri_kp]=z;

	int northog = ((kp)>(mmax-1))?(mmax-1):(kp);  // if more than mmax done, we orthog all mmax history.
	for(int back=0;back<northog;back++){

	  int peri_back=(k-back)%mmax;   	  assert((k-back)>=0);

	  b=-real(innerProduct(q[peri_back],Az))/qq[peri_back];
	  p[peri_kp]=p[peri_kp]+b*p[peri_back];
	  q[peri_kp]=q[peri_kp]+b*q[peri_back];

	}
	qq[peri_kp]=norm2(q[peri_kp]); // could use axpy_norm


      }
      assert(0); // never reached
      return cp;
    }
  };
}
#endif
