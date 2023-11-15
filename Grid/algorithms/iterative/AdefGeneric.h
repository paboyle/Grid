    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/AdefGeneric.h

    Copyright (C) 2015

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
#ifndef GRID_ALGORITHMS_ITERATIVE_GENERIC_PCG
#define GRID_ALGORITHMS_ITERATIVE_GENERIC_PCG

  /*
   * Compared to Tang-2009:  P=Pleft. P^T = PRight Q=MssInv. 
   * Script A = SolverMatrix 
   * Script P = Preconditioner
   *
   * Implement ADEF-2
   *
   * Vstart = P^Tx + Qb
   * M1 = P^TM + Q
   * M2=M3=1
   */
NAMESPACE_BEGIN(Grid);

template<class Field>
class TwoLevelCG : public LinearFunction<Field>
{
 public:
  RealD   Tolerance;
  Integer MaxIterations;
  GridBase *grid;

  // Fine operator, Smoother, CoarseSolver
  LinearOperatorBase<Field>   &_FineLinop;
  LinearFunction<Field>   &_Smoother;
  
  // more most opertor functions
  TwoLevelCG(RealD tol,
	     Integer maxit,
	     LinearOperatorBase<Field>   &FineLinop,
	     LinearFunction<Field>       &Smoother,
	     GridBase *fine) : 
      Tolerance(tol), 
      MaxIterations(maxit),
      _FineLinop(FineLinop),
      _Smoother(Smoother)
  {
    grid       = fine;
  };
  
  virtual void operator() (const Field &src, Field &x)
  {
#if 0
    Field resid(grid);
    RealD f;
    RealD rtzp,rtz,a,d,b;
    RealD rptzp;
    
    Field p(grid);
    Field z(grid);
    Field tmp(grid);
    Field mmp(grid);
    Field r  (grid);
    Field mu (grid);
    Field rp (grid);
    
    //Initial residual computation & set up
    double tn;

    GridStopWatch HDCGTimer;
    HDCGTimer.Start();
    //////////////////////////
    // x0 = Vstart -- possibly modify guess
    //////////////////////////
    x=Zero();
    Vstart(x,src);

    // r0 = b -A x0
    _FineLinop.HermOp(x,mmp);

    axpy(r, -1.0, mmp, src);    // Recomputes r=src-x0
    rp=r;

    //////////////////////////////////
    // Compute z = M1 x
    //////////////////////////////////
    PcgM1(r,z);
    rtzp =real(innerProduct(r,z));

    ///////////////////////////////////////
    // Except Def2, M2 is trivial
    ///////////////////////////////////////
    p=z;

    RealD ssq =  norm2(src);
    RealD rsq =  ssq*Tolerance*Tolerance;

    GRID_TRACE("MultiGrid TwoLevel ");
    std::cout<<GridLogMessage<<"HDCG: k=0 residual "<<rtzp<<" target rsq "<<rsq<<" ssq "<<ssq<<std::endl;
    
    for (int k=1;k<=MaxIterations;k++){

      rtz=rtzp;
      d= PcgM3(p,mmp);
      a = rtz/d;

      axpy(x,a,p,x);
      RealD rn = axpy_norm(r,-a,mmp,r);

      PcgM1(r,z);

      rtzp =real(innerProduct(r,z));

      int ipcg=1; // almost free inexact preconditioned CG
      if (ipcg) {
	rptzp =real(innerProduct(rp,z));
      } else {
	rptzp =0;
      }
      b = (rtzp-rptzp)/rtz;

      PcgM2(z,mu); // ADEF-2 this is identity. Axpy possible to eliminate

      axpy(p,b,p,mu);  // mu = A r

      RealD rrn=sqrt(rn/ssq);
      RealD rtn=sqrt(rtz/ssq);
      std::cout<<GridLogMessage<<"HDCG: Pcg k= "<<k<<" residual = "<<rrn<<std::endl;

      if ( ipcg ) {
	axpy(rp,0.0,r,r);
      }

      // Stopping condition
      if ( rn <= rsq ) { 

	HDCGTimer.Stop();
	std::cout<<GridLogMessage<<"HDCG: Pcg converged in "<<k<<" iterations and "<<HDCGTimer.Elapsed()<<std::endl;;

	_FineLinop.HermOp(x,mmp);			  
	axpy(tmp,-1.0,src,mmp);

	RealD  mmpnorm = sqrt(norm2(mmp));
	RealD  xnorm   = sqrt(norm2(x));
	RealD  srcnorm = sqrt(norm2(src));
	RealD  tmpnorm = sqrt(norm2(tmp));
	RealD  true_residual = tmpnorm/srcnorm;
	std::cout<<GridLogMessage
		 <<"HDCG: true residual is "<<true_residual
		 <<" solution "<<xnorm
		 <<" source "<<srcnorm
		 <<" mmp "<<mmpnorm	  
		 <<std::endl;

	return;
      }

    }
    std::cout<<GridLogMessage<<"HDCG: not converged"<<std::endl;
    RealD  xnorm   = sqrt(norm2(x));
    RealD  srcnorm = sqrt(norm2(src));
    std::cout<<GridLogMessage<<"HDCG: non-converged solution "<<xnorm<<" source "<<srcnorm<<std::endl;
    
    return ;
#else
  RealD f;
  RealD rtzp,rtz,a,d,b;
  RealD rptzp;

  /////////////////////////////
  // Set up history vectors
  /////////////////////////////
  int mmax = 20;
  std::vector<Field> p(mmax,grid);
  std::vector<Field> mmp(mmax,grid);
  std::vector<RealD> pAp(mmax);
  Field z(grid);
  Field tmp(grid);
  Field  mp (grid);
  Field  r  (grid);
  Field  mu (grid);

  //Initial residual computation & set up
  RealD guess   = norm2(x);
  RealD src_nrm = norm2(src);

  if ( src_nrm == 0.0 ) {
    std::cout << GridLogMessage<<"HDCG: fPcg given trivial source norm "<<src_nrm<<std::endl;
    x=Zero();
  }
  RealD tn;

  GridStopWatch HDCGTimer;
  HDCGTimer.Start();
  //////////////////////////
  // x0 = Vstart -- possibly modify guess
  //////////////////////////
  Vstart(x,src);

  // r0 = b -A x0
  _FineLinop.HermOp(x,mmp[0]);
  axpy (r, -1.0,mmp[0], src);    // Recomputes r=src-Ax0
  {
    double n1 = norm2(x);
    double n2 = norm2(mmp[0]);
    double n3 = norm2(r);
    std::cout<<GridLogMessage<<"x,vstart,r = "<<n1<<" "<<n2<<" "<<n3<<std::endl;
  }

  //////////////////////////////////
  // Compute z = M1 x
  //////////////////////////////////
  PcgM1(r,z);
  rtzp =real(innerProduct(r,z));

  ///////////////////////////////////////
  // Solve for Mss mu = P A z and set p = z-mu
  // Def2: p = 1 - Q Az = Pright z 
  // Other algos M2 is trivial
  ///////////////////////////////////////
  PcgM2(z,p[0]);

  RealD ssq =  norm2(src);
  RealD rsq =  ssq*Tolerance*Tolerance;

  std::cout << GridLogMessage<<"HDCG: k=0 residual "<<rtzp<<" rsq "<<rsq<<"\n";

  Field pp(grid);

  for (int k=0;k<=MaxIterations;k++){
    
    int peri_k  = k % mmax;
    int peri_kp = (k+1) % mmax;

    rtz=rtzp;
    d= PcgM3(p[peri_k],mmp[peri_k]);
    a = rtz/d;
    
    // Memorise this
    pAp[peri_k] = d;

    
    axpy(x,a,p[peri_k],x);
    RealD rn = axpy_norm(r,-a,mmp[peri_k],r);

    // Compute z = M x
    PcgM1(r,z);

    {
      RealD n1,n2;
      n1=norm2(r);
      n2=norm2(z);
      std::cout << GridLogMessage<<"HDCG::fPcg iteration "<<k<<" : vector r,z "<<n1<<" "<<n2<<"\n";
    }
    rtzp =real(innerProduct(r,z));
    std::cout << GridLogMessage<<"HDCG::fPcg iteration "<<k<<" : inner rtzp "<<rtzp<<"\n";

    //    PcgM2(z,p[0]);
    PcgM2(z,mu); // ADEF-2 this is identity. Axpy possible to eliminate

    p[peri_kp]=mu;

    // Standard search direction  p -> z + b p    ; b = 
    b = (rtzp)/rtz;

    int northog;
    // k=zero  <=> peri_kp=1;        northog = 1
    // k=1     <=> peri_kp=2;        northog = 2
    // ...               ...                  ...
    // k=mmax-2<=> peri_kp=mmax-1;   northog = mmax-1
    // k=mmax-1<=> peri_kp=0;        northog = 1

    //    northog     = (peri_kp==0)?1:peri_kp; // This is the fCG(mmax) algorithm
    northog     = (k>mmax-1)?(mmax-1):k;        // This is the fCG-Tr(mmax-1) algorithm
    
    std::cout<<GridLogMessage<<"HDCG::fPcg iteration "<<k<<" : orthogonalising to last "<<northog<<" vectors\n";
    for(int back=0; back < northog; back++){
      int peri_back = (k-back)%mmax;
      RealD pbApk= real(innerProduct(mmp[peri_back],p[peri_kp]));
      RealD beta = -pbApk/pAp[peri_back];
      axpy(p[peri_kp],beta,p[peri_back],p[peri_kp]);
    }

    RealD rrn=sqrt(rn/ssq);
    RealD rtn=sqrt(rtz/ssq);
    RealD rtnp=sqrt(rtzp/ssq);

    std::cout<<GridLogMessage<<"HDCG: fPcg k= "<<k<<" residual = "<<rrn<<"\n";

    // Stopping condition
    if ( rn <= rsq ) { 

      HDCGTimer.Stop();
      std::cout<<GridLogMessage<<"HDCG: fPcg converged in "<<k<<" iterations and "<<HDCGTimer.Elapsed()<<std::endl;;
      
      _FineLinop.HermOp(x,mmp[0]);			  
      axpy(tmp,-1.0,src,mmp[0]);
      
      RealD  mmpnorm = sqrt(norm2(mmp[0]));
      RealD  xnorm   = sqrt(norm2(x));
      RealD  srcnorm = sqrt(norm2(src));
      RealD  tmpnorm = sqrt(norm2(tmp));
      RealD  true_residual = tmpnorm/srcnorm;
      std::cout<<GridLogMessage
	       <<"HDCG: true residual is "<<true_residual
	       <<" solution "<<xnorm
	       <<" source "<<srcnorm
	       <<" mmp "<<mmpnorm	  
	       <<std::endl;
      
      return;
    }

  }
  HDCGTimer.Stop();
  std::cout<<GridLogMessage<<"HDCG: not converged "<<HDCGTimer.Elapsed()<<std::endl;
  RealD  xnorm   = sqrt(norm2(x));
  RealD  srcnorm = sqrt(norm2(src));
  std::cout<<GridLogMessage<<"HDCG: non-converged solution "<<xnorm<<" source "<<srcnorm<<std::endl;
#endif
  }
  

 public:

  virtual void PcgM1(Field & in, Field & out)     =0;
  virtual void Vstart(Field & x,const Field & src)=0;

  virtual void PcgM2(const Field & in, Field & out) {
    out=in;
  }

  virtual RealD PcgM3(const Field & p, Field & mmp){
    RealD dd;
    _FineLinop.HermOp(p,mmp);
    ComplexD dot = innerProduct(p,mmp);
    dd=real(dot);
    return dd;
  }

  /////////////////////////////////////////////////////////////////////
  // Only Def1 has non-trivial Vout.
  /////////////////////////////////////////////////////////////////////

};
  
template<class Field, class CoarseField, class Aggregation>
class TwoLevelADEF2 : public TwoLevelCG<Field>
{
 public:
  ///////////////////////////////////////////////////////////////////////////////////
  // Need something that knows how to get from Coarse to fine and back again
  //  void ProjectToSubspace(CoarseVector &CoarseVec,const FineField &FineVec){
  //  void PromoteFromSubspace(const CoarseVector &CoarseVec,FineField &FineVec){
  ///////////////////////////////////////////////////////////////////////////////////
  GridBase *coarsegrid;
  Aggregation &_Aggregates;                    
  LinearFunction<CoarseField> &_CoarseSolver;
  LinearFunction<CoarseField> &_CoarseSolverPrecise;
  ///////////////////////////////////////////////////////////////////////////////////
  
  // more most opertor functions
  TwoLevelADEF2(RealD tol,
		Integer maxit,
		LinearOperatorBase<Field>    &FineLinop,
		LinearFunction<Field>        &Smoother,
		LinearFunction<CoarseField>  &CoarseSolver,
		LinearFunction<CoarseField>  &CoarseSolverPrecise,
		Aggregation &Aggregates
		) :
      TwoLevelCG<Field>(tol,maxit,FineLinop,Smoother,Aggregates.FineGrid),
      _CoarseSolver(CoarseSolver),
      _CoarseSolverPrecise(CoarseSolverPrecise),
      _Aggregates(Aggregates)
  {
    coarsegrid = Aggregates.CoarseGrid;
  };

  virtual void PcgM1(Field & in, Field & out)
  {
    GRID_TRACE("MultiGridPreconditioner ");
    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]

    Field tmp(this->grid);
    Field Min(this->grid);
    CoarseField PleftProj(this->coarsegrid);
    CoarseField PleftMss_proj(this->coarsegrid);

    GridStopWatch SmootherTimer;
    GridStopWatch MatrixTimer;
    SmootherTimer.Start();
    this->_Smoother(in,Min);
    SmootherTimer.Stop();

    MatrixTimer.Start();
    this->_FineLinop.HermOp(Min,out);
    MatrixTimer.Stop();
    axpy(tmp,-1.0,out,in);          // tmp  = in - A Min

    GridStopWatch ProjTimer;
    GridStopWatch CoarseTimer;
    GridStopWatch PromTimer;
    ProjTimer.Start();
    this->_Aggregates.ProjectToSubspace(PleftProj,tmp);     
    ProjTimer.Stop();
    CoarseTimer.Start();
    this->_CoarseSolver(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    CoarseTimer.Stop();
    PromTimer.Start();
    this->_Aggregates.PromoteFromSubspace(PleftMss_proj,tmp);// tmp = Q[in - A Min]  
    PromTimer.Stop();
    std::cout << GridLogPerformance << "PcgM1 breakdown "<<std::endl;
    std::cout << GridLogPerformance << "\tSmoother   " << SmootherTimer.Elapsed() <<std::endl;
    std::cout << GridLogPerformance << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
    std::cout << GridLogPerformance << "\tProj       " << ProjTimer.Elapsed() <<std::endl;
    std::cout << GridLogPerformance << "\tCoarse     " << CoarseTimer.Elapsed() <<std::endl;
    std::cout << GridLogPerformance << "\tProm       " << PromTimer.Elapsed() <<std::endl;

    axpy(out,1.0,Min,tmp); // Min+tmp
  }

  virtual void Vstart(Field & x,const Field & src)
  {
    ///////////////////////////////////
    // Choose x_0 such that 
    // x_0 = guess +  (A_ss^inv) r_s = guess + Ass_inv [src -Aguess]
    //                               = [1 - Ass_inv A] Guess + Assinv src
    //                               = P^T guess + Assinv src 
    //                               = Vstart  [Tang notation]
    // This gives:
    // W^T (src - A x_0) = src_s - A guess_s - r_s
    //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
    //                   = 0 
    ///////////////////////////////////
    Field r(this->grid);
    Field mmp(this->grid);
    CoarseField PleftProj(this->coarsegrid);
    CoarseField PleftMss_proj(this->coarsegrid);

    this->_Aggregates.ProjectToSubspace(PleftProj,src);     
    this->_CoarseSolverPrecise(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    this->_Aggregates.PromoteFromSubspace(PleftMss_proj,x);  

  }

};

template<class Field>
class TwoLevelADEF1defl : public TwoLevelCG<Field>
{
public:
  const std::vector<Field> &evec;
  const std::vector<RealD> &eval;
  
  TwoLevelADEF1defl(RealD tol,
		   Integer maxit,
		   LinearOperatorBase<Field>   &FineLinop,
		   LinearFunction<Field>   &Smoother,
		   std::vector<Field> &_evec,
		   std::vector<RealD> &_eval) : 
    TwoLevelCG<Field>(tol,maxit,FineLinop,Smoother,_evec[0].Grid()),
    evec(_evec),
    eval(_eval)
  {};

  // Can just inherit existing M2
  // Can just inherit existing M3

  // Simple vstart - do nothing
  virtual void Vstart(Field & x,const Field & src){
    x=src; // Could apply Q
  };

  // Override PcgM1
  virtual void PcgM1(Field & in, Field & out)
  {
    GRID_TRACE("EvecPreconditioner ");
    int N=evec.size();
    Field Pin(this->grid);
    Field Qin(this->grid);

    //MP  + Q = M(1-AQ) + Q = M
    // // If we are eigenvector deflating in coarse space
    // // Q   = Sum_i |phi_i> 1/lambda_i <phi_i|
    // // A Q = Sum_i |phi_i> <phi_i|
    // // M(1-AQ) = M(1-proj) + Q
    Qin.Checkerboard()=in.Checkerboard();
    Qin = Zero();
    Pin = in;
    for (int i=0;i<N;i++) {
      const Field& tmp = evec[i];
      auto ip = TensorRemove(innerProduct(tmp,in));
      axpy(Qin, ip / eval[i],tmp,Qin);
      axpy(Pin, -ip ,tmp,Pin);
    }

    this->_Smoother(Pin,out);

    out = out + Qin;
  }
};

NAMESPACE_END(Grid);

#endif
