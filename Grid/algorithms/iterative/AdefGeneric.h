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
   * Vout = x
   */
NAMESPACE_BEGIN(Grid);

template<class Field, class CoarseField, class Aggregation>
class TwoLevelFlexiblePcg : public LinearFunction<Field>
{
 public:
  RealD   Tolerance;
  Integer MaxIterations;
  const int mmax = 1;
  GridBase *grid;
  GridBase *coarsegrid;

  // Fine operator, Smoother, CoarseSolver
  LinearOperatorBase<Field>   &_FineLinop;
  LinearFunction<Field>   &_Smoother;
  LinearFunction<CoarseField> &_CoarseSolver;
  LinearFunction<CoarseField> &_CoarseSolverPrecise;

  // Need something that knows how to get from Coarse to fine and back again
  //  void ProjectToSubspace(CoarseVector &CoarseVec,const FineField &FineVec){
  //  void PromoteFromSubspace(const CoarseVector &CoarseVec,FineField &FineVec){
  Aggregation &_Aggregates;                    
  
  // more most opertor functions
  TwoLevelFlexiblePcg(RealD tol,
		      Integer maxit,
		      LinearOperatorBase<Field>   &FineLinop,
		      LinearFunction<Field>   &Smoother,
		      LinearFunction<CoarseField>  &CoarseSolver,
		      LinearFunction<CoarseField>  &CoarseSolverPrecise,
		      Aggregation &Aggregates
		     ) : 
      Tolerance(tol), 
      MaxIterations(maxit),
      _FineLinop(FineLinop),
      _Smoother(Smoother),
      _CoarseSolver(CoarseSolver),
      _CoarseSolverPrecise(CoarseSolverPrecise),
      _Aggregates(Aggregates)
  {
    coarsegrid = Aggregates.CoarseGrid;
    grid       = Aggregates.FineGrid;
  };

  void Inflexible(Field &src,Field &psi)
  {
    Field resid(grid);
    RealD f;
    RealD rtzp,rtz,a,d,b;
    RealD rptzp;
    
    Field x(grid); 
    Field p(grid);
    Field z(grid);
    Field tmp(grid);
    Field mmp(grid);
    Field r  (grid);
    Field mu (grid);
    Field rp (grid);

    //Initial residual computation & set up
    RealD guess = norm2(psi);
    double tn;

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

	std::cout<<GridLogMessage<<"HDCG: Pcg converged in "<<k<<" iterations"<<std::endl;;

	_FineLinop.HermOp(x,mmp);			  
	axpy(tmp,-1.0,src,mmp);

	RealD  mmpnorm = sqrt(norm2(mmp));
	RealD  psinorm = sqrt(norm2(x));
	RealD  srcnorm = sqrt(norm2(src));
	RealD  tmpnorm = sqrt(norm2(tmp));
	RealD  true_residual = tmpnorm/srcnorm;
	std::cout<<GridLogMessage<<"HDCG: true residual is "<<true_residual
		 <<" solution "<<psinorm<<" source "<<srcnorm<<std::endl;

	return;
      }

    }
    std::cout << "HDCG: Pcg not converged"<<std::endl;
    return ;
  }
  
  // The Pcg routine is common to all, but the various matrices differ from derived 
  // implementation to derived implmentation
  void operator() (const Field &src, Field &psi){

    psi.Checkerboard() = src.Checkerboard();
    grid             = src.Grid();

    RealD f;
    RealD rtzp,rtz,a,d,b;
    RealD rptzp;
    RealD tn;
    RealD guess = norm2(psi);
    RealD ssq   = norm2(src);
    RealD rsq   = ssq*Tolerance*Tolerance;
    
    /////////////////////////////
    // Set up history vectors
    /////////////////////////////
    std::vector<Field> p  (mmax,grid);
    std::vector<Field> mmp(mmax,grid);
    std::vector<RealD> pAp(mmax);

    Field x  (grid);
    Field z  (grid);
    Field tmp(grid);
    Field r  (grid);
    Field mu (grid);
  
    //////////////////////////
    // x0 = Vstart -- possibly modify guess
    //////////////////////////
    x=Zero();
    Vstart(x,src);

    // r0 = b -A x0
    _FineLinop.HermOp(x,mmp[0]); // Fine operator
    axpy (r, -1.0,mmp[0], src);    // Recomputes r=src-Ax0

    //////////////////////////////////
    // Compute z = M1 r
    //////////////////////////////////
    PcgM1(r,z);
    rtzp =real(innerProduct(r,z));

    ///////////////////////////////////////
    // Solve for Mss mu = P A z and set p = z-mu
    ///////////////////////////////////////
    PcgM2(z,p[0]);

    for (int k=0;k<=MaxIterations;k++){
    
      int peri_k  = k % mmax;
      int peri_kp = (k+1) % mmax;

      rtz=rtzp;
      d= PcgM3(p[peri_k],mmp[peri_k]);
      a = rtz/d;
    
      // Memorise this
      pAp[peri_k] = d;
      std::cout << GridLogMessage << " pCG d "<< d<<std::endl;

      axpy(x,a,p[peri_k],x);
      //      std::cout << GridLogMessage << " pCG x "<< norm2(x)<<std::endl;
      RealD rn = axpy_norm(r,-a,mmp[peri_k],r);

      std::cout << GridLogMessage << " pCG rn "<< rn<<std::endl;

      // Compute z = M x
      PcgM1(r,z);
      //      std::cout << GridLogMessage << " pCG z "<< norm2(z)<<std::endl;

      rtzp =real(innerProduct(r,z));
      std::cout << GridLogMessage << " pCG rtzp "<<rtzp<<std::endl;
      //      std::cout << GridLogMessage << " pCG r "<<norm2(r)<<std::endl;

      PcgM2(z,mu); // ADEF-2 this is identity. Axpy possible to eliminate

      //      std::cout << GridLogMessage << " pCG mu "<<norm2(mu)<<std::endl;
      
      p[peri_kp]=mu;

      //      std::cout << GridLogMessage << " pCG p[peri_kp] "<<norm2(p[peri_kp])<<std::endl;

      // Standard search direction  p -> z + b p 
      b = (rtzp)/rtz;
      std::cout << GridLogMessage << " pCG b "<< b<<std::endl;

      int northog;
      //    northog     = (peri_kp==0)?1:peri_kp; // This is the fCG(mmax) algorithm
      northog     = (k>mmax-1)?(mmax-1):k;        // This is the fCG-Tr(mmax-1) algorithm
    
      for(int back=0; back < northog; back++){
	int peri_back = (k-back)%mmax;
	RealD pbApk= real(innerProduct(mmp[peri_back],p[peri_kp]));
	RealD beta = -pbApk/pAp[peri_back];
	axpy(p[peri_kp],beta,p[peri_back],p[peri_kp]);
      }
      //      std::cout << GridLogMessage << " pCG p[peri_kp] orthog "<< norm2(p[peri_kp])<<std::endl;

      RealD rrn=sqrt(rn/ssq);
      std::cout<<GridLogMessage<<"TwoLevelfPcg: k= "<<k<<" residual = "<<rrn<<std::endl;

      // Stopping condition
      if ( rn <= rsq ) { 

	_FineLinop.HermOp(x,mmp[0]); // Shouldn't this be something else?
	axpy(tmp,-1.0,src,mmp[0]);
	
	RealD psinorm = sqrt(norm2(x));
	RealD srcnorm = sqrt(norm2(src));
	RealD tmpnorm = sqrt(norm2(tmp));
	RealD true_residual = tmpnorm/srcnorm;
	std::cout<<GridLogMessage<<"TwoLevelfPcg:   true residual is "<<true_residual<<std::endl;
	std::cout<<GridLogMessage<<"TwoLevelfPcg: target residual was"<<Tolerance<<std::endl;
	return;
      }
    }
    // Non-convergence
    assert(0);
  }

 public:

  virtual void PcgM1(Field & in, Field & out)
  {
    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]

    Field tmp(grid);
    Field Min(grid);
    CoarseField PleftProj(coarsegrid);
    CoarseField PleftMss_proj(coarsegrid);

    _Smoother(in,Min);

    _FineLinop.HermOp(Min,out);
    axpy(tmp,-1.0,out,in);          // tmp  = in - A Min

    _Aggregates.ProjectToSubspace(PleftProj,tmp);     
    _CoarseSolver(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    _Aggregates.PromoteFromSubspace(PleftMss_proj,tmp);// tmp = Q[in - A Min]  

    axpy(out,1.0,Min,tmp); // Min+tmp

  }

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
    Field r(grid);
    Field mmp(grid);
    CoarseField PleftProj(coarsegrid);
    CoarseField PleftMss_proj(coarsegrid);

    _Aggregates.ProjectToSubspace(PleftProj,src);     
    _CoarseSolverPrecise(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    _Aggregates.PromoteFromSubspace(PleftMss_proj,x);  

  }

  /////////////////////////////////////////////////////////////////////
  // Only Def1 has non-trivial Vout.
  /////////////////////////////////////////////////////////////////////
  virtual void   Vout  (Field & in, Field & out,Field & src){
    out = in;
  }
};

NAMESPACE_END(Grid);
#endif
