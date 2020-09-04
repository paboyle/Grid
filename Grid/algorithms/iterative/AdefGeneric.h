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

NAMESPACE_BEGIN(Grid);
  /*
   * Compared to Tang-2009:  P=Pleft. P^T = PRight Q=MssInv. 
   * Script A = SolverMatrix 
   * Script P = Preconditioner
   *
   * Deflation methods considered
   *      -- Solve P A x = P b        [ like Luscher ]
   * DEF-1        M P A x = M P b     [i.e. left precon]
   * DEF-2        P^T M A x = P^T M b
   * ADEF-1       Preconditioner = M P + Q      [ Q + M + M A Q]
   * ADEF-2       Preconditioner = P^T M + Q
   * BNN          Preconditioner = P^T M P + Q
   * BNN2         Preconditioner = M P + P^TM +Q - M P A M 
   * 
   * Implement ADEF-2
   *
   * Vstart = P^Tx + Qb
   * M1 = P^TM + Q
   * M2=M3=1
   * Vout = x
   */


template<class Field, class CoarseField, class Aggregates>
class TwoLevelFlexiblePcg : public LinearFunction<Field>
{
 public:

  int verbose;

  RealD   Tolerance;
  Integer MaxIterations;
  const int mmax = 4;
  GridBase *FineGrid;
  GridBase *CoarseGrid;

  LinearOperatorBase<Field>   &_Linop;
  LinearFunction<Field>     &_Smoother;
  LinearFunction<CoarseField> &_CoarseSolver;
  Aggregates                  &_Aggregates;
  
  // more most opertor functions
  TwoLevelFlexiblePcg(RealD tol,
		      Integer maxit,
		      LinearOperatorBase<Field> *Linop,
		      LinearFunction<Field>   *Smoother,
		      LinearFunction<CoarseField> *CoarseSolver,
		      Aggregates *AggP
		      ) : 
  Tolerance(tol), 
    MaxIterations(maxit),
    _Linop(*Linop),
    _Smoother(*Smoother),
    _CoarseSolver(*CoarseSolver),
    _Aggregates(*AggP)
  { 
    CoarseGrid=_Aggregates.CoarseGrid;
    FineGrid=_Aggregates.FineGrid;
    verbose=0;
  };

  // The Pcg routine is common to all, but the various matrices differ from derived 
  // implementation to derived implmentation
  void operator() (const Field &src, Field &psi){

    psi.Checkerboard() = src.Checkerboard();

    RealD rtzp,rtz,a,d,b;
    //    RealD rptzp;
    //    RealD tn;
    RealD guess = norm2(psi);
    RealD ssq   = norm2(src);
    RealD rsq   = ssq*Tolerance*Tolerance;
    
    /////////////////////////////
    // Set up history vectors
    /////////////////////////////
    std::vector<Field> p  (mmax,FineGrid);
    std::vector<Field> mmp(mmax,FineGrid);
    std::vector<RealD> pAp(mmax);

    Field x  (FineGrid); x = psi;
    Field z  (FineGrid);
    Field tmp(FineGrid);
    Field r  (FineGrid);
    Field mu (FineGrid);
  
    //////////////////////////
    // x0 = Vstart -- possibly modify guess
    //////////////////////////
    x=src;
    Vstart(x,src);

    // r0 = b -A x0
    _Linop.HermOp(x,mmp[0]); // Shouldn't this be something else?
    axpy (r, -1.0,mmp[0], src);    // Recomputes r=src-Ax0

    //////////////////////////////////
    // Compute z = M1 x
    //////////////////////////////////
    M1(r,z);
    rtzp =real(innerProduct(r,z));

    ///////////////////////////////////////
    // Solve for Mss mu = P A z and set p = z-mu
    // Def2: p = 1 - Q Az = Pright z 
    // Other algos M2 is trivial
    ///////////////////////////////////////
    M2(z,p[0]);

    for (int k=0;k<=MaxIterations;k++){
    
      int peri_k  = k % mmax;
      int peri_kp = (k+1) % mmax;

      rtz=rtzp;
      d= M3(p[peri_k],mmp[peri_k]);
      a = rtz/d;
    
      // Memorise this
      pAp[peri_k] = d;

      axpy(x,a,p[peri_k],x);
      RealD rn = axpy_norm(r,-a,mmp[peri_k],r);

      // Compute z = M x
      M1(r,z);

      rtzp =real(innerProduct(r,z));

      M2(z,mu); // ADEF-2 this is identity. Axpy possible to eliminate

      p[peri_kp]=mu;

      // Standard search direction  p -> z + b p    ; b = 
      b = (rtzp)/rtz;

      int northog;
      //    northog     = (peri_kp==0)?1:peri_kp; // This is the fCG(mmax) algorithm
      northog     = (k>mmax-1)?(mmax-1):k;        // This is the fCG-Tr(mmax-1) algorithm
    
      for(int back=0; back < northog; back++){
	int peri_back = (k-back)%mmax;
	RealD pbApk= real(innerProduct(mmp[peri_back],p[peri_kp]));
	RealD beta = -pbApk/pAp[peri_back];
	axpy(p[peri_kp],beta,p[peri_back],p[peri_kp]);
      }

      RealD rrn=sqrt(rn/ssq);
      std::cout<<GridLogMessage<<"TwoLevelfPcg: k= "<<k<<" residual = "<<rrn<<std::endl;

      // Stopping condition
      if ( rn <= rsq ) { 

	_Linop.HermOp(x,mmp[0]); // Shouldn't this be something else?
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

  virtual void M1(Field & in, Field & out) 
  {// the smoother
    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    Field tmp(FineGrid);
    Field Min(FineGrid);

    CoarseField PleftProj(CoarseGrid);
    CoarseField PleftMss_proj(CoarseGrid);

    _Smoother(in,Min); // Smoother call

    _Linop.HermOp(Min,out);
    axpy(tmp,-1.0,out,in);          // tmp  = in - A Min

    _Aggregates.ProjectToSubspace(PleftProj,tmp);     
    _CoarseSolver(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    _Aggregates.PromoteFromSubspace(PleftMss_proj,tmp);// tmp = Q[in - A Min]  
    axpy(out,1.0,Min,tmp); // Min+tmp
  }

  virtual void M2(const Field & in, Field & out) 
  {
    out=in;
  }

  virtual RealD M3(const Field & p, Field & mmp)
  {
    double d,dd;
    _Linop.HermOpAndNorm(p,mmp,d,dd);
    return dd;
  }

  virtual void Vstart(Field & x,const Field & src)
  {
    //case PcgDef2:
    //case PcgAdef2: 
    //case PcgAdef2f:
    //case PcgV11f:
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
    Field r(FineGrid);
    Field mmp(FineGrid);

    CoarseField PleftProj(CoarseGrid);
    CoarseField PleftMss_proj(CoarseGrid);
    
    _Linop.HermOp(x,mmp);
    axpy (r, -1.0, mmp, src);        // r_{-1} = src - A x
    _Aggregates.ProjectToSubspace(PleftProj,r);     
    _CoarseSolver(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    _Aggregates.PromoteFromSubspace(PleftMss_proj,mmp);  
    x=x+mmp;

  }

  /////////////////////////////////////////////////////////////////////
  // Only Def1 has non-trivial Vout. Override in Def1
  /////////////////////////////////////////////////////////////////////
  virtual void   Vout  (Field & in, Field & out,Field & src){
    out = in;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Pright and Pleft are common to all implementations
  ////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void Pright(Field & in,Field & out)
  {
    // P_R  = [ 1              0 ] 
    //        [ -Mss^-1 Msb    0 ] 
    Field in_sbar(FineGrid);

    CoarseField PleftProj(CoarseGrid);
    CoarseField PleftMss_proj(CoarseGrid);

    _Aggregates.ProjectToSubspace(PleftProj,in);     
    _Aggregates.PromoteFromSubspace(PleftProj,out);  
    axpy(in_sbar,-1.0,out,in);       // in_sbar = in - in_s 

    _Linop.HermOp(in_sbar,out);
    _Aggregates.ProjectToSubspace(PleftProj,out);           // Mssbar in_sbar  (project)

    _CoarseSolver(PleftProj,PleftMss_proj); // Mss^{-1} Mssbar 
    _Aggregates.PromoteFromSubspace(PleftMss_proj,out);     // 

    axpy(out,-1.0,out,in_sbar);     // in_sbar - Mss^{-1} Mssbar in_sbar
  }
  virtual void Pleft (Field & in,Field & out)
  {
    // P_L  = [ 1  -Mbs Mss^-1] 
    //        [ 0   0         ] 
    Field in_sbar(FineGrid);
    Field    tmp2(FineGrid);
    Field    Mtmp(FineGrid);

    CoarseField PleftProj(CoarseGrid);
    CoarseField PleftMss_proj(CoarseGrid);

    _Aggregates.ProjectToSubspace(PleftProj,in);     
    _Aggregates.PromoteFromSubspace(PleftProj,out);  
    axpy(in_sbar,-1.0,out,in);      // in_sbar = in - in_s

    _CoarseSolver(PleftProj,PleftMss_proj); // Mss^{-1} in_s
    _Aggregates.PromoteFromSubspace(PleftMss_proj,out);

    _Linop.HermOp(out,Mtmp);

    _Aggregates.ProjectToSubspace(PleftProj,Mtmp);      // Msbar s Mss^{-1}
    _Aggregates.PromoteFromSubspace(PleftProj,tmp2);

    axpy(out,-1.0,tmp2,Mtmp);
    axpy(out,-1.0,out,in_sbar);     // in_sbar - Msbars Mss^{-1} in_s
  }
};
NAMESPACE_END(Grid);

#endif
