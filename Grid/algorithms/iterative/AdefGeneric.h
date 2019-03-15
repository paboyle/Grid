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

// abstract base
template<class Field, class CoarseField>
class TwoLevelFlexiblePcg : public LinearFunction<Field>
{
 public:
  int verbose;
  RealD   Tolerance;
  Integer MaxIterations;
  const int mmax = 5;
  GridBase *grid;
  GridBase *coarsegrid;

  LinearOperatorBase<Field>   *_Linop
  OperatorFunction<Field>     *_Smoother,
  LinearFunction<CoarseField> *_CoarseSolver;

  // Need somthing that knows how to get from Coarse to fine and back again
  
  // more most opertor functions
  TwoLevelFlexiblePcg(RealD tol,
		     Integer maxit,
		     LinearOperatorBase<Field> *Linop,
		     LinearOperatorBase<Field> *SmootherLinop,
		     OperatorFunction<Field>   *Smoother,
		     OperatorFunction<CoarseField>  CoarseLinop
		     ) : 
      Tolerance(tol), 
      MaxIterations(maxit),
      _Linop(Linop),
      _PreconditionerLinop(PrecLinop),
      _Preconditioner(Preconditioner)
  { 
    verbose=0;
  };

  // The Pcg routine is common to all, but the various matrices differ from derived 
  // implementation to derived implmentation
  void operator() (const Field &src, Field &psi){
  void operator() (const Field &src, Field &psi){

    psi.checkerboard = src.checkerboard;
    grid             = src._grid;

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

    Field x  (grid); x = psi;
    Field z  (grid);
    Field tmp(grid);
    Field r  (grid);
    Field mu (grid);
  
    //////////////////////////
    // x0 = Vstart -- possibly modify guess
    //////////////////////////
    x=src;
    Vstart(x,src);

    // r0 = b -A x0
    HermOp(x,mmp); // Shouldn't this be something else?
    axpy (r, -1.0,mmp[0], src);    // Recomputes r=src-Ax0

    //////////////////////////////////
    // Compute z = M1 x
    //////////////////////////////////
    M1(r,z,tmp,mp,SmootherMirs);
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
      d= M3(p[peri_k],mp,mmp[peri_k],tmp);
      a = rtz/d;
    
      // Memorise this
      pAp[peri_k] = d;

      axpy(x,a,p[peri_k],x);
      RealD rn = axpy_norm(r,-a,mmp[peri_k],r);

      // Compute z = M x
      M1(r,z,tmp,mp);

      rtzp =real(innerProduct(r,z));

      M2(z,mu); // ADEF-2 this is identity. Axpy possible to eliminate

      p[peri_kp]=p[peri_k];

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

	HermOp(x,mmp); // Shouldn't this be something else?
	axpy(tmp,-1.0,src,mmp[0]);
	
	RealD psinorm = sqrt(norm2(x));
	RealD srcnorm = sqrt(norm2(src));
	RealD tmpnorm = sqrt(norm2(tmp));
	RealD true_residual = tmpnorm/srcnorm;
	std::cout<<GridLogMessage<<"TwoLevelfPcg:   true residual is "<<true_residual<<std::endl;
	std::cout<<GridLogMessage<<"TwoLevelfPcg: target residual was"<<Tolerance<<std::endl;
	return k;
      }
    }
    // Non-convergence
    assert(0);
  }

 public:

  virtual void M(Field & in,Field & out,Field & tmp) {

  }

  virtual void M1(Field & in, Field & out) {// the smoother

    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    Field tmp(grid);
    Field Min(grid);

    PcgM(in,Min); // Smoother call

    HermOp(Min,out);
    axpy(tmp,-1.0,out,in);          // tmp  = in - A Min

    ProjectToSubspace(tmp,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    PromoteFromSubspace(PleftMss_proj,tmp);// tmp = Q[in - A Min]  
    axpy(out,1.0,Min,tmp); // Min+tmp
  }

  virtual void M2(const Field & in, Field & out) {
    out=in;
    // Must override for Def2 only
    //  case PcgDef2:
    //    Pright(in,out);
    //    break;
  }

  virtual RealD M3(const Field & p, Field & mmp){
    double d,dd;
    HermOpAndNorm(p,mmp,d,dd);
    return dd;
    // Must override for Def1 only
    //  case PcgDef1:
    //    d=linop_d->Mprec(p,mmp,tmp,0,1);// Dag no
    //      linop_d->Mprec(mmp,mp,tmp,1);// Dag yes
    //    Pleft(mp,mmp);
    //    d=real(linop_d->inner(p,mmp));
  }

  virtual void VstartDef2(Field & xconst Field & src){
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
    Field r(grid);
    Field mmp(grid);
    
    HermOp(x,mmp);
    axpy (r, -1.0, mmp, src);        // r_{-1} = src - A x
    ProjectToSubspace(r,PleftProj);     
    ApplyInverseCG(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,mmp);  
    x=x+mmp;

  }

  virtual void Vstart(Field & x,const Field & src){
    return;
  }

  /////////////////////////////////////////////////////////////////////
  // Only Def1 has non-trivial Vout. Override in Def1
  /////////////////////////////////////////////////////////////////////
  virtual void   Vout  (Field & in, Field & out,Field & src){
    out = in;
    //case PcgDef1:
    //    //Qb + PT x
    //    ProjectToSubspace(src,PleftProj);     
    //    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    //    PromoteFromSubspace(PleftMss_proj,tmp);  
    //    
    //    Pright(in,out);
    //    
    //    linop_d->axpy(out,tmp,out,1.0);
    //    break;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Pright and Pleft are common to all implementations
  ////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void Pright(Field & in,Field & out){
    // P_R  = [ 1              0 ] 
    //        [ -Mss^-1 Msb    0 ] 
    Field in_sbar(grid);

    ProjectToSubspace(in,PleftProj);     
    PromoteFromSubspace(PleftProj,out);  
    axpy(in_sbar,-1.0,out,in);       // in_sbar = in - in_s 

    HermOp(in_sbar,out);
    ProjectToSubspace(out,PleftProj);           // Mssbar in_sbar  (project)

    ApplyInverse     (PleftProj,PleftMss_proj); // Mss^{-1} Mssbar 
    PromoteFromSubspace(PleftMss_proj,out);     // 

    axpy(out,-1.0,out,in_sbar);     // in_sbar - Mss^{-1} Mssbar in_sbar
  }
  virtual void Pleft (Field & in,Field & out){
    // P_L  = [ 1  -Mbs Mss^-1] 
    //        [ 0   0         ] 
    Field in_sbar(grid);
    Field    tmp2(grid);
    Field    Mtmp(grid);

    ProjectToSubspace(in,PleftProj);     
    PromoteFromSubspace(PleftProj,out);  
    axpy(in_sbar,-1.0,out,in);      // in_sbar = in - in_s

    ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} in_s
    PromoteFromSubspace(PleftMss_proj,out);

    HermOp(out,Mtmp);

    ProjectToSubspace(Mtmp,PleftProj);      // Msbar s Mss^{-1}
    PromoteFromSubspace(PleftProj,tmp2);

    axpy(out,-1.0,tmp2,Mtmp);
    axpy(out,-1.0,out,in_sbar);     // in_sbar - Msbars Mss^{-1} in_s
  }
}

template<class Field>
class TwoLevelFlexiblePcgADef2 : public TwoLevelFlexiblePcg<Field> {
 public:
  virtual void M(Field & in,Field & out,Field & tmp){

  } 
  virtual void M1(Field & in, Field & out,Field & tmp,Field & mp){

  }
  virtual void M2(Field & in, Field & out){

  }
  virtual RealD M3(Field & p, Field & mp,Field & mmp, Field & tmp){

  }
  virtual void Vstart(Field & in, Field & src, Field & r, Field & mp, Field & mmp, Field & tmp){

  }
}
/*
template<class Field>
class TwoLevelFlexiblePcgAD : public TwoLevelFlexiblePcg<Field> {
 public:
  virtual void M(Field & in,Field & out,Field & tmp); 
  virtual void M1(Field & in, Field & out,Field & tmp,Field & mp);
  virtual void M2(Field & in, Field & out);
  virtual RealD M3(Field & p, Field & mp,Field & mmp, Field & tmp);
  virtual void Vstart(Field & in, Field & src, Field & r, Field & mp, Field & mmp, Field & tmp);
}

template<class Field>
class TwoLevelFlexiblePcgDef1 : public TwoLevelFlexiblePcg<Field> {
 public:
  virtual void M(Field & in,Field & out,Field & tmp); 
  virtual void M1(Field & in, Field & out,Field & tmp,Field & mp);
  virtual void M2(Field & in, Field & out);
  virtual RealD M3(Field & p, Field & mp,Field & mmp, Field & tmp);
  virtual void Vstart(Field & in, Field & src, Field & r, Field & mp, Field & mmp, Field & tmp);
  virtual void   Vout  (Field & in, Field & out,Field & src,Field & tmp);
}

template<class Field>
class TwoLevelFlexiblePcgDef2 : public TwoLevelFlexiblePcg<Field> {
 public:
  virtual void M(Field & in,Field & out,Field & tmp); 
  virtual void M1(Field & in, Field & out,Field & tmp,Field & mp);
  virtual void M2(Field & in, Field & out);
  virtual RealD M3(Field & p, Field & mp,Field & mmp, Field & tmp);
  virtual void Vstart(Field & in, Field & src, Field & r, Field & mp, Field & mmp, Field & tmp);
}

template<class Field>
class TwoLevelFlexiblePcgV11: public TwoLevelFlexiblePcg<Field> {
 public:
  virtual void M(Field & in,Field & out,Field & tmp); 
  virtual void M1(Field & in, Field & out,Field & tmp,Field & mp);
  virtual void M2(Field & in, Field & out);
  virtual RealD M3(Field & p, Field & mp,Field & mmp, Field & tmp);
  virtual void Vstart(Field & in, Field & src, Field & r, Field & mp, Field & mmp, Field & tmp);
}
*/
#endif
