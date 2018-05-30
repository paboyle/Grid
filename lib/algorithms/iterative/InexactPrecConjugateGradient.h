    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/InexactPrecConjugateGradient.h

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@phys.columbia.edu>

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
#ifndef GRID_INEXACT_PREC_CONJUGATE_GRADIENT_H_
#define GRID_INEXACT_PREC_CONJUGATE_GRADIENT_H_

namespace Grid {

//Inexact preconditioned CG based on Golub, Ye, SIAM J. Sci. Comput., 21(4), 13051320. 
//(https://pdfs.semanticscholar.org/d2a9/d5bab02146a7fe3a244677432d21e33a2d98.pdf)
template <class Field>
class InexactPreconditionedConjugateGradient : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  
  LinearOperatorBase<Field> &Prec;

  InexactPreconditionedConjugateGradient(LinearOperatorBase<Field> &_Prec, RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Prec(_Prec),
    Tolerance(tol),
    MaxIterations(maxit),
    ErrorOnNoConverge(err_on_no_conv){};

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {

    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    Real ssq = norm2(src);
    RealD rsq = Tolerance * Tolerance * ssq; //inner stopping condition

    Field p(src);
    Field r(src);
    Field rnm1(src);
    Field mmp(src);
    Field z(src);

    //Initialize variables
    Linop.HermOp(psi, mmp);   
    r = src - mmp;
    
    Real cp = norm2(r);
    
    p = zero;
    Real alpha = 0, beta = 0;
    
    Real z_nm1_dot_r_nm1;

    int n;
    for(n=1; n <= MaxIterations; n++) {
      //Check stopping condition
      if (cp <= rsq) {
        Linop.HermOp(psi, mmp);
        r = mmp - src;

        RealD srcnorm = sqrt(norm2(src));
        RealD resnorm = sqrt(norm2(r));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage << "InexactPreconditionedConjugateGradient Converged on iteration " << n << std::endl;
        std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq)<<std::endl;
  	std::cout << GridLogMessage << "\tTrue residual " << true_residual<<std::endl;
  	std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);

  	IterationsToComplete = n;	

        return;
      }
      
      std::cout << GridLogIterative << std::setprecision(8)
		<< "InexactPreconditionedConjugateGradient: n=" << n << " residual " << cp << " target " << rsq << std::endl;

      //Apply preconditioner to current residual
      Prec.HermOp(r, z);

      //Update beta and store appropriate variables for next iteration
      Real z_n_dot_r_n  = sqrt(norm(innerProduct(z,r)));

      if(n>1){
	//  z^T_n ( r_n - r_{n-1} )
	//  -----------------------
	//     z^T_{n-1} r_{n-1}

	Real z_n_dot_r_nm1 = sqrt(norm(innerProduct(z,rnm1)));
	beta = ( z_n_dot_r_n - z_n_dot_r_nm1 ) / z_nm1_dot_r_nm1;
	std::cout << GridLogIterative << "beta " << beta << std::endl;
      }

      z_nm1_dot_r_nm1 = z_n_dot_r_n; //for next iteration
      rnm1 = r;

      axpy(p, beta, p, z); //p = beta * p + z

      //Compute alpha
      Linop.HermOp(p, mmp);      
      alpha = z_n_dot_r_n / sqrt(norm(innerProduct(p, mmp)));

      std::cout << GridLogIterative << "alpha " << alpha << std::endl;

      //Update residual and solution
      cp = axpy_norm(r, -alpha, mmp, r);
      axpy(psi, alpha, p, psi);
    }
    std::cout << GridLogMessage << "InexactPreconditionedConjugateGradient did NOT converge"
              << std::endl;

    if (ErrorOnNoConverge) assert(0);
    IterationsToComplete = n;
  }
};


template<class Field>
class PolynomialPreconditioner :  public LinearOperatorBase<Field> {
  Chebyshev<Field> Cheby;
  LinearOperatorBase<Field> &linop;
public:
  int InnerIterations;
  int order;
  PolynomialPreconditioner(LinearOperatorBase<Field> &_linop,RealD lo, RealD hi, int _order) 
    : linop(_linop), Cheby(lo,hi,_order,__InverseApproximation) 
    {
      InnerIterations=0;
      order = _order;
    };

  void OpDiag (const Field &in, Field &out){ assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp){ assert(0); }
  void Op     (const Field &in, Field &out){ assert(0); }
  void AdjOp  (const Field &in, Field &out){ assert(0); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ 
    HermOp(in,out);
    n1 = 0; n2 = norm2(out); 
  }
  void HermOp(const Field &in, Field &out){ 
    Cheby(linop,in,out); 
    InnerIterations+=order;
  }
};

template<class Field>
class DoNothingLinearOperator :  public LinearOperatorBase<Field> {
public:
  void OpDiag (const Field &in, Field &out){ assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp){ assert(0); }
  void Op     (const Field &in, Field &out){ assert(0); }
  void AdjOp  (const Field &in, Field &out){ assert(0); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ out = in; n1 = 0; n2 = norm2(out); }
  void HermOp(const Field &in, Field &out){ out = in; }
};

template<class Field>
class FixedIterConjugateGradientPreconditioner :  public LinearOperatorBase<Field> {
public:
  LinearOperatorBase<Field> &linop;
  ConjugateGradient<Field> CG;

  FixedIterConjugateGradientPreconditioner (LinearOperatorBase<Field> &_linop, Integer _iter): linop(_linop), CG(1e-20, _iter){
    CG.ErrorOnNoConverge = false;
  }

  void OpDiag (const Field &in, Field &out){ assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp){ assert(0); }
  void Op     (const Field &in, Field &out){ assert(0); }
  void AdjOp  (const Field &in, Field &out){ assert(0); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ 
    out = zero;
    CG(linop,in,out);    
    n2 = norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    out = zero;
    CG(linop,in,out);    
  }
};

template<class Field>
class SloppyConjugateGradientPreconditioner :  public LinearOperatorBase<Field> {
public:
  LinearOperatorBase<Field> &linop;
  ConjugateGradient<Field> CG;
  int InnerIterations;

  SloppyConjugateGradientPreconditioner (LinearOperatorBase<Field> &_linop, Real _resid, Integer max_iter): linop(_linop), CG(_resid, max_iter), InnerIterations(0){
  }

  void ResetCounters(){ InnerIterations = 0; }

  void OpDiag (const Field &in, Field &out){ assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp){ assert(0); }
  void Op     (const Field &in, Field &out){ assert(0); }
  void AdjOp  (const Field &in, Field &out){ assert(0); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ 
    out = zero;
    CG(linop,in,out);    
    InnerIterations += CG.IterationsToComplete;
    n2 = norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    out = zero;
    CG(linop,in,out);
    InnerIterations += CG.IterationsToComplete;
  }
};


template<class FieldH, class FieldL>
class SloppyConjugateGradientLowerPrecPreconditioner :  public LinearOperatorBase<FieldH> {
public:
  LinearOperatorBase<FieldL> &linop;
  ConjugateGradient<FieldL> CG;
  GridBase* L_grid; //lower-prec Grid
  int InnerIterations;
  FieldL tmp_l1;
  FieldL tmp_l2;

  SloppyConjugateGradientLowerPrecPreconditioner (LinearOperatorBase<FieldL> &_linop, GridBase* _L_grid, Real _resid, Integer max_iter): linop(_linop), CG(_resid, max_iter), InnerIterations(0), L_grid(_L_grid), tmp_l1(_L_grid), tmp_l2(_L_grid){
    CG.ErrorOnNoConverge = false;
  }

  void ResetCounters(){ InnerIterations = 0; }

  void OpDiag (const FieldH &in, FieldH &out){ assert(0); }
  void OpDir  (const FieldH &in, FieldH &out,int dir,int disp){ assert(0); }
  void Op     (const FieldH &in, FieldH &out){ assert(0); }
  void AdjOp  (const FieldH &in, FieldH &out){ assert(0); }
  void HermOpAndNorm(const FieldH &in, FieldH &out,RealD &n1,RealD &n2){ 
    precisionChange(tmp_l1, in);
    tmp_l2 = zero;
    CG(linop,tmp_l1,tmp_l2);    
    InnerIterations += CG.IterationsToComplete;
    precisionChange(out, tmp_l2);
    n2 = norm2(tmp_l2);
  }
  void HermOp(const FieldH &in, FieldH &out){
    precisionChange(tmp_l1, in);
    tmp_l2 = zero;
    CG(linop,tmp_l1,tmp_l2);    
    InnerIterations += CG.IterationsToComplete;
    precisionChange(out, tmp_l2);
  }
};

}



#endif
