/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/MinimalResidual.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_MINIMAL_RESIDUAL_H
#define GRID_MINIMAL_RESIDUAL_H

namespace Grid {

/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

template <class Field>
class MinimalResidual : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge;  // throw an assert when the MR fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the MR took to finish. Filled in upon completion
  
  MinimalResidual(RealD tol, Integer maxit, bool err_on_no_conv = true)
      : Tolerance(tol),
        MaxIterations(maxit),
        ErrorOnNoConverge(err_on_no_conv){};

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src,
                  Field &psi) {
    psi.checkerboard = src.checkerboard; // Check
    conformable(psi, src);

    /////
    RealD cp, c, a, d, b, ssq, qq, b_pred;

    Field p(src);
    Field mmp(src);
    Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);
    /////

    Field p {src};
    Field matrixTimesPsi {src};
    Field r {src};

    RealD alpha {};

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    Linop.HermOp(psi, matrixTimesPsi);

    r = src - matrixTimesPsi;

    Linop.HermOp(r, p);

    alpha = innerProduct(p,r) / innerProduct(p,p);
    psi = psi + alpha * r;
    r   = r   - alpha * p;

    Linop.HermOp(r, p);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    // RealD cp, c, a, d, b, ssq, qq, b_pred;

    Field p(src);
    Field matrixTimesPsi(src);
    // Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    
    Linop.HermOpAndNorm(psi, matrixTimesPsi, d, b);
    

    r = src - matrixTimesPsi;
    p = matrixTimesPsi;

    a = norm2(p);
    cp = a;
    ssq = norm2(src);

    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:   src " << ssq << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:   matrixTimesPsi " << b << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:  cp,r " << cp << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:     p " << a << std::endl;

    RealD rsq = Tolerance * Tolerance * ssq;

    // Check if guess is really REALLY good :)
    if (cp <= rsq) {
      return;
    }

    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual: k=0 residual " << cp << " target " << rsq
              << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {
      c = cp;

      MatrixTimer.Start();
      Linop.HermOpAndNorm(p, matrixTimesPsi, d, qq);
      MatrixTimer.Stop();

      LinalgTimer.Start();
      //  RealD    qqck = norm2(matrixTimesPsi);
      //  ComplexD dck  = innerProduct(p,matrixTimesPsi);

      a = c / d;
      b_pred = a * (a * qq - d) / c;

      cp = axpy_norm(r, -a, matrixTimesPsi, r);
      b = cp / c;

      // Fuse these loops ; should be really easy
      psi = a * p + psi;
      p = p * b + r;

      LinalgTimer.Stop();
      std::cout << GridLogIterative << "MinimalResidual: Iteration " << k
                << " residual " << cp << " target " << rsq << std::endl;

      // Stopping condition
      if (cp <= rsq) {
        SolverTimer.Stop();
        Linop.HermOpAndNorm(psi, matrixTimesPsi, d, qq);
        p = matrixTimesPsi - src;

        RealD matrixTimesPsiNorm = sqrt(norm2(matrixTimesPsi));
        RealD psinorm = sqrt(norm2(psi));
        RealD srcnorm = sqrt(norm2(src));
        RealD resnorm = sqrt(norm2(p));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage
                  << "MinimalResidual: Converged on iteration " << k << std::endl;
        std::cout << GridLogMessage << "Computed residual " << sqrt(cp / ssq)
                  << " true residual " << true_residual << " target "
                  << Tolerance << std::endl;
        std::cout << GridLogMessage << "Time elapsed: Iterations "
                  << SolverTimer.Elapsed() << " Matrix  "
                  << MatrixTimer.Elapsed() << " Linalg "
                  << LinalgTimer.Elapsed();
        std::cout << std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);
	IterationsToComplete = k;	
        return;
      }
    }
    std::cout << GridLogMessage << "MinimalResidual did NOT converge"
              << std::endl;
    if (ErrorOnNoConverge) assert(0);
    IterationsToComplete = k;
  }

  //! Minimal-residual (MR) algorithm for a generic Linear Operator
  /*! \ingroup invert
   * This subroutine uses the Minimal Residual (MR) algorithm to determine
   * the solution of the set of linear equations. Here we allow M to be nonhermitian.
   *
   *    M . Psi  =  src
   *
   * Algorithm:
   *
   *  Psi[0]                                      Argument
   *  r[0]    :=  src  -  M . Psi[0] ;            Initial residual
   *  IF |r[0]| <= RsdCG |src| THEN RETURN;       Converged?
   *  FOR k FROM 1 TO MaxCG DO                    MR iterations
   *      a[k-1]  := <M.r[k-1],r[k-1]> / <M.r[k-1],M.r[k-1]> ;
   *      ap[k-1] := MRovpar * a[k] ;             Overrelaxtion step
   *      Psi[k]  += ap[k-1] r[k-1] ;                   New solution vector
   *      r[k]    -= ap[k-1] A . r[k-1] ;         New residual
   *      IF |r[k]| <= RsdCG |src| THEN RETURN;   Converged?

   * Arguments:

   *  \param M       Linear Operator             (Read)
   *  \param src     Source                      (Read)
   *  \param psi     Solution                    (Modify)
   *  \param RsdCG   MR residual accuracy        (Read)
   *  \param MRovpar Overrelaxation parameter    (Read)
   *  \param MaxIterations   Maximum MR iterations       (Read)

   * Local Variables:

   *  r         Residual vector
   *  cp        | r[k] |**2
   *  c         | r[k-1] |**2
   *  k         MR iteration counter
   *  a         a[k]
   *  d         < M.r[k], M.r[k] >
   *  R_Aux     Temporary for  M.Psi
   *  Mr        Temporary for  M.r

   * Global Variables:

   *  MaxIterations       Maximum number of MR iterations allowed
   *  RsdCG       Maximum acceptable MR residual (relative to source)
   *
   * Subroutines:
   *
   *  M           Apply matrix to vector
   *
   * @{
   */

  // TODO: figure out what isign from chroma is supposed to do
  void tmpImplFromChroma(LinearOperatorBase<Field> &Linop, const Field &src,
                  Field &psi) {
    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    Complex a, c;
    Complex c;
    RealD d;

    Field Mr(src);
    Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD ssq = norm2(src); // flopcount.addSiteFlops(4*Nc*Ns,s); // stands for "source squared"
    RealD rsd_sq = Tolerance * Tolerance * ssq; // flopcount.addSiteFlops(4*Nc*Ns,s); // stands for "residual squared"

    /*  r[0]  :=  src - M . Psi[0] */
    /*  r  :=  M . Psi  */
    M(Mr, psi, isign); // flopcount.addFlops(M.nFlops());

    r = src - Mr; // flopcount.addSiteFlops(2*Nc*Ns,s);

    RealD cp = norm2(r); /*  Cp = |r[0]|^2 */ /* 2 Nc Ns  flops */ // flopcount.addSiteFlops(4*Nc*Ns, s);

    if (cp <= rsd_sq) { /*  IF |r[0]| <= Tolerance|src| THEN RETURN; */
      return;
    }

    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual: k=0 residual " << cp << " target " << rsq_sq << std::endl;

    /*  FOR k FROM 1 TO MaxIterations DO */
    auto k = 0;
    while( (k < MaxIterations) && (cp > rsd_sq) )
    {
      ++k;

      /*  a[k-1] := < M.r[k-1], r[k-1] >/ < M.r[k-1], M.r[k-1] > ; */

      M(Mr, r, isign); /*  Mr = M * r  */  // flopcount.addFlops(M.nFlops());

      c = innerProduct(Mr, r); /*  c = < M.r, r > */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      d = norm2(Mr); /*  d = | M.r | ** 2  */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      a = c / d;  /*  a = c / d */

      a = a * MRovpar; /*  a[k-1] *= MRovpar ; */


      psi = psi + r * a;  /*  Psi[k] += a[k-1] r[k-1] ; */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      r = r - Mr * a; /*  r[k] -= a[k-1] M . r[k-1] ; */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      cp = norm2(r); /*  cp  =  | r[k] |**2 */ // flopcount.addSiteFlops(4*Nc*Ns,s);

//    std::cout << "InvMR: k = " << k << "  cp = " << cp << endl;
    }

    IterationsToComplete = k;

    res.resid   = sqrt(cp);
    swatch.stop();
    std::cout << "InvMR: k = " << k << "  cp = " << cp << endl;
    // flopcount.report("invmr", swatch.getTimeInSeconds());

    // Compute the actual residual
    {
      M(Mr, psi, isign);
      RealD actual_res = norm2(src- Mr);
      res.resid = sqrt(actual_res);
    }

    if ( IterationsToComplete == MaxIterations )
      std::cerr << "Nonconvergence Warning" << endl;

    END_CODE();
    return res;

  }
};
}
#endif
