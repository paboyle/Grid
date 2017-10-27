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

template<class Field> class MinimalResidual : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge; // throw an assert when the MR fails to converge.
                          // Defaults true.
  RealD   Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; // Number of iterations the MR took to finish. Filled in upon completion

  MinimalResidual(RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Tolerance(tol),
      MaxIterations(maxit),
      ErrorOnNoConverge(err_on_no_conv){};

  //! Minimal-residual (MR) algorithm for a generic Linear Operator
  /*! \ingroup invert
   * This subroutine uses the Minimal Residual (MR) algorithm to determine
   * the solution of the set of linear equations. Here we allow M to be
   nonhermitian.
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

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {

    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    Complex a, c;
    RealD   d;

    Field Mr(src);
    Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD ssq = norm2(src); // flopcount.addSiteFlops(4*Nc*Ns,s); // stands for "source squared"
    RealD rsd_sq = Tolerance * Tolerance * ssq; // flopcount.addSiteFlops(4*Nc*Ns,s); //
                                     // stands for "residual squared"

    /*  r[0]  :=  src - M . Psi[0] */
    /*  r  :=  M . Psi  */
    // M(Mr, psi, isign); // flopcount.addFlops(M.nFlops());
    Linop.Op(psi, Mr); // flopcount.addFlops(M.nFlops());

    r = src - Mr; // flopcount.addSiteFlops(2*Nc*Ns,s);

    RealD cp = norm2(r);   /*  Cp = |r[0]|^2 */
      /* 2 Nc Ns  flops */ // flopcount.addSiteFlops(4*Nc*Ns, s);
    // auto cp = norm2(r); /*  Cp = |r[0]|^2 */ /* 2 Nc Ns  flops */ //
    // flopcount.addSiteFlops(4*Nc*Ns, s);

    if(cp <= rsd_sq) { /*  IF |r[0]| <= Tolerance|src| THEN RETURN; */
      return;
    }

    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual: k=0 residual " << cp << " target " << rsd_sq << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    auto k = 0;
    while((k < MaxIterations) && (cp > rsd_sq)) {
      ++k;

      /*  a[k-1] := < M.r[k-1], r[k-1] >/ < M.r[k-1], M.r[k-1] > ; */

      MatrixTimer.Start();
      // M(Mr, r, isign); /*  Mr = M * r  */  // flopcount.addFlops(M.nFlops());
      Linop.Op(r, Mr); /*  Mr = M * r  */ // flopcount.addFlops(M.nFlops());
      MatrixTimer.Stop();

      LinalgTimer.Start();

      c = innerProduct(Mr, r); /*  c = < M.r, r > */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      d = norm2(Mr); /*  d = | M.r | ** 2  */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      a = c / d;

      // a = a * MRovpar; /*  a[k-1] *= MRovpar ; */

      psi = psi + r * a; /*  Psi[k] += a[k-1] r[k-1] ; */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      r = r - Mr * a; /*  r[k] -= a[k-1] M . r[k-1] ; */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      cp = norm2(r); /*  cp  =  | r[k] |**2 */ // flopcount.addSiteFlops(4*Nc*Ns,s);

      LinalgTimer.Stop();

      std::cout << GridLogIterative << "MinimalResidual: Iteration " << k
                << " residual " << cp << " target " << rsd_sq << std::endl;
    }
    SolverTimer.Stop();

    IterationsToComplete = k;

    // res.resid   = sqrt(cp);
    std::cout << "InvMR: k = " << k << "  cp = " << cp << std::endl;
    // flopcount.report("invmr", swatch.getTimeInSeconds());

    std::cout << GridLogMessage << "MinimalResidual Converged on iteration " << k << std::endl;
    std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq)<<std::endl;
    // std::cout << GridLogMessage << "\tTrue residual " << true_residual<<std::endl;
    // std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

    std::cout << GridLogMessage << "Time breakdown "<<std::endl;
    std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed() <<std::endl;
    std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
    std::cout << GridLogMessage << "\tLinalg     " << LinalgTimer.Elapsed() <<std::endl;

    // Compute the actual residual
    {
      // M(Mr, psi, isign);
      Linop.Op(psi, Mr);
      Field tmp = src - Mr;
      // RealD actual_res = norm2(src-Mr);
      RealD actual_res = norm2(tmp);
      // res.resid = sqrt(actual_res);
    }

    if(IterationsToComplete == MaxIterations)
      std::cerr << "Nonconvergence Warning" << std::endl;

    // return res;
  }
};
} // namespace Grid
#endif
