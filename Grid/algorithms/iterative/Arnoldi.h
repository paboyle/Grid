/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./lib/algorithms/iterative/Arnoldi.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Patrick Oare <poare@bnl.gov>

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
#ifndef GRID_ARNOLDI_H
#define GRID_ARNOLDI_H

NAMESPACE_BEGIN(Grid); 

/**
 * Implementation of the Arnoldi algorithm.
 */
template<class Field> 
class Arnoldi {

  private:
  
    std::string cname = std::string("Arnoldi");
    int MaxIter;   // Max iterations
    RealD Tolerance;
    RealD ssq;
    RealD rtol;
    int Nm;           // Number of basis vectors to track (equals MaxIter if no restart)
    int Nk;           // Number of basis vectors to keep every restart (equals -1 if no restart)
    int Nstop;       // Stop after converging Nstop eigenvectors.

    LinearOperatorBase<Field> &Linop;
    GridBase *Grid;

    RealD approxLambdaMax;
    RealD beta_k;
    Field f;
    std::vector<Field> basis;               // orthonormal Arnoldi basis
    Eigen::MatrixXcd Hess;                  // Hessenberg matrix of size Nbasis (after construction)
    Eigen::MatrixXcd Qt;                    // Transpose of basis rotation which projects out high modes.

    Eigen::VectorXcd evals;                 // evals of Hess
    Eigen::MatrixXcd littleEvecs;           // Nm x Nm evecs matrix
    std::vector<Field> evecs;               // Vector of evec fields

    RitzFilter ritzFilter;                        // how to sort evals

  public:       

    Arnoldi(LinearOperatorBase<Field> &_Linop, GridBase *_Grid, RealD _Tolerance, RitzFilter filter = EvalReSmall)
      : Linop(_Linop), Grid(_Grid), Tolerance(_Tolerance), ritzFilter(filter), f(_Grid), MaxIter(-1), Nm(-1), Nk(-1), 
          Nstop (-1), evals (0), evecs (), ssq (0.0), rtol (0.0), beta_k (0.0), approxLambdaMax (0.0)
    {
      f = Zero();
    };

    /**
     * Runs the Arnoldi loop with(out) implicit restarting. For each iteration:
     *   - Runs an Arnoldi step.
     *   - Computes the eigensystem of the Hessenberg matrix.
     *   - Performs implicit restarting.
     */
    void operator()(const Field& v0, int _maxIter, int _Nm, int _Nk, int _Nstop, bool doubleOrthog = false) {
      MaxIter = _maxIter;
      Nm = _Nm; Nk = _Nk;
      Nstop = _Nstop;
      
      ssq = norm2(v0);
      RealD approxLambdaMax = approxMaxEval(v0);
      rtol = Tolerance * approxLambdaMax;

      ComplexComparator compareComplex (ritzFilter);
      std::cout << GridLogMessage << "Comparing Ritz values with: " << ritzFilter << std::endl;

      int start = 1;
      Field startVec = v0;
      littleEvecs = Eigen::MatrixXcd::Zero(Nm, Nm);
      for (int i = 0; i < MaxIter; i++) {
        std::cout << GridLogMessage << "Restart Iteration " << i << std::endl;

        // Perform Arnoldi steps to compute Krylov basis and Rayleigh quotient (Hess)
        arnoldiIteration(startVec, Nm, start, doubleOrthog);
        startVec = f;

        // compute eigensystem and sort evals
        // compute_eigensystem();
        compute_eigensystem(Hess);
        std::cout << GridLogMessage << "Eigenvalues after Arnoldi step: " << std::endl << evals << std::endl;

        std::sort(evals.begin(), evals.end(), compareComplex);
        std::cout << GridLogMessage << "Ritz values after sorting (first Nk preserved): " << std::endl << evals << std::endl;
        // SU(N)::tepidConfiguration

        // Implicit restart to de-weight unwanted eigenvalues
        implicitRestart(_Nm, _Nk);      // probably can delete _Nm and _Nk from function args
        start = Nk;

        // check convergence and return if needed.
        int Nconv = converged();
        std::cout << GridLogMessage << "Number of evecs converged: " << Nconv << std::endl;
        if (Nconv >= Nstop || i == MaxIter - 1) {
          std::cout << GridLogMessage << "Converged with " << Nconv << " / " << Nstop << " eigenvectors on iteration " 
                        << i << "." << std::endl;
          basisRotate(evecs, Qt, 0, Nk, 0, Nk, Nm);
          std::cout << GridLogMessage << "Eigenvalues [first " << Nconv << " converged]: " << std::endl << evals << std::endl;
          return;
        }
      }      
    }

    /**
     * Approximates the maximum eigenvalue of Linop.Op to normalize the residual and test for convergence. 
     * 
     * Parameters
     * ----------
     * Field& v0
     *  Source field to start with. Must have non-zero norm.
     * int MAX_ITER (default = 50)
     *  Maximum number of iterations for power approximation. 
     * 
     * Returns
     * -------
     * RealD lamApprox
     *  Approximation of largest eigenvalue. 
     */
    RealD approxMaxEval(const Field& v0, int MAX_ITER = 50) {
      assert (norm2(v0) > 1e-8);                        // must have relatively large source norm to start
      RealD lamApprox = 0.0;
      RealD denom = 1.0; RealD num = 1.0;
      Field v0cp (Grid); Field tmp (Grid);
      v0cp = v0;
      denom = std::sqrt(norm2(v0cp));
      for (int i = 0; i < MAX_ITER; i++) {
        Linop.Op(v0cp, tmp);                               // CAREFUL: do not do Op(tmp, tmp)
        v0cp = tmp;
        num = std::sqrt(norm2(v0cp));                      // num = |A^{n+1} v0|
        lamApprox = num / denom;                           // lam = |A^{n+1} v0| / |A^n v0|
        std::cout << GridLogDebug << "Approx for max eval: " << lamApprox << std::endl;
        denom = num;                                       // denom = |A^{n} v0|
      }
      return lamApprox;
    }

    /**
     * Constructs the Arnoldi basis for the Krylov space K_n(D, src). (TODO make private)
     * 
     * Parameters
     * ----------
     * v0 : Field&
     *  Source to generate Krylov basis. 
     * Nm : int
     *  Final size of the basis desired. If the basis becomes complete before a basis of size Nm is constructed 
     *  (determined by relative tolerance Tolerance), stops iteration there. 
     * doubleOrthog : bool (default = false)
     *  Whether to double orthogonalize the basis (for numerical cancellations) or not. 
     * start        : int (default = 0)
     *  If non-zero, assumes part of the Arnoldi basis has already been constructed. 
     */
    void arnoldiIteration(const Field& v0, int Nm, int start = 1, bool doubleOrthog = false)
    {

      ComplexD coeff;
      Field w (Grid);           // A acting on last Krylov vector. 

      if (start == 1) {       // initialize everything that we need.
        RealD v0Norm = 1 / std::sqrt(ssq);
        basis.push_back(v0Norm * v0);                // normalized source

        Hess = Eigen::MatrixXcd::Zero(Nm, Nm);
        f = Zero();
      } else {
        assert( start == basis.size() );      // should be starting at the end of basis (start = Nk)
        Eigen::MatrixXcd HessCp = Hess;
        Hess = Eigen::MatrixXcd::Zero(Nm, Nm);
        Hess(Eigen::seqN(0, Nk), Eigen::seqN(0, Nk)) = HessCp;
      }

      // Construct next Arnoldi vector by normalizing w_i = Dv_i - \sum_j v_j h_{ji}
      for (int i = start - 1; i < Nm; i++) {

        Linop.Op(basis.back(), w);
        for (int j = 0; j < basis.size(); j++) {
          coeff = innerProduct(basis[j], w);       // coeff = h_{ij}. Note that since {vi} is ONB it's OK to subtract it off after. 
          Hess(j, i) = coeff;
          w -= coeff * basis[j];
        }

        if (doubleOrthog) {
          // TODO implement
        }

        // add w_i to the pile
        if (i < Nm - 1) {
          coeff = std::sqrt(norm2(w));
          Hess(i+1, i) = coeff;
          basis.push_back(
            (1.0/coeff) * w
          );
        }

        // after iterations, update f and beta_k = ||f||
        f = w;                                // make sure f is not normalized
        beta_k = std::sqrt(norm2(f));         // beta_k = ||f_k|| determines convergence.
      }

      std::cout << GridLogMessage << "|f|^2 after Arnoldi step = " << norm2(f) << std::endl;
      std::cout << GridLogDebug << "Computed Hessenberg matrix = " << std::endl << Hess << std::endl;

      return;
    }

    /**
     * Approximates the eigensystem of the linear operator by computing the eigensystem of 
     * the Hessenberg matrix. Assumes that the Hessenberg matrix has already been constructed (by 
     * calling the operator() function).
     * 
     * TODO implement in parent class eventually.
     * 
     * Parameters
     * ----------
     * Eigen::MatrixXcd& S
     *  Schur matrix (upper triangular) similar to original Rayleigh quotient.
     */
    void compute_eigensystem(Eigen::MatrixXcd& S)
    {

      std::cout << GridLogMessage << "Computing eigenvalues." << std::endl;

      evecs.clear();

      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es;
      es.compute(S);
      evals = es.eigenvalues();
      littleEvecs = es.eigenvectors();

      // Convert evecs to lattice fields
      for (int k = 0; k < evals.size(); k++) {
        Eigen::VectorXcd vec = littleEvecs.col(k);
        Field tmp (basis[0].Grid());
        tmp = Zero();
        for (int j = 0; j < basis.size(); j++) {
          tmp = tmp + vec[j] * basis[j];
        }
        evecs.push_back(tmp);
      }

      std::cout << GridLogMessage << "Eigenvalues: " << std::endl << evals << std::endl;

    }

    /**
     * Verifies the factorization DV = V^\dag H + f e^\dag with the last-computed 
     * V, H, f. 
     */
    // RealD verifyFactorization() {
    //   int k = basis.size();         // number of basis vectors, also the size of H.
    //   std::vector<Field> factorized (k, Zero());
    //   Field tmp (FGrid); tmp = Zero();
    //   for (int i = 0; i < basis.size(); i++) {
    //     Linop.Op(basis[i], tmp);
    //   }
    //   // basisRotate(basis, Q, 0, Nk, 0, Nk, Nm);
    //   // Linop.Op(, )
    // }

    /* Getters */
    Eigen::MatrixXcd    getHessenbergMat()  { return Hess; }
    Field               getF()              { return f; }
    std::vector<Field>  getBasis()          { return basis; }
    Eigen::VectorXcd    getEvals()          { return evals; }
    std::vector<Field>  getEvecs()          { return evecs; }

    /**
     * Implements implicit restarting for Arnoldi. Assumes eigenvalues are sorted. 
     * 
     * Parameters
     * ----------
     * int _Nm
     *  Size of basis to keep (Hessenberg is MxM).
     * int Nk
     *  Number of basis vectors to keep at each restart.
     */
    void implicitRestart(int _Nm, int _Nk) {
      assert ( _Nk <= _Nm );
      Nm = _Nm; Nk = _Nk;
      int Np = Nm - Nk;       // keep Nk smallest (or largest, depends on sort function) evecs
      
      std::cout << GridLogMessage << "Computing QR Factorizations." << std::endl;

      Eigen::MatrixXcd Q = Eigen::MatrixXcd::Identity(Nm, Nm);
      Eigen::MatrixXcd Qi (Nm, Nm);
      Eigen::MatrixXcd R (Nm, Nm);

      for (int i = Nk; i < Nm; i++) {        // keep the first Nk eigenvalues and iterate through the last Np. Should loop Np times

        // Useful debugging output
        std::cout << GridLogDebug << "Computing QR factorization for i = " << i << std::endl;
        std::cout << GridLogDebug << "Eval shift = " << evals[i] << std::endl;
        std::cout << GridLogDebug << "Hess before rotation: " << Hess << std::endl;

        // QR factorize 
        Eigen::HouseholderQR<Eigen::MatrixXcd> QR (Hess - evals[i] * Eigen::MatrixXcd::Identity(Nm, Nm));
        Qi = QR.householderQ();
        Q = Q * Qi;
        Hess = Qi.adjoint() * Hess * Qi;

        std::cout << GridLogDebug << "Qt up to i = " << Q.transpose() << std::endl;

      }

      std::cout << GridLogDebug << "Hess after all rotations: " << std::endl << Hess << std::endl; 

      // form Arnoldi vector f: f is normal to the basis vectors and its norm \beta is used to determine the Ritz estimate. 
      std::complex<double> beta = Hess(Nk, Nk-1);
      std::complex<double> sigma = Q(Nm-1, Nk-1);
      f = basis[Nk] * beta + f * sigma;
      RealD betak = std::sqrt(norm2(f));
      std::cout << GridLogMessage << "|f|^2 after implicit restart = " << norm2(f) << std::endl;

      // Rotate basis by Qt
      Qt = Q.transpose();
      basisRotate(basis, Qt, 0, Nk + 1, 0, Nm, Nm);

      // rotate
      basisRotate(evecs, Qt, 0, Nk + 1, 0, Nm, Nm);

      // Truncate the basis and restart
      basis = std::vector<Field> (basis.begin(), basis.begin() + Nk);
      // evecs = std::vector<Field> (evecs.begin(), evecs.begin() + Nk);
      Hess = Hess(Eigen::seqN(0, Nk), Eigen::seqN(0, Nk));

      std::cout << "evecs size: " << evecs.size() << std::endl;

    }
  
    /**
     * Computes the number of Arnoldi eigenvectors that have converged. An eigenvector s is considered converged 
     * for a tolerance epsilon if 
     *    r(s) := |\beta e_m^T s| < epsilon
     * where beta is the norm of f_{m+1}.
     * 
     * Parameters
     * ----------
     * 
     * Returns
     * -------
     * int : Number of converged eigenvectors.
     */
    int converged() {
      int Nconv = 0;
      for (int k = 0; k < evecs.size(); k++) {
        RealD emTs = std::abs(littleEvecs(Nm - 1, k));           // e_m^T s
        RealD ritzEstimate = beta_k * emTs;
        // TODO should be ritzEstimate < Tolerance * lambda_max
        std::cout << GridLogMessage << "Ritz estimate for evec " << k << " = " << ritzEstimate << std::endl;
        if (ritzEstimate < rtol) {
          Nconv++;
        }
      }
      return Nconv;
    }

};
    
NAMESPACE_END(Grid);
#endif
