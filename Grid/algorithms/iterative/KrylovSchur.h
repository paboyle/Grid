/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./lib/algorithms/iterative/KrylovSchur.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Chulwoo Jung <chulwoo@bnl.gov>
Author: Christoph Lehner <clehner@bnl.gov>
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
#ifndef GRID_KRYLOVSCHUR_H
#define GRID_KRYLOVSCHUR_H

NAMESPACE_BEGIN(Grid); 

/**
 * Computes a complex Schur decomposition of a complex matrix A using Eigen's matrix library. The Schur decomposition,
 *    A = Q^\dag S Q
 * factorizes A into a unitary matrix Q and an upper triangular matrix S. The eigenvalues of A lie on the diagonal of the upper triangular matrix S. 
 * The Schur decomposition is not unique: in particular, any ordering of the eigenvalues of A can be used as the diagonal of the matrix S. 
 * This class supports eigenvalue reordering by swapping two adjacent eigenvalues with a unitary transformation. 
 */
class ComplexSchurDecomposition {

  private:

    typedef Eigen::MatrixXcd CMat;

    CMat A;           // Matrix to decompose, A = Q^\dag S Q
    CMat Q;           // Unitary matrix Q
    CMat S;           // Upper triangular matrix S

    // Placeholders for Givens rotation
    CMat Givens;      // Givens rotation
    ComplexD s;       // off-diagonal element
    ComplexD lam1;    // First eval for swap
    ComplexD lam2;    // Second eval for swap
    ComplexD phi;     // phase of s
    RealD r;          // norm of s and lam2 - lam1
    
    int Nm;                               // size of matrix problem

    ComplexComparator cCompare;           // function to sort the Schur matrix. 
  
  public:

    /**
     * If the input matrix _A is in Hessenberg form (upper triangular + first subdiagonal non-zero), then the Schur 
     * decomposition is easier to compute. 
     */
    ComplexSchurDecomposition(CMat _A, bool isHess, RitzFilter ritzFilter = EvalReSmall) : A(_A), Nm (_A.rows()), cCompare (ritzFilter)
    {
      Eigen::ComplexSchur<CMat> schur (Nm);
      if (isHess) {
        schur.computeFromHessenberg(_A, CMat::Identity(Nm, Nm), true);
      } else {
        schur.compute(_A, true);
      }
      S = schur.matrixT();
      Q = schur.matrixU().adjoint();      // Eigen computes A = Q S Q^\dag, we want A = Q^\dag S Q
    }

    // Getters
    int  getNm()      { return Nm; }         // size of matrix problem
    CMat getMatrixA() { return A;  }         // matrix for decomposition
    CMat getMatrixQ() { return Q;  }         // unitary matrix Q
    CMat getMatrixS() { return S;  }         // Schur matrix (upper triangular) S
    CMat getRitz()    { return S.diagonal(); }

    /**
     * Checks the Schur decomposition A = Q^\dag S Q holds for the computed matrices. Returns if the relative 
     * Frobenius norm || A - Q^\dag S Q || / || A || is less than rtol. 
     */
    bool checkDecomposition(RealD rtol = 1e-8) {
      RealD Anorm = A.norm();
      if (Anorm < rtol) {
        std::cout << GridLogMessage << "Zero matrix" << std::endl;
        return true;
      }

      std::cout << GridLogDebug << "S = " << std::endl << S << std::endl;
      std::cout << GridLogDebug << "Q = " << std::endl << Q << std::endl;

      CMat A2 = Q.adjoint() * S * Q;
      std::cout << GridLogDebug << "Q^dag S Q = " << std::endl << A2 << std::endl;
      RealD dA = (A - A2).norm() / Anorm;
      return (dA < rtol);
    }

    /**
     * Swaps the components on the diagonal of the Schur matrix at index i with index i + 1. 
     * Updates the orthogonal matrix Q accordingly. 
     */
    void swapEvals(int i) {

      assert(0 <= i && i <= Nm - 1);        // can only swap blocks with upper left index between 0 and Nm - 1

      // get parameters for rotation
      s    = S(i, i+1);
      lam1 = S(i, i);
      lam2 = S(i+1, i+1);
      phi  = s / std::abs(s);
      r    = std::sqrt(std::pow(std::abs(s), 2) + std::pow(std::abs(lam2 - lam1), 2));

      // compute Givens rotation corresponding to these parameters
      Givens = CMat::Identity(Nm, Nm);
      Givens(i, i)     = std::abs(s) / r;
      Givens(i+1, i+1) = Givens(i, i);
      Givens(i, i+1)   = (phi / r) * std::conj(lam2 - lam1);
      Givens(i+1, i)   = -std::conj(Givens(i, i+1));
      
      // rotate Schur matrix and unitary change of basis matrix Q
      S = Givens * S * Givens.adjoint();
      Q = Givens * Q;

      return;
    }

    /**
     * Reorders a Schur matrix &Schur to have the Ritz values that we would like to keep for 
     * restart as the first Nk elements on the diagonal. 
     * 
     * This algorithm is implemented as Nk iterations of a a reverse bubble sort with comparator compare.
     * TODO: pass in compare function as an argument, default to compare with <. 
     */
    // void schurReorder(int Nk, std::function compare) {
    void schurReorder(int Nk) {
      for (int i = 0; i < Nk; i++) {
        for (int k = 0; k <= Nm - 2; k++) {
          int idx = Nm - 2 - k;
          // TODO use RitzFilter enum here
          // if (std::abs(S(idx, idx)) < std::abs(S(idx+1, idx+1))) {            // sort by largest modulus of eigenvalue
          // if (std::real(S(idx, idx)) > std::real(S(idx+1, idx+1))) {        // sort by smallest real eigenvalue
          if ( cCompare(S(idx+1, idx+1), S(idx, idx)) ) {            // sort by largest modulus of eigenvalue
            swapEvals(idx);
          }
        }
      }
      return;
    }

    void schurReorderBlock() {
      // TODO method stub
      return;
    }

};

/**
 * Implementation of the Krylov-Schur algorithm.
 */
template<class Field> 
class KrylovSchur {

  private:
  
    std::string cname = std::string("KrylovSchur");
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
    Field u;                                // Residual vector perpendicular to Krylov space (u_{k+1} in notes)
    Eigen::VectorXcd b;                     // b vector in Schur decomposition (e_{k+1} in Arnoldi).
    std::vector<Field> basis;               // orthonormal Krylov basis
    Eigen::MatrixXcd Rayleigh;              // Rayleigh quotient of size Nbasis (after construction)
    Eigen::MatrixXcd Qt;                    // Transpose of basis rotation which projects out high modes.

    Eigen::VectorXcd evals;                 // evals of Rayleigh quotient
    Eigen::MatrixXcd littleEvecs;           // Nm x Nm evecs matrix
    std::vector<Field> evecs;               // Vector of evec fields

    RitzFilter ritzFilter;                        // how to sort evals

  public:       

    KrylovSchur(LinearOperatorBase<Field> &_Linop, GridBase *_Grid, RealD _Tolerance, RitzFilter filter = EvalReSmall)
      : Linop(_Linop), Grid(_Grid), Tolerance(_Tolerance), ritzFilter(filter), u(_Grid), MaxIter(-1), Nm(-1), Nk(-1), Nstop (-1),
        evals (0), evecs (), ssq (0.0), rtol (0.0), beta_k (0.0), approxLambdaMax (0.0)
    {
      u = Zero();
    };


    /* Getters */
    Eigen::MatrixXcd    getRayleighQuotient()  { return Rayleigh; }
    Field               getU()                 { return u;        }
    std::vector<Field>  getBasis()             { return basis;    }
    Eigen::VectorXcd    getEvals()             { return evals;    }
    std::vector<Field>  getEvecs()             { return evecs;    }

    /**
     * Runs the Krylov-Schur loop.
     *   - Runs an Arnoldi step to generate the Rayleigh quotient and Krylov basis.
     *   - Schur decompose the Rayleigh quotient. 
     *   - Permutes the Rayleigh quotient according to the eigenvalues. 
     *   - Truncate the Krylov-Schur expansion. 
     */
    void operator()(const Field& v0, int _maxIter, int _Nm, int _Nk, int _Nstop, bool doubleOrthog = true) {
      MaxIter = _maxIter;
      Nm = _Nm; Nk = _Nk;
      Nstop = _Nstop;
      
      ssq = norm2(v0);
      RealD approxLambdaMax = approxMaxEval(v0);
      rtol = Tolerance * approxLambdaMax;
      // rtol = Tolerance;

      b = Eigen::VectorXcd::Zero(Nm);       // start as e_{k+1}
      b(Nm-1) = 1.0;

      int start = 0;
      Field startVec = v0;
      littleEvecs = Eigen::MatrixXcd::Zero(Nm, Nm);
      for (int i = 0; i < MaxIter; i++) {
        std::cout << GridLogMessage << "Restart Iteration " << i << std::endl;

        // Perform Arnoldi steps to compute Krylov basis and Rayleigh quotient (Hess)
        arnoldiIteration(startVec, Nm, start, doubleOrthog);
        startVec = u;        // original code
        start = Nk;

        // checkKSDecomposition();

        // Perform a Schur decomposition on Rayleigh
        // ComplexSchurDecomposition schur (Rayleigh, false);
        ComplexSchurDecomposition schur (Rayleigh, false, ritzFilter);
        std::cout << GridLogDebug << "Schur decomp holds? " << schur.checkDecomposition() << std::endl;

        // Rearrange Schur matrix so wanted evals are on top left (like MATLAB's ordschur)
        std::cout << GridLogMessage << "Reordering Schur eigenvalues" << std::endl;
        schur.schurReorder(Nk);
        Eigen::MatrixXcd Q = schur.getMatrixQ();
        Qt = Q.adjoint();                           // TODO should Q be real?
        Eigen::MatrixXcd S = schur.getMatrixS();
        // std::cout << GridLogDebug << "Schur decomp holds after reorder? " << schur.checkDecomposition() << std::endl;

        std::cout << GridLogMessage << "*** ROTATING TO SCHUR BASIS *** " << std::endl;

        // Rotate Krylov basis, b vector, redefine Rayleigh quotient and evecs, and truncate.
        Rayleigh = schur.getMatrixS();
        b = Q * b;            // b^\dag = b^\dag * Q^\dag <==> b = Q*b
        // basisRotate(basis, Q, 0, Nm, 0, Nm, Nm);
        // basisRotate(evecs, Q, 0, Nm, 0, Nm, Nm);

        std::vector<Field> basis2; 
        constructUR(basis2, basis, Qt, Nm);
        basis = basis2;

        // std::vector<Field> evecs2;
        // constructUR(evecs2, evecs, Qt, Nm);
        // constructRU(evecs2, evecs, Q, Nm);
        // evecs = evecs2;
        // littleEvecs = littleEvecs * Q.adjoint();      // TODO try this and see if it works
        // littleEvecs = Q * littleEvecs;      // TODO try this and see if it works
        // std::cout << GridLogDebug << "Ritz vectors rotated correctly? " << checkEvecRotation() << std::endl;

        // checkKSDecomposition();

        std::cout << GridLogMessage << "*** TRUNCATING FOR RESTART *** " << std::endl;

        Rayleigh = Rayleigh(Eigen::seqN(0, Nk), Eigen::seqN(0, Nk));
        basis = std::vector<Field> (basis.begin(), basis.begin() + Nk);
        // evecs = std::vector<Field> (evecs.begin(), evecs.begin() + Nk);
        // littleEvecs = littleEvecs(Eigen::seqN(0, Nk), Eigen::seqN(0, Nk));
        Eigen::VectorXcd btmp = b.head(Nk);             // careful about how you slice Eigen::VectorXcd!
        b = btmp;

        checkKSDecomposition();

        // Compute eigensystem of Rayleigh. Note the eigenvectors correspond to the sorted eigenvalues.
        std::cout << GridLogDebug << "S: " << std::endl << S << std::endl;
        // computeEigensystem(S);
        computeEigensystem(Rayleigh);
        std::cout << GridLogMessage << "Eigenvalues (first Nk sorted): " << std::endl << evals << std::endl;

        // check convergence and return if needed. 
        int Nconv = converged();
        std::cout << GridLogMessage << "Number of evecs converged: " << Nconv << std::endl;
        if (Nconv >= Nstop || i == MaxIter - 1) {
          std::cout << GridLogMessage << "Converged with " << Nconv << " / " << Nstop << " eigenvectors on iteration " 
                        << i << "." << std::endl;
          basisRotate(evecs, Qt, 0, Nk, 0, Nk, Nm);
          std::cout << GridLogMessage << "Eigenvalues: " << evals << std::endl;
          return;
        }
      }
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
    void arnoldiIteration(const Field& v0, int Nm, int start = 0, bool doubleOrthog = true)
    {

      ComplexD coeff;
      Field w (Grid);           // A acting on last Krylov vector. 

      if (start == 0) {       // initialize everything that we need.
        RealD v0Norm = 1 / std::sqrt(ssq);
        basis.push_back(v0Norm * v0);                // normalized source

        Rayleigh = Eigen::MatrixXcd::Zero(Nm, Nm);
        u = Zero();
      } else {
        assert( start == basis.size() );      // should be starting at the end of basis (start = Nk)
        std::cout << GridLogMessage << "Resetting Rayleigh and b" << std::endl;
        Eigen::MatrixXcd RayleighCp = Rayleigh;
        Rayleigh = Eigen::MatrixXcd::Zero(Nm, Nm);
        Rayleigh(Eigen::seqN(0, Nk), Eigen::seqN(0, Nk)) = RayleighCp;

        // append b^\dag to Rayleigh, add u to basis
        Rayleigh(Nk, Eigen::seqN(0, Nk)) = b.adjoint();
        basis.push_back(u);
        b = Eigen::VectorXcd::Zero(Nm);
      }

      // Construct next Arnoldi vector by normalizing w_i = Dv_i - \sum_j v_j h_{ji}
      for (int i = start; i < Nm; i++) {
        Linop.Op(basis.back(), w);
        for (int j = 0; j < basis.size(); j++) {
          coeff = innerProduct(basis[j], w);       // coeff = h_{ij}. Note that since {vi} is ONB it's OK to subtract it off after. 
          Rayleigh(j, i) = coeff;
          w -= coeff * basis[j];
        }

        if (doubleOrthog) {
          for (int j = 0; j < basis.size(); j++) {
            coeff = innerProduct(basis[j], w);      // see if there is any residual component in basis[j] direction
            Rayleigh(j, i) += coeff;                // if coeff is non-zero, adjust Rayleigh
            w -= coeff * basis[j];
          }
        }

        // add w_i to the pile
        if (i < Nm - 1) {
          coeff = std::sqrt(norm2(w));
          Rayleigh(i+1, i) = coeff;
          basis.push_back(
            (1.0/coeff) * w
          );
        }

        // after iterations, update u and beta_k = ||u|| before norm
        u = w;                                // make sure u is normalized
        beta_k = std::sqrt(norm2(u));         // beta_k = ||f_k|| determines convergence.
        u = (1/beta_k) * u;
      }

      b(Nm - 1) = beta_k;

      // std::cout << GridLogMessage << "|f|^2 after Arnoldi step = " << norm2(f) << std::endl;
      std::cout << GridLogMessage << "beta_k = |u| (before norm) after Arnoldi step = " << beta_k << std::endl;
      std::cout << GridLogDebug << "Computed Rayleigh quotient = " << std::endl << Rayleigh << std::endl;

      return;
    }

    /**
     * Approximates the eigensystem of the linear operator by computing the eigensystem of 
     * the Rayleigh quotient. Assumes that the Rayleigh quotient has already been constructed (by 
     * calling the operator() function).
     * 
     * TODO implement in parent class eventually.
     * 
     * Parameters
     * ----------
     * Eigen::MatrixXcd& S
     *  Schur matrix (upper triangular) similar to original Rayleigh quotient.
     */
    void computeEigensystem(Eigen::MatrixXcd& S)
    {
      std::cout << GridLogMessage << "Computing eigenvalues." << std::endl;

      evecs.clear();
      evals = S.diagonal();

      // TODO: is there a faster way to get the eigenvectors of a triangular matrix?
      // Rayleigh.triangularView<Eigen::Upper> tri;

      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es;
      // es.compute(Rayleigh);
      es.compute(S);
      littleEvecs = es.eigenvectors();

      std::cout << GridLogDebug << "Little evecs: " << littleEvecs << std::endl;

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
    }

    /**
     * Approximates the maximum eigenvalue of Linop.Op to normalize the residual and test for convergence. 
     * 
     * TODO implement in parent class eventually
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
     * Computes the number of Krylov-Schur eigenvectors that have converged. An eigenvector s is considered converged 
     * for a tolerance epsilon if 
     *    r(s) := |\beta e_m^T s| < epsilon
     * where beta is the norm of f_{m+1}.
     * 
     * TODO implement in parent class eventually
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
      std::cout << GridLogDebug << "b: " << b << std::endl;
      Field tmp (Grid); Field fullEvec (Grid);
      for (int k = 0; k < evecs.size(); k++) {
        Eigen::VectorXcd evec_k = littleEvecs.col(k);
        RealD ritzEstimate = std::abs(b.dot(evec_k));           // b^\dagger s
        std::cout << GridLogMessage << "Ritz estimate for evec " << k << " = " << ritzEstimate << std::endl;
        if (ritzEstimate < rtol) {
          Nconv++;
        }

      }
      // Check that Ritz estimate is explicitly || D (Uy) - lambda (Uy) ||
      // checkRitzEstimate();
      return Nconv;
    }

    /**
     * Checks the Krylov-Schur decomposition DU = UR + f b^\dag with the last-computed 
     * U, R, f, and b. 
     */
    bool checkKSDecomposition(RealD tol = 1e-8) {

      std::cout << GridLogMessage << "*** CHECKING KRYLOV-SCHUR DECOMPOSITION *** " << std::endl;

      int k = basis.size();         // number of basis vectors, also the size of Rayleigh.

      // rotate basis by Rayleigh to construct UR
      std::vector<Field> rotated = basis;
      Eigen::MatrixXcd Rt = Rayleigh.adjoint();
      constructUR(rotated, basis, Rayleigh, k);             // manually rotate
      // basisRotate(rotated, Rt, 0, k, 0, k, k);           // UR

      // TODO: make a new function that I'm positive does what this is doing
      // just take the basis U = (u1 u2 ... uNm) and form the linear combination UR from R
      
      // For each i, form D u(i) and subtract off (US - u b^\dag)(i)
      RealD delta = 0.0; RealD deltaSum = 0;
      Field tmp (Grid); tmp = Zero();
      for (int i = 0; i < k; i++) {
        Linop.Op(basis[i], tmp);          // tmp = D u(i)
        delta = norm2(tmp - rotated[i] - u * std::conj(b(i)));
        delta = delta / norm2(tmp);               // relative tolerance
        deltaSum += delta;

        std::cout << GridLogDebug << "Iteration " << i << std::endl;
        std::cout << GridLogDebug << "Du = " << norm2(tmp) << std::endl;
        std::cout << GridLogDebug << "rotated = " << norm2(rotated[i]) << std::endl;
        std::cout << GridLogDebug << "b[i] = " << b(i) << std::endl;
        std::cout << GridLogMessage << "Deviation in decomp, column " << i << ": " << delta << std::endl;
      }
      std::cout << GridLogMessage << "Squared sum of relative deviations in decomposition: " << deltaSum << std::endl;

      // std::cout << "[DEBUG] testing basis rotate" << std::endl;
      // std::vector<Field> rotated2;
      // constructUR(rotated2, basis, Rayleigh, k);
      // for (int i = 0; i < k; i++) {
      //   std::cout << "rotated[i] - UR[i] = " << norm2(rotated[i] - rotated2[i]) << std::endl;
      // }

      return deltaSum < tol;
    }

    /**
     * Checks the Ritz vector s was rotated correctly by explicitly recomputing the 
     * eigenvectors of the rotated Rayleigh quotient. 
     */
    bool checkRitzRotation(RealD tol = 1e-8) {
      std::cout << GridLogMessage << "*** CHECKING RITZ VECTOR ROTATION *** " << std::endl;

      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es;
      es.compute(Rayleigh);
      Eigen::MatrixXcd littleEvecs2 = es.eigenvectors();
      RealD dLittle = (littleEvecs2 - littleEvecs).norm() / littleEvecs.norm();
      std::cout << GridLogMessage << "|littleEvecs2 - littleEvecs| / |littleEvecs| = " << dLittle << std::endl;

      std::cout << GridLogMessage << "Forming full eigenvectors" << std::endl;
      RealD delta = 0.0; RealD deltaSum = 0;
      for (int k = 0; k < evals.size(); k++) {
        Eigen::VectorXcd vec = littleEvecs.col(k);
        Field tmpEvec (Grid);
        tmpEvec = Zero();
        for (int j = 0; j < basis.size(); j++) {
          tmpEvec = tmpEvec + vec[j] * basis[j];
        }
        delta = norm2(tmpEvec - evecs[k]) / norm2(evecs[k]);
        std::cout << GridLogMessage << "Deviation in evec " << k << ": " << delta << std::endl;
        deltaSum += delta;
      }
      return deltaSum < tol;
    }

    /**
     * Checks the Ritz estimate R(s) is indeed the deviation of a Ritz eigenvector from being a true eigenvector. 
     */
    void checkRitzEstimate(RealD tol = 1e-8) {
      std::cout << GridLogMessage << "*** CHECKING RITZ ESTIMATE *** " << std::endl;
      Field tmp (Grid);
      for (int k = 0; k < evecs.size(); k++) {
        tmp = Zero();
        Linop.Op(evecs[k], tmp);        // D evec
        RealD ritz = std::sqrt(norm2(tmp - evals[k] * evecs[k]));
        std::cout << "Ritz estimate " << k << " = " << ritz << std::endl;
      }
      return;
    }

    /**
     * Given a vector of fields U (equivalently, a LxN matrix, where L is the number of degrees of 
     * freedom on the lattice field) and an NxN matrix R, forms the product UR. 
     * 
     * Note that I believe this is equivalent to basisRotate(U, R.adjoint(), 0, N, 0, N, N), but I'm 
     * not 100% sure (this will be slower and unoptimized though).
     */
    void constructUR(std::vector<Field>& UR, std::vector<Field> &U, Eigen::MatrixXcd& R, int N) {
      Field tmp (Grid);
      UR.clear();
      for (int i = 0; i < N; i++) {
        tmp = Zero();
        for (int j = 0; j < N; j++) {
          tmp = tmp + U[j] * R(j, i);
        }
        UR.push_back(tmp);
      }
      return;
    }

    /**
     * Same as constructUR but for the product order RU. 
     */
    void constructRU(std::vector<Field>& RU, std::vector<Field> &U, Eigen::MatrixXcd& R, int N) {
      Field tmp (Grid);
      RU.clear();
      for (int i = 0; i < N; i++) {
        tmp = Zero();
        for (int j = 0; j < N; j++) {
          tmp = tmp + R(i, j) * U[j];
        }
        RU.push_back(tmp);
      }
      return;
    }

};
    
NAMESPACE_END(Grid);
#endif
