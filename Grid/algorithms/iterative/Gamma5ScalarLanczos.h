/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    γ5-Scalar Lanczos algorithm for γ5-Hermitian operators (Wilson Dirac).

    Single-vector Lanczos using the γ5-inner product (u,v)_{γ5} = u†γ5v.
    Builds the Krylov space {q, D_W q, D_W² q, ...} and projects D_W onto it,
    producing a real tridiagonal matrix T_m whose eigenvalues approximate those
    of D_W.  Complex conjugate pairs of D_W eigenvalues appear as complex
    eigenvalues of the real non-symmetric T_m.

    Recurrence (G_k = sign(q_k†γ5 q_k) = ±1):
        r_k      = D_W q_k − α_k q_k − β_sup[k-1] q_{k-1}
        β_sub[k] = √|r_k†γ5 r_k|           (lower subdiag T_{k+1,k})
        β_sup[k] = β_sub[k] G_{k+1}/G_k    (upper superdiag T_{k,k+1})
        q_{k+1}  = r_k / β_sub[k]

    T_m is real but not symmetric in general; eigenvalues via EigenSolver.

*************************************************************************************/
#ifndef GRID_GAMMA5_SCALAR_LANCZOS_H
#define GRID_GAMMA5_SCALAR_LANCZOS_H

#include <functional>
#include <numeric>
#include <iomanip>

NAMESPACE_BEGIN(Grid);

template<class Field>
class Gamma5ScalarLanczos {
public:
  using Gamma5Func = std::function<void(const Field&, Field&)>;

private:
  typedef Eigen::MatrixXcd CMat;
  typedef Eigen::VectorXcd CVec;
  typedef Eigen::MatrixXd  RMat;

  LinearOperatorBase<Field>& Linop;
  GridBase*                  Grid_;
  Gamma5Func                 applyGamma5;
  RealD                      Tolerance;
  int                        nSteps;

  std::vector<Field>  basis;
  std::vector<RealD>  alpha_;    // diagonal of T_m
  std::vector<RealD>  beta_sub_; // lower subdiagonal: T_{k+1,k}
  std::vector<RealD>  beta_sup_; // upper superdiagonal: T_{k,k+1}
  std::vector<int>    gsign_;    // G_k = ±1

  CVec               evals_;
  std::vector<Field> evecs_;
  std::vector<RealD> residuals_;

public:
  bool doEvalCheck = false;

  Gamma5ScalarLanczos(LinearOperatorBase<Field>& op, GridBase* grid,
                      Gamma5Func g5, RealD tol = 1e-8)
    : Linop(op), Grid_(grid), applyGamma5(g5), Tolerance(tol), nSteps(0)
  {}

  CVec               getEvals()     { return evals_;     }
  std::vector<Field> getEvecs()     { return evecs_;     }
  std::vector<RealD> getResiduals() { return residuals_; }

  // v0       : starting vector
  // maxSteps : maximum Lanczos steps
  // Nstop    : target Ritz pairs to return (0 = all)
  // reorthog : full γ5-reorthogonalisation at each step
  // filter   : eigenvalue selection criterion
  void operator()(const Field& v0, int maxSteps, int Nstop,
                  bool reorthog = false, RitzFilter filter = EvalImNormSmall)
  {
    basis.clear(); alpha_.clear(); beta_sub_.clear(); beta_sup_.clear(); gsign_.clear();
    nSteps = 0;

    // Normalise starting vector in γ5-norm: q_0 = v0 / √|v0†γ5 v0|
    Field q(Grid_); q = v0;
    Field g5q(Grid_);
    applyGamma5(q, g5q);
    RealD omega = std::real(toStdCmplx(innerProduct(g5q, q)));
    RealD beta0 = std::sqrt(std::abs(omega));
    assert(beta0 > 1e-14 && "starting vector has zero γ5-norm");
    q *= (1.0 / beta0);
    gsign_.push_back(omega >= 0 ? +1 : -1);
    basis.push_back(q);

    std::cout << GridLogMessage
              << "Gamma5ScalarLanczos: start  (v0†γ5 v0)=" << omega
              << "  G[0]=" << gsign_[0] << std::endl;

    for (int k = 0; k < maxSteps; k++) {
      // ── Apply operator ─────────────────────────────────────────────────
      Field p(Grid_);
      Linop.Op(basis[k], p);

      // ── α_k = (q_k†γ5 D_W q_k) / G_k  (real) ─────────────────────────
      applyGamma5(basis[k], g5q);
      RealD a = std::real(toStdCmplx(innerProduct(g5q, p))) / RealD(gsign_[k]);
      alpha_.push_back(a);

      // ── r = D_W q_k − α_k q_k − β_sup[k-1] q_{k-1} ───────────────────
      Field r(Grid_); r = p;
      r -= basis[k] * a;
      if (k > 0)
        r -= basis[k-1] * beta_sup_[k-1];

      // ── Optional γ5-reorthogonalisation ────────────────────────────────
      if (reorthog) {
        for (int i = 0; i <= k; i++) {
          applyGamma5(basis[i], g5q);
          ComplexD c = toStdCmplx(innerProduct(g5q, r)) / RealD(gsign_[i]);
          r -= basis[i] * c;
        }
      }

      // ── β_sub[k] = √|r†γ5 r| ───────────────────────────────────────────
      applyGamma5(r, g5q);
      omega = std::real(toStdCmplx(innerProduct(g5q, r)));
      RealD beta = std::sqrt(std::abs(omega));

      std::cout << GridLogMessage
                << "Gamma5ScalarLanczos: k=" << k
                << "  alpha=" << a
                << "  beta=" << beta
                << "  G[k]=" << gsign_[k]
                << "  (r†γ5r)=" << omega << std::endl;

      beta_sub_.push_back(beta);

      if (beta < Tolerance) {
        std::cout << GridLogMessage
                  << "Gamma5ScalarLanczos: invariant subspace at step " << k << std::endl;
        nSteps = k + 1;
        break;
      }

      int G_next = (omega >= 0) ? +1 : -1;
      gsign_.push_back(G_next);
      // β_sup[k] = T_{k,k+1} = β_sub[k] * G_{k+1} / G_k
      beta_sup_.push_back(beta * RealD(G_next) / RealD(gsign_[k]));

      Field qnew(Grid_); qnew = r; qnew *= (1.0 / beta);
      basis.push_back(qnew);
      nSteps = k + 1;
    }

    if (nSteps == 0) return;
    computeRitzPairs(nSteps, Nstop, filter);
  }

private:
  void computeRitzPairs(int m, int Nstop, RitzFilter filter)
  {
    // Assemble real tridiagonal T_m
    RMat Tm = RMat::Zero(m, m);
    for (int k = 0; k < m; k++) {
      Tm(k, k) = alpha_[k];
      if (k + 1 < m) {
        Tm(k+1, k) = beta_sub_[k];
        Tm(k,   k+1) = beta_sup_[k];
      }
    }

    std::cout << GridLogMessage
              << "Gamma5ScalarLanczos: T_m (" << m << "×" << m << "):" << std::endl;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++)
        std::cout << "  " << std::setprecision(8) << std::setw(16) << Tm(i,j);
      std::cout << std::endl;
    }

    Eigen::EigenSolver<RMat> es(Tm);
    CVec lambdas = es.eigenvalues();
    CMat Y       = es.eigenvectors();

    // Sort by filter criterion
    ComplexComparator cComp(filter);
    std::vector<int> idx(m);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      return cComp(toStdCmplx(lambdas(a)), toStdCmplx(lambdas(b)));
    });

    int Nout = (Nstop > 0 && Nstop < m) ? Nstop : m;
    evals_.resize(Nout);
    evecs_.clear();   evecs_.reserve(Nout);
    residuals_.clear(); residuals_.reserve(Nout);

    // Trailing β for residual estimate: β_sub[m-1] * |y_j[m-1]|
    RealD beta_trail = (m - 1 < (int)beta_sub_.size()) ? beta_sub_[m-1] : 0.0;

    for (int ji = 0; ji < Nout; ji++) {
      int  j  = idx[ji];
      evals_(ji) = lambdas(j);
      CVec yj = Y.col(j);

      // Ritz vector: u_j = Σ_k q_k y_j[k]
      Field uj(Grid_); uj = Zero();
      for (int k = 0; k < m && k < (int)basis.size(); k++)
        uj += basis[k] * toStdCmplx(yj(k));
      evecs_.push_back(uj);

      RealD res = beta_trail * std::abs(toStdCmplx(yj(m-1)));
      if (!std::isfinite(res)) res = std::numeric_limits<RealD>::infinity();
      residuals_.push_back(res);

      std::cout << GridLogMessage
                << "  Ritz[" << std::setw(3) << ji << "]"
                << "  lambda=" << evals_(ji)
                << "  |res|=" << res;

      if (doEvalCheck) {
        Field w(Grid_);
        Linop.Op(uj, w);
        w -= uj * evals_(ji);
        RealD actual = std::sqrt(norm2(w));
        std::cout << "  ||D_W u - λu||=" << actual;
      }
      std::cout << std::endl;
    }
  }
};

NAMESPACE_END(Grid);

#endif // GRID_GAMMA5_SCALAR_LANCZOS_H
