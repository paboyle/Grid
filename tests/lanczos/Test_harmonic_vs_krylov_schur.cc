/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Comparison test: HarmonicBlockedKrylovSchur vs BlockedKrylovSchur.

    Both algorithms are run on the same diagonal Hermitian operator and the
    resulting eigenvalues are compared.  doVerify=true is used so the KS
    decomposition check  max|H-M|  and the per-column residuals are printed
    at each step.  For BKS these should be O(machine epsilon) at all times.
    For HBKS they should be O(machine epsilon) AFTER Arnoldi, but may show
    large per-column deviations AFTER restart+truncation (because the rotation
    Q from Schur(Hhat) does not give an upper-triangular H_new, so the
    truncated KS relation is only approximate).

*************************************************************************************/
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

// Diagonal real Hermitian operator  out = scale * in
template<class Field>
class DumbOperator : public LinearOperatorBase<Field> {
public:
  LatticeComplex scale;
  DumbOperator(GridBase* grid) : scale(grid) {
    GridParallelRNG pRNG(grid);
    pRNG.SeedFixedIntegers({5,6,7,8});
    random(pRNG, scale);
    scale = exp(-Grid::real(scale) * 3.0);
  }
  void OpDirAll(const Field& in, std::vector<Field>& out) {}
  void OpDiag(const Field& in, Field& out) {}
  void OpDir(const Field& in, Field& out, int dir, int disp) {}
  void Op(const Field& in, Field& out) { out = scale * in; }
  void AdjOp(const Field& in, Field& out) { out = scale * in; }
  void HermOp(const Field& in, Field& out) { out = scale * in; }
  void HermOpAndNorm(const Field& in, Field& out, double& n1, double& n2) {
    out = scale * in;
    ComplexD d = innerProduct(in, out); n1 = real(d);
    d = innerProduct(out, out);         n2 = real(d);
  }
};

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  GridCartesian* grid = SpaceTimeGrid::makeFourDimGrid(
    GridDefaultLatt(),
    GridDefaultSimd(Nd, vComplex::Nsimd()),
    GridDefaultMpi());

  GridParallelRNG RNG(grid);
  RNG.SeedFixedIntegers({1,2,3,4});

  typedef LatticeComplex Field;
  DumbOperator<Field> op(grid);

  //----------------------------------------------------------------------
  // Parameters (kept small so output is readable)
  //----------------------------------------------------------------------
  const int Nblock  = 4;
  const int Nm      = 20;   // total vectors (= 5 blocks * Nblock=4)
  const int Nk      = 8;    // total kept   (= 2 blocks * Nblock=4)
  const int Nstop   = 4;
  const int maxIter = 8;
  const RealD tol   = 1e-6;

  // Two identical starting blocks
  std::vector<Field> v0(Nblock, Field(grid));
  std::vector<Field> v0b(Nblock, Field(grid));
  for (int t = 0; t < Nblock; t++) {
    random(RNG, v0[t]);
    v0b[t] = v0[t];   // identical start for fair comparison
  }

  //----------------------------------------------------------------------
  // Run BlockedKrylovSchur with doVerify=true
  //----------------------------------------------------------------------
  std::cout << GridLogMessage
            << "\n========================================" << std::endl;
  std::cout << GridLogMessage
            << " BlockKrylovSchur  (Nblock=" << Nblock
            << " Nm=" << Nm << " Nk=" << Nk << ")" << std::endl;
  std::cout << GridLogMessage
            << "========================================\n" << std::endl;

  BlockKrylovSchur<Field> bks(op, grid, tol, EvalImNormSmall);
  bks(v0, maxIter, Nm, Nk, Nstop, Nblock,
      /*doubleOrthog=*/true, /*doVerify=*/true);

  auto bks_evals = bks.getEvals();
  std::cout << GridLogMessage
            << "BKS eigenvalues (" << bks_evals.size() << "):" << std::endl;
  for (int k = 0; k < (int)bks_evals.size(); k++)
    std::cout << GridLogMessage << "  [" << k << "] " << bks_evals[k] << std::endl;

  //----------------------------------------------------------------------
  // Run HarmonicBlockedKrylovSchur with doVerify=true
  //----------------------------------------------------------------------
  std::cout << GridLogMessage
            << "\n========================================" << std::endl;
  std::cout << GridLogMessage
            << " HarmonicBlockKrylovSchur  (Nblock=" << Nblock
            << " Nm=" << Nm << " Nk=" << Nk << " shift=0)" << std::endl;
  std::cout << GridLogMessage
            << "========================================\n" << std::endl;

  HarmonicBlockKrylovSchur<Field> hbks(op, grid, tol, 0.0, EvalImNormSmall);
  hbks(v0b, maxIter, Nm, Nk, Nstop, Nblock,
       /*doubleOrthog=*/true, /*doVerify=*/true);

  auto hbks_evals = hbks.getEvals();
  std::cout << GridLogMessage
            << "HBKS eigenvalues (" << hbks_evals.size() << "):" << std::endl;
  for (int k = 0; k < (int)hbks_evals.size(); k++)
    std::cout << GridLogMessage << "  [" << k << "] " << hbks_evals[k] << std::endl;

  //----------------------------------------------------------------------
  // Compare
  //----------------------------------------------------------------------
  std::cout << GridLogMessage
            << "\n========================================" << std::endl;
  std::cout << GridLogMessage << " Eigenvalue comparison" << std::endl;
  std::cout << GridLogMessage
            << "========================================" << std::endl;

  // Sort both sets by real part for comparison
  std::vector<ComplexD> bvec(bks_evals.data(),
                             bks_evals.data() + bks_evals.size());
  std::vector<ComplexD> hvec(hbks_evals.data(),
                             hbks_evals.data() + hbks_evals.size());
  auto cmpRe = [](const ComplexD& a, const ComplexD& b){ return a.real() < b.real(); };
  std::sort(bvec.begin(), bvec.end(), cmpRe);
  std::sort(hvec.begin(), hvec.end(), cmpRe);

  int nCmp = std::min(bvec.size(), hvec.size());
  double maxDiff = 0.0;
  for (int k = 0; k < nCmp; k++) {
    double diff = std::abs(bvec[k].real() - hvec[k].real()) + std::abs(bvec[k].imag() - hvec[k].imag());
    maxDiff = std::max(maxDiff, diff);
    std::cout << GridLogMessage
              << "  k=" << k
              << "  BKS=" << bvec[k]
              << "  HBKS=" << hvec[k]
              << "  |diff|=" << diff << std::endl;
  }
  std::cout << GridLogMessage << "  max |BKS - HBKS| = " << maxDiff << std::endl;

  Grid_finalize();
  return 0;
}
