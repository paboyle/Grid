/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Test for Gamma5BlockLanczos on a simple diagonal operator.

    The operator is D = diag(scale_i) where scale is complex random.
    γ5 is taken to be the identity (scalar field has no spin structure),
    so the γ5-inner product reduces to the standard Euclidean inner product
    and the algorithm should find the eigenvalues of D directly.

    For a genuine Wilson Dirac test, pass the actual γ5 functor.

*************************************************************************************/
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

// Diagonal complex operator: out = scale * in
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

  // For LatticeComplex (scalar field) γ5 = identity
  auto gamma5 = [](const Field& in, Field& out){ out = in; };

  Field v0(grid);
  random(RNG, v0);

  const int maxSteps = 20;
  const int Nstop    = 4;
  const RealD tol    = 1e-6;

  std::cout << GridLogMessage
            << "\n========================================" << std::endl;
  std::cout << GridLogMessage
            << " Gamma5BlockLanczos  (maxSteps=" << maxSteps
            << " Nstop=" << Nstop << ")" << std::endl;
  std::cout << GridLogMessage
            << "========================================\n" << std::endl;

  Gamma5BlockLanczos<Field> g5bl(op, grid, gamma5, tol);
  g5bl.doEvalCheck = true;
  g5bl(v0, maxSteps, Nstop, /*reorthog=*/true);

  auto evals = g5bl.getEvals();
  std::cout << GridLogMessage
            << "Gamma5BlockLanczos eigenvalues (" << evals.size() << "):" << std::endl;
  for (int k = 0; k < (int)evals.size(); k++)
    std::cout << GridLogMessage << "  [" << k << "] " << evals(k) << std::endl;

  Grid_finalize();
  return 0;
}
