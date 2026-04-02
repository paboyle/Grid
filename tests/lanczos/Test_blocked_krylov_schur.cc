/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Test for BlockedKrylovSchur: verifies the KS decomposition A V = V H + F B^dag
    by explicit operator applies, before and after each restart.

    Uses DumbOperator (diagonal, real, Hermitian) from Test_synthetic_lanczos.

    Tests Nblock=1 (scalar, regression) and Nblock=2,3 (exercises the B^H fix).

*************************************************************************************/
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

// Diagonal real Hermitian operator (eigenvalues = scale lattice sites)
template<class Field>
class DumbOperator : public LinearOperatorBase<Field> {
public:
  LatticeComplex scale;

  DumbOperator(GridBase* grid) : scale(grid) {
    GridParallelRNG pRNG(grid);
    std::vector<int> seeds({5,6,7,8});
    pRNG.SeedFixedIntegers(seeds);
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
    ComplexD d = innerProduct(in, out);  n1 = real(d);
    d = innerProduct(out, out);          n2 = real(d);
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

  int nFail = 0;

  //--------------------------------------------------------------------
  // Helper lambda: run BKS with doVerify and check it doesn't crash
  //--------------------------------------------------------------------
  auto runTest = [&](const std::string& label, int Nblock, int Nm, int Nk,
                     int maxIter, int Nstop) {
    std::cout << GridLogMessage << "===== " << label << " =====" << std::endl;

    BlockKrylovSchur<Field> bks(op, grid, 1e-6, EvalReSmall);

    std::vector<Field> v0(Nblock, Field(grid));
    for (int t = 0; t < Nblock; t++) random(RNG, v0[t]);

    bks(v0, maxIter, Nm, Nk, Nstop, Nblock,
        /*doubleOrthog=*/true, /*doVerify=*/true);

    std::cout << GridLogMessage << label << " done." << std::endl;
  };

  // Test 1: Nblock=1 — scalar case, regression (Nm,Nk now total vectors)
  runTest("Nblock=1  Nm=10 Nk=5  maxIter=3", 1, 10, 5, 3, 5);

  // Test 2: Nblock=2 — exercises the B^H fix for off-diagonal elements
  runTest("Nblock=2  Nm=16 Nk=8  maxIter=3", 2, 16, 8, 3, 4);

  // Test 3: Nblock=3 — further stress-test the B^H fix
  runTest("Nblock=3  Nm=27 Nk=9  maxIter=3", 3, 27, 9, 3, 3);

  // Test 4: Nblock=2, larger cycle — more restarts
  runTest("Nblock=2  Nm=24 Nk=12  maxIter=5", 2, 24, 12, 5, 6);

  if (nFail == 0)
    std::cout << GridLogMessage << "All BlockedKrylovSchur tests completed." << std::endl;

  Grid_finalize();
  return nFail;
}
