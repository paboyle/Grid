/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_reduction.cc

    Copyright (C) 2024

Author: Peter Boyle <pboyle@bnl.gov>

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

static int passed = 0;
static int failed = 0;

static void check(bool ok, const std::string &msg)
{
  if (ok) {
    std::cout << GridLogMessage << "PASS  " << msg << std::endl;
    passed++;
  } else {
    std::cout << GridLogMessage << "FAIL  " << msg << std::endl;
    failed++;
  }
}

// Squared magnitude of a Grid scalar tensor aggregate: innerProduct(a,a).
// For iScalar:      real(conj(a)*a)
// For iMatrix<T,N>: sum_{i,j} real(conj(a_ij)*a_ij)  (Frobenius)
// Named squaredSum to make clear the squaring is applied to the aggregate
// (the sum), not to individual site values before summing.
template<class T>
RealD squaredSum(const T &a)
{
  return (RealD)real(TensorRemove(innerProduct(a, a)));
}

template<class Field>
void testReduction(GridCartesian *grid, GridParallelRNG &rng,
                   const std::string &name, int Ncomp)
{
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object  sobj;
  typedef typename vobj::scalar_type    scalar_type;

  const Integer V      = grid->_gsites;
  const Integer osites = grid->oSites();

  // Detect single vs double precision by comparing fundamental scalar sizes.
  const bool isFloat = (sizeof(scalar_type) < sizeof(ComplexD));

  std::cout << GridLogMessage << "=== " << name << " ===" << std::endl;

  Field field(grid);

  //--------------------------------------------------------------------
  // a) Gaussian random field: sum_gpu (new CUB path) vs sum_gpu_old
  //    (preserved hand-rolled shared-memory path).  Both promote lanes
  //    to double internally, so results should agree to near-roundoff.
  //--------------------------------------------------------------------
#if defined(GRID_CUDA) || defined(GRID_HIP) || defined(GRID_SYCL)
  {
    gaussian(rng, field);

    autoView(v, field, AcceleratorRead);
    sobj new_result = sum_gpu    (&v[0], osites);
    sobj old_result = sum_gpu_old(&v[0], osites);

    sobj  diff  = new_result - old_result;
    RealD diffn = squaredSum(diff);
    RealD refn  = squaredSum(old_result);
    RealD reldiff = (refn > 0.0) ? std::sqrt(diffn / refn) : std::sqrt(diffn);

    // Float fields: both paths cast from double to float, expect O(eps_float).
    // Double fields: ordering differences at most O(V * eps_double).
    RealD tol = isFloat ? 1e-6 : 1e-10;

    std::cout << GridLogMessage
              << name << " random reldiff = " << reldiff << std::endl;
    check(reldiff < tol, name + " random: sum_gpu agrees with sum_gpu_old");
  }
#endif

  //--------------------------------------------------------------------
  // b) Constant field via field = 1.0.
  //
  //    Grid's iMatrix::operator=(scalar) sets only the diagonal, so:
  //      LatticeComplex       -> scalar 1.0        (Ncomp = 1  nonzero per site)
  //      LatticeColourMatrix  -> Nc x Nc identity  (Ncomp = Nc nonzero per site)
  //      LatticePropagator    -> (Ns*Nc)^2 identity (Ncomp = Ns*Nc nonzero per site)
  //
  //    After GlobalSum: sum_result has Ncomp diagonal entries each equal to V,
  //    all off-diagonal entries zero.  Grid's recursive innerProduct computes
  //    the Frobenius inner product (sum of |element|^2 over all indices), giving
  //
  //        innerProduct(sum_result, sum_result) = Ncomp * V^2
  //--------------------------------------------------------------------
  {
    field = 1.0;
    sobj sum_result = sum(field);   // uses new GPU path + GlobalSum

    RealD got      = squaredSum(sum_result);
    RealD expected = (RealD)Ncomp * (RealD)V * (RealD)V;
    RealD reldiff  = std::abs(got - expected) / expected;

    std::cout << GridLogMessage
              << name << " const: got " << got
              << "  expected " << expected
              << "  reldiff "  << reldiff << std::endl;
    check(reldiff < 1e-8, name + " const: innerProduct(sum,sum) = Ncomp*V^2");
  }
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  Coordinate latt = GridDefaultLatt();
  Coordinate simd = GridDefaultSimd(Nd, vComplexD::Nsimd());
  Coordinate mpi  = GridDefaultMpi();

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(latt, simd, mpi);
  GridParallelRNG rng(grid);
  rng.SeedFixedIntegers({1, 2, 3, 4});

  std::cout << GridLogMessage << "Lattice : " << latt << std::endl;
  std::cout << GridLogMessage << "Volume  : " << grid->_gsites << std::endl;

  testReduction<LatticeComplexF>      (grid, rng, "LatticeComplexF",      1      );
  testReduction<LatticeComplexD>      (grid, rng, "LatticeComplexD",      1      );
  testReduction<LatticeColourMatrixF> (grid, rng, "LatticeColourMatrixF", Nc     );
  testReduction<LatticeColourMatrixD> (grid, rng, "LatticeColourMatrixD", Nc     );
  testReduction<LatticePropagatorF>   (grid, rng, "LatticePropagatorF",   Ns*Nc  );
  testReduction<LatticePropagatorD>   (grid, rng, "LatticePropagatorD",   Ns*Nc  );

  std::cout << GridLogMessage << "==============================" << std::endl;
  std::cout << GridLogMessage << passed << " PASSED   " << failed << " FAILED" << std::endl;

  Grid_finalize();
  return (failed > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
