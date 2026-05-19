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
  // a) Timing: Niter timed calls reporting us/call and GB/s.
  //--------------------------------------------------------------------
#if defined(GRID_CUDA) || defined(GRID_HIP) || defined(GRID_SYCL)
  {
    const int Nwarm = 5;
    const int Niter = 100;

    gaussian(rng, field);

    {
      autoView(v, field, AcceleratorRead);
      for (int i = 0; i < Nwarm; i++) sum_gpu(&v[0], osites);
    }

    RealD t_new;
    {
      autoView(v, field, AcceleratorRead);
      t_new = -usecond();
      for (int i = 0; i < Niter; i++) sum_gpu(&v[0], osites);
      t_new += usecond();
    }

    RealD bytes = (RealD)osites * sizeof(vobj);
    RealD GBs   = bytes / (t_new / Niter) * 1e-3;

    std::cout << GridLogMessage << name << " timing (" << Niter << " calls):" << std::endl;
    std::cout << GridLogMessage
              << "  sum_gpu " << t_new/Niter << " us   " << GBs << " GB/s" << std::endl;
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
  Coordinate mpi  = GridDefaultMpi();

  GridCartesian *UGrid   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), mpi);
  GridCartesian *UGrid_f = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexF::Nsimd()), mpi);

  GridParallelRNG rng(UGrid);
  rng.SeedFixedIntegers({1, 2, 3, 4});
  GridParallelRNG rng_f(UGrid_f);
  rng_f.SeedFixedIntegers({1, 2, 3, 4});

  std::cout << GridLogMessage << "Lattice : " << latt << std::endl;
  std::cout << GridLogMessage << "Volume  : " << UGrid->_gsites << std::endl;

  testReduction<LatticeComplexF>      (UGrid_f, rng_f, "LatticeComplexF",      1      );
  testReduction<LatticeComplexD>      (UGrid,   rng,   "LatticeComplexD",      1      );
  testReduction<LatticeColourMatrixF> (UGrid_f, rng_f, "LatticeColourMatrixF", Nc     );
  testReduction<LatticeColourMatrixD> (UGrid,   rng,   "LatticeColourMatrixD", Nc     );
  testReduction<LatticePropagatorF>   (UGrid_f, rng_f, "LatticePropagatorF",   Ns*Nc  );
  testReduction<LatticePropagatorD>   (UGrid,   rng,   "LatticePropagatorD",   Ns*Nc  );

  std::cout << GridLogMessage << "==============================" << std::endl;
  std::cout << GridLogMessage << passed << " PASSED   " << failed << " FAILED" << std::endl;

  Grid_finalize();
  return (failed > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
