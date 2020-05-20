/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_blocking_operators.cc

    Copyright (C) 2015 - 2020

Author: Daniel Richtmann <daniel.richtmann@ur.de>

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
#include <Benchmark_helpers.h>

using namespace Grid;

#ifndef NBASIS
#define NBASIS 40
#endif

template<class CComplex, int nbasis>
void convertLayout(std::vector<Lattice<CComplex>> const& in, Lattice<iVector<CComplex, nbasis>>& out) {
  assert(in.size() == nbasis);
  for(auto const& elem : in) conformable(elem, out);

  auto  out_v = out.View();
  auto  in_vc = getViewContainer(in);
  auto* in_va = &in_vc[0];
  accelerator_for(ss, out.Grid()->oSites(), CComplex::Nsimd(), {
    for(int i = 0; i < nbasis; i++) { coalescedWrite(out_v[ss](i), in_va[i](ss)); }
  });
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  /////////////////////////////////////////////////////////////////////////////
  //                          Read from command line                         //
  /////////////////////////////////////////////////////////////////////////////

  // clang-format off
  const int  nBasis    = NBASIS; static_assert((nBasis & 0x1) == 0, "");
  const int  nB        = nBasis / 2;
  Coordinate blockSize = Commandline::readCoordinate(&argc, &argv, "--blocksize", Coordinate({4, 4, 4, 4}));
  uint64_t   nIterMin  = Commandline::readInt(&argc, &argv, "--miniter", 1000);
  uint64_t   nSecMin   = Commandline::readInt(&argc, &argv, "--minsec", 5);
  int        Ls        = Commandline::readInt(&argc, &argv, "--Ls", 1);
  // clang-format on

  std::cout << GridLogMessage << "Compiled with nBasis = " << nBasis << " -> nB = " << nB << std::endl;
  std::cout << GridLogMessage << "Using Ls = " << Ls << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  //                              General setup                              //
  /////////////////////////////////////////////////////////////////////////////

  Coordinate clatt = calcCoarseSize(GridDefaultLatt(), blockSize);

  GridCartesian*         UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian*         UGrid_c   = SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid_c = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_c);
  GridCartesian*         FGrid_f   = nullptr;
  GridRedBlackCartesian* FrbGrid_f = nullptr;
  GridCartesian*         FGrid_c   = nullptr;
  GridRedBlackCartesian* FrbGrid_c = nullptr;

  if(Ls != 1) {
    FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid_f);
    FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid_f);
    FGrid_c   = SpaceTimeGrid::makeFiveDimGrid(1, UGrid_c);
    FrbGrid_c = SpaceTimeGrid::makeFiveDimRedBlackGrid(1, UGrid_c);
  } else {
    FGrid_f   = UGrid_f;
    FrbGrid_f = UrbGrid_f;
    FGrid_c   = UGrid_c;
    FrbGrid_c = UrbGrid_c;
  }

  UGrid_f->show_decomposition();
  FGrid_f->show_decomposition();
  UGrid_c->show_decomposition();
  FGrid_c->show_decomposition();

  GridParallelRNG UPRNG_f(UGrid_f);
  GridParallelRNG FPRNG_f(FGrid_f);
  GridParallelRNG UPRNG_c(UGrid_c);
  GridParallelRNG FPRNG_c(FGrid_c);

  std::vector<int> seeds({1, 2, 3, 4});

  UPRNG_f.SeedFixedIntegers(seeds);
  FPRNG_f.SeedFixedIntegers(seeds);
  UPRNG_c.SeedFixedIntegers(seeds);
  FPRNG_c.SeedFixedIntegers(seeds);

  RealD tol = getPrecision<vComplex>::value == 2 ? 1e-15 : 1e-7;

  /////////////////////////////////////////////////////////////////////////////
  //                             Type definitions                            //
  /////////////////////////////////////////////////////////////////////////////

  typedef Aggregation<vSpinColourVector, vTComplex, nBasis>     Aggregates;
  typedef CoarsenedMatrix<vSpinColourVector, vTComplex, nBasis> CoarseDiracMatrix;
  typedef CoarseDiracMatrix::CoarseVector                       CoarseFermionField;
  typedef CoarseDiracMatrix::CoarseMatrix                       CoarseLinkField;
  typedef CoarseDiracMatrix::FineField::vector_object           Fobj;
  typedef Lattice<iScalar<vInteger>>                            Coor;
  typedef Lattice<typename Fobj::tensor_reduced>                FineComplexField;
  typedef typename Fobj::scalar_type                            ScalarType;
  typedef CoarseDiracMatrix::CoarseComplexField                 CoarseComplexField;

  /////////////////////////////////////////////////////////////////////////////
  //                           Setup of Aggregation                          //
  /////////////////////////////////////////////////////////////////////////////

  const int cb = 0;
  const int checkOrthog = 1;
  const int gsPasses = 1;

  Aggregates Aggs(FGrid_c, FGrid_f, cb);
  Aggs.CreateSubspaceRandom(FPRNG_f);
  if(Ls != 1) performChiralDoublingG5R5(Aggs.subspace);
  else        performChiralDoublingG5C(Aggs.subspace);
  Aggs.Orthogonalise(checkOrthog, gsPasses);

  /////////////////////////////////////////////////////////////////////////////
  //             Calculate numbers needed for performance figures            //
  /////////////////////////////////////////////////////////////////////////////

  auto siteElems_f = Ns * Nc;
  auto siteElems_c = nBasis;

  std::cout << GridLogDebug << "siteElems_f = " << siteElems_f << std::endl;
  std::cout << GridLogDebug << "siteElems_c = " << siteElems_c << std::endl;

  double UVolume_f = std::accumulate(UGrid_f->_fdimensions.begin(), UGrid_f->_fdimensions.end(), 1, std::multiplies<double>());
  double FVolume_f = std::accumulate(FGrid_f->_fdimensions.begin(), FGrid_f->_fdimensions.end(), 1, std::multiplies<double>());
  double UVolume_c = std::accumulate(UGrid_c->_fdimensions.begin(), UGrid_c->_fdimensions.end(), 1, std::multiplies<double>());
  double FVolume_c = std::accumulate(FGrid_c->_fdimensions.begin(), FGrid_c->_fdimensions.end(), 1, std::multiplies<double>());

  /////////////////////////////////////////////////////////////////////////////
  //                 Set up stuff required for the benchmark                 //
  /////////////////////////////////////////////////////////////////////////////

  Geometry geom(FGrid_c->_ndimension);

  std::vector<FineComplexField> masks(geom.npoint, FGrid_f);
  std::vector<CoarseningLookupTable> lut(geom.npoint);

  { // original setup code taken from CoarsenedMatrix
    FineComplexField one(FGrid_f); one = ScalarType(1.0, 0.0);
    FineComplexField zero(FGrid_f); zero = ScalarType(0.0, 0.0);
    Coor             coor(FGrid_f);

    for(int p=0; p<geom.npoint; ++p) {
      int     dir   = geom.directions[p];
      int     disp  = geom.displacements[p];
      Integer block = (FGrid_f->_rdimensions[dir] / FGrid_c->_rdimensions[dir]);

      LatticeCoordinate(coor, dir);

      if(disp == 0) {
        masks[p] = Zero();
      } else if(disp == 1) {
        masks[p] = where(mod(coor,block) == (block-1), one, zero);
      } else if(disp == -1) {
        masks[p] = where(mod(coor,block) == (Integer)0, one, zero);
      }

      lut[p].populate(FGrid_c, masks[p]);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                           Start of benchmarks                           //
  /////////////////////////////////////////////////////////////////////////////

  {
    std::cout << GridLogMessage << "***************************************************************************" << std::endl;
    std::cout << GridLogMessage << "Running benchmark for blockMaskedInnerProduct" << std::endl;
    std::cout << GridLogMessage << "***************************************************************************" << std::endl;

    LatticeFermion src(FGrid_f);
    std::vector<CoarseComplexField> res_mask(nBasis, FGrid_c);
    CoarseFermionField res_lut(FGrid_c);

    random(FPRNG_f, src);

    for(int p=0; p<geom.npoint; ++p) {
      res_lut = Zero();
      for(auto& elem : res_mask) elem = Zero();

      auto workSites_f = real(TensorRemove(sum(masks[p]))); // remove sites with zero in its mask from flop and byte counting

      std::cout << GridLogMessage << "point = " << p << ": workSites, volume, ratio = " << workSites_f << ", "
                << FVolume_f << ", " << workSites_f / FVolume_f << std::endl;

      double flop = workSites_f * (8 * siteElems_f) * nBasis;
      double byte = workSites_f * (2 * 1 + 2 * siteElems_f) * nBasis * sizeof(Complex);

      BenchmarkFunctionMRHS(blockMaskedInnerProduct,
                            flop, byte, nIterMin, nSecMin, nBasis,
                            res_mask[rhs], masks[p], Aggs.subspace[rhs], src);

      BenchmarkFunction(blockLutedInnerProduct,
                        flop, byte, nIterMin, nSecMin,
                        res_lut, src, Aggs.subspace, lut[p]);

      if(p != geom.npoint-1) { // does nothing for self stencil point
        CoarseFermionField tmp_lut(FGrid_c);
        convertLayout(res_mask, tmp_lut);
        assertResultMatchesReference(tol, res_lut, tmp_lut);
      }
    }
    std::cout << GridLogMessage << "Results are equal" << std::endl;
  }

  {
    std::cout << GridLogMessage << "***************************************************************************" << std::endl;
    std::cout << GridLogMessage << "Running benchmark for blockProject" << std::endl;
    std::cout << GridLogMessage << "***************************************************************************" << std::endl;

    LatticeFermion src(FGrid_f);
    CoarseFermionField res_project(FGrid_c);
    std::vector<CoarseComplexField> res_mask(nBasis, FGrid_c);
    CoarseFermionField res_lut(FGrid_c);

    random(FPRNG_f, src);

    res_project = Zero();
    for(auto& elem : res_mask) elem = Zero();
    res_lut = Zero();

    FineComplexField      fullMask(FGrid_f); fullMask = ScalarType(1.0, 0.0);
    CoarseningLookupTable fullLut(FGrid_c, fullMask);

    double flop = FVolume_f * (8 * siteElems_f) * nBasis;
    double byte = FVolume_f * (2 * 1 + 2 * siteElems_f) * nBasis * sizeof(Complex);

    BenchmarkFunction(blockProject,
                      flop, byte, nIterMin, nSecMin,
                      res_project, src, Aggs.subspace);

    BenchmarkFunctionMRHS(blockMaskedInnerProduct,
                          flop, byte, nIterMin, nSecMin, nBasis,
                          res_mask[rhs], fullMask, Aggs.subspace[rhs], src);

    BenchmarkFunction(blockLutedInnerProduct,
                      flop, byte, nIterMin, nSecMin,
                      res_lut, src, Aggs.subspace, fullLut);

    assertResultMatchesReference(tol, res_project, res_lut);
    CoarseFermionField tmp_lut(FGrid_c);
    convertLayout(res_mask, tmp_lut);
    assertResultMatchesReference(tol, res_lut, tmp_lut);

    std::cout << GridLogMessage << "Results are equal" << std::endl;
  }

  Grid_finalize();
}
