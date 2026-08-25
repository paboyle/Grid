/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_multi_reduction.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */

//////////////////////////////////////////////////////////////////////////////
// Gate for the batched Krylov linear algebra in Lattice_reduction.h:
//   innerProductMulti  vs  innerProduct, one j at a time
//   axpyMulti          vs  sequential axpy
//   axpyMultiNorm      vs  sequential axpy + norm2
// on TWO field types through the SAME code path:
//   fine  : LatticeFermionD on the 4D grid           (vSpinColourVector)
//   coarse: Lattice<iVector<CComplex,nbasis>> on a 5D grid with the extra
//           dimension standing in for nrhs, as the coarse mRHS solver does.
// Batch sizes cover every template width and the >16 chunking path.
//
//   mpirun -n 1 ./Test_multi_reduction --grid 8.8.8.8 --mpi 1.1.1.1
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>

using namespace Grid;

static int failures = 0;
static void Report(const std::string &name, bool pass, const std::string &detail="")
{
  std::cout << GridLogMessage << "  " << name << (pass ? "  PASS" : "  ** FAIL **");
  if ( detail.size() ) std::cout << "   " << detail;
  std::cout << std::endl;
  if ( !pass ) failures++;
}
static double Cabs(const ComplexD &z){ double re=z.real(), im=z.imag(); return std::sqrt(re*re+im*im); }

template<class Field>
void Check(const std::string &tag, GridBase *grid, GridParallelRNG &RNG, const std::vector<int> &batches)
{
  const int Mmax = 20;
  std::vector<Field> x; for(int j=0;j<Mmax;j++){ x.emplace_back(grid); gaussian(RNG,x[j]); }
  Field r(grid); gaussian(RNG,r);
  Field z0(grid); gaussian(RNG,z0);
  Field za(grid), zb(grid), d(grid);

  for(int m : batches){
    std::vector<const Field*> win(m);
    std::vector<ComplexD> b(m);
    for(int j=0;j<m;j++){ win[j]=&x[j]; b[j]=ComplexD(0.3*j-1.0, 0.7-0.2*j); }

    // dots
    std::vector<ComplexD> ip;
    innerProductMulti(ip,win,r);
    double worst=0.0;
    for(int j=0;j<m;j++){
      ComplexD ref = innerProduct(x[j],r);
      worst = std::max(worst, Cabs(ip[j]-ref)/std::max(Cabs(ref),1.0e-300));
    }
    Report(tag+" innerProductMulti m="+std::to_string(m), worst<1.0e-12, "worst rel "+std::to_string(worst));

    // update
    za = z0; axpyMulti(za,b,win);
    zb = z0; for(int j=0;j<m;j++) axpy(zb,b[j],x[j],zb);
    d = za - zb;
    double rel = std::sqrt(norm2(d)/norm2(zb));
    Report(tag+" axpyMulti m="+std::to_string(m), rel<1.0e-12, "rel "+std::to_string(rel));

    // update + norm
    za = z0; RealD nn = axpyMultiNorm(za,b,win);
    RealD nref = norm2(zb);
    d = za - zb;
    rel = std::sqrt(norm2(d)/norm2(zb));
    double nrel = std::fabs(nn-nref)/nref;
    Report(tag+" axpyMultiNorm m="+std::to_string(m), rel<1.0e-12 && nrel<1.0e-12,
           "rel "+std::to_string(rel)+" norm rel "+std::to_string(nrel));
  }
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                             GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  const int nrhs = 3;   // stands in for the coarse mRHS dimension
  GridCartesian *CGrid = SpaceTimeGrid::makeFiveDimGrid(nrhs, UGrid);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers(seeds);
  GridParallelRNG RNG5(CGrid); RNG5.SeedFixedIntegers(seeds);

  std::vector<int> batches({1,2,3,4,5,8,9,16,17,20});

  Check<LatticeFermionD>("fine  ", UGrid, RNG4, batches);

  const int nbasis = 8;
  typedef iVector<iScalar<iScalar<iScalar<vComplexD>>>,nbasis> CoarseSiteVector;
  typedef Lattice<CoarseSiteVector> CoarseField;
  Check<CoarseField>("coarse", CGrid, RNG5, batches);

  std::cout << GridLogMessage << (failures ? "Test_multi_reduction: FAILURES"
                                           : "Test_multi_reduction: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
