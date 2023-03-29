/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/solver/Test_eofa_inv.cc

Copyright (C) 2017

Author: Christopher Kelly <ckelly@bnl.gov>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

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
 ;

int main (int argc, char** argv)
{
  Grid_init(&argc, &argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  const int Ls = 8;

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         *FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  // Want a different conf at every run
  // First create an instance of an engine.
  std::random_device rnd_device;
  // Specify the engine and distribution.
  std::mt19937 mersenne_engine(rnd_device());
  std::uniform_int_distribution<int> dist(1, 100);

  auto gen = std::bind(dist, mersenne_engine);
  std::vector<int> seeds4(4);
  generate(begin(seeds4), end(seeds4), gen);

  //std::vector<int> seeds4({1,2,3,5});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  LatticeFermion phi        (FGrid);  gaussian(RNG5, phi);
  LatticeFermion Mphi       (FGrid);
  LatticeFermion MphiPrime  (FGrid);

  LatticeGaugeField U(UGrid);
  SU<Nc>::HotConfiguration(RNG4,U);

  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD b  = 2.5;
  RealD c  = 1.5;
  RealD mf = 0.01;
  RealD mb = 1.0;
  RealD M5 = 1.8;
  MobiusEOFAFermionD Lop(U, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mf, mf, mb, 0.0, -1, M5, b, c);
  MobiusEOFAFermionD Rop(U, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mb, mf, mb, -1.0, 1, M5, b, c);
  OneFlavourRationalParams Params(0.95, 100.0, 5000, 1.0e-10, 12);
  ConjugateGradient<LatticeFermion> CG(1.0e-10, 5000);
  ExactOneFlavourRatioPseudoFermionAction<WilsonImplR> Meofa(Lop, Rop, CG, CG, CG, CG, CG, Params, false);

  GridSerialRNG  sRNG; sRNG.SeedFixedIntegers(seeds4);


  //Random field
  LatticeFermion eta(FGrid);
  gaussian(RNG5,eta);
  
  //Check left inverse
  LatticeFermion Meta(FGrid);
  Meofa.Meofa(U, eta, Meta);

  LatticeFermion MinvMeta(FGrid);
  Meofa.MeofaInv(U, Meta, MinvMeta);

  LatticeFermion diff = MinvMeta - eta;

  std::cout << GridLogMessage << "eta: " << norm2(eta) << " M*eta: " << norm2(Meta) << " M^{-1}*M*eta: " << norm2(MinvMeta) << "  M^{-1}*M*eta - eta: " << norm2(diff) << " (expect 0)" << std::endl;
  assert(norm2(diff) < 1e-8);

  //Check right inverse
  LatticeFermion MinvEta(FGrid);
  Meofa.MeofaInv(U, eta, MinvEta);

  LatticeFermion MMinvEta(FGrid);
  Meofa.Meofa(U, MinvEta, MMinvEta);

  diff = MMinvEta - eta;
  
  std::cout << GridLogMessage << "eta: " << norm2(eta) << " M^{-1}*eta: " << norm2(MinvEta) << " M*M^{-1}*eta: " << norm2(MMinvEta) << "  M*M^{-1}*eta - eta: " << norm2(diff) << " (expect 0)" << std::endl;
  assert(norm2(diff) < 1e-8);

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
}
