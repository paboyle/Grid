/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/debug/Test_heatbath_dwf_eofa.cc

Copyright (C) 2017

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
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

//////////////////////////////////////////////////////////////////////////////////////////
// This program sets up the initial pseudofermion field |Phi> = Meofa^{-1/2}*|eta>, and
// then uses this Phi to compute the action <Phi|Meofa|Phi>.
// If all is working, one should find that <eta|eta> = <Phi|Meofa|Phi>.
//////////////////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
 ;

typedef GparityWilsonImplR FermionImplPolicy;
typedef GparityMobiusEOFAFermionR FermionAction;
typedef typename FermionAction::FermionField FermionField;

// Parameters for test
const std::vector<int> grid_dim = { 8, 8, 8, 8 };
const int              Ls       = 8;
const int              Npoles   = 12;
const RealD            b        = 2.5;
const RealD            c        = 1.5;
const RealD            mf       = 0.01;
const RealD            mpv      = 1.0;
const RealD            M5       = 1.8;

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is set up to use " << threads << " threads" << std::endl;

  // Initialize spacetime grid
  std::cout << GridLogMessage << "Lattice dimensions: " << grid_dim << "  Ls: " << Ls << std::endl;
  GridCartesian*         UGrid   = SpaceTimeGrid::makeFourDimGrid(grid_dim,
                                    GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian*         FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  // Set up RNGs
  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid);
  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);

  // Random gauge field
  LatticeGaugeField Umu(UGrid);
  SU3::HotConfiguration(RNG4, Umu);

  FermionAction::ImplParams params;
  FermionAction Lop(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mf,  mf, mpv,  0.0, -1, M5, b, c, params);
  FermionAction Rop(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mpv, mf, mpv, -1.0,  1, M5, b, c, params);

  // Construct the action and test the heatbath (zero initial guess)
  {
    OneFlavourRationalParams Params(0.95, 100.0, 5000, 1.0e-12, Npoles);
    ConjugateGradient<FermionField> CG(1.0e-12, 5000);
    ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> Meofa(Lop, Rop, CG, Params, false);

    Meofa.refresh(Umu, RNG5);
    printf("<Phi|Meofa|Phi> = %1.15e\n", Meofa.S(Umu));
  }

  // Construct the action and test the heatbath (forecasted initial guesses)
  {
    OneFlavourRationalParams Params(0.95, 100.0, 5000, 1.0e-12, Npoles);
    ConjugateGradient<FermionField> CG(1.0e-12, 5000);
    ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> Meofa(Lop, Rop, CG, Params, true);

    Meofa.refresh(Umu, RNG5);
    printf("<Phi|Meofa|Phi> = %1.15e\n", Meofa.S(Umu));
  }

  return 0;
}
