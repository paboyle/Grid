/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonAdjointFermionGauge.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include "Grid/Grid.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid {
namespace QCD {

// Here change the allowed (higher) representations
typedef Representations< FundamentalRepresentation, AdjointRepresentation , TwoIndexSymmetricRepresentation> TheRepresentations;

class HmcRunner : public NerscHmcRunnerHirep< TheRepresentations > {
 public:
  void BuildTheAction(int argc, char **argv)

  {
    typedef WilsonAdjImplR AdjImplPolicy; // gauge field implemetation for the pseudofermions
    typedef WilsonAdjFermionR AdjFermionAction; // type of lattice fermions (Wilson, DW, ...)
    typedef WilsonTwoIndexSymmetricImplR SymmImplPolicy; 
    typedef WilsonTwoIndexSymmetricFermionR SymmFermionAction; 

 
    typedef typename AdjFermionAction::FermionField AdjFermionField;
    typedef typename SymmFermionAction::FermionField SymmFermionField;

    UGrid = SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
        GridDefaultMpi());
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

    FGrid = UGrid;
    FrbGrid = UrbGrid;

    // temporarily need a gauge field
    //LatticeGaugeField U(UGrid);
    AdjointRepresentation::LatticeField UA(UGrid);
    TwoIndexSymmetricRepresentation::LatticeField US(UGrid);

    // Gauge action
    WilsonGaugeActionR Waction(5.6);

    Real adjoint_mass = -0.1;
    Real symm_mass = -0.5;
    AdjFermionAction AdjFermOp(UA, *FGrid, *FrbGrid, adjoint_mass);
    SymmFermionAction SymmFermOp(US, *FGrid, *FrbGrid, symm_mass);

    ConjugateGradient<AdjFermionField> CG_adj(1.0e-8, 10000, false);
    ConjugateGradient<SymmFermionField> CG_symm(1.0e-8, 10000, false);

    // Pass two solvers: one for the force computation and one for the action
    TwoFlavourPseudoFermionAction<AdjImplPolicy> Nf2_Adj(AdjFermOp, CG_adj, CG_adj);
    TwoFlavourPseudoFermionAction<SymmImplPolicy> Nf2_Symm(SymmFermOp, CG_symm, CG_symm);

    // Collect actions
    ActionLevel<LatticeGaugeField, TheRepresentations > Level1(1);
    Level1.push_back(&Nf2_Adj);
    Level1.push_back(&Nf2_Symm);

    ActionLevel<LatticeGaugeField, TheRepresentations > Level2(4);
    Level2.push_back(&Waction);

    TheAction.push_back(Level1);
    TheAction.push_back(Level2);

    Run(argc, argv);
  };
};
}
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads
            << " threads" << std::endl;

  HmcRunner TheHMC;

  TheHMC.BuildTheAction(argc, argv);
}
