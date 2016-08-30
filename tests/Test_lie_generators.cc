/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_lie_generators.cc

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid.h>

#include <qcd/utils/CovariantCshift.h>
#include <qcd/utils/SUn.h>
#include <qcd/utils/SUnAdjoint.h>
#include <qcd/representations/adjoint.h>
#include <qcd/utils/WilsonLoops.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  std::vector<int> latt({4, 4, 4, 8});
  GridCartesian* grid = SpaceTimeGrid::makeFourDimGrid(
      latt, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());

  GridRedBlackCartesian* rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for SU(2)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  SU2::printGenerators();
  std::cout << "Dimension of adjoint representation: "<< SU2Adjoint::Dimension << std::endl;
  SU2Adjoint::printGenerators();
  SU2::testGenerators();
  SU2Adjoint::testGenerators();

  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for SU(3)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  SU3::printGenerators();
  std::cout << "Dimension of adjoint representation: "<< SU3Adjoint::Dimension << std::endl;
  SU3Adjoint::printGenerators();
  SU3::testGenerators();
  SU3Adjoint::testGenerators();

  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  std::cout<<GridLogMessage<<"* Generators for SU(4)"<<std::endl;
  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  SU4::printGenerators();
  std::cout << "Dimension of adjoint representation: "<< SU4Adjoint::Dimension << std::endl;
  SU4Adjoint::printGenerators();
  SU4::testGenerators();
  SU4Adjoint::testGenerators();

  //  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  //  std::cout<<GridLogMessage<<"* Generators for SU(5)"<<std::endl;
  //  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  //  SU5::printGenerators();
  //  SU5::testGenerators();


  // Projectors 
  GridParallelRNG gridRNG(grid);
  gridRNG.SeedRandomDevice();
  SU3Adjoint::LatticeAdjMatrix Gauss(grid);
  SU3::LatticeAlgebraVector ha(grid);
  SU3::LatticeAlgebraVector hb(grid);
  random(gridRNG,Gauss);

  std::cout << GridLogMessage << "Start projectOnAlgebra" << std::endl;
  SU3Adjoint::projectOnAlgebra(ha, Gauss);
  std::cout << GridLogMessage << "end projectOnAlgebra" << std::endl;
  std::cout << GridLogMessage << "Start projector" << std::endl;
  SU3Adjoint::projector(hb, Gauss);
  std::cout << GridLogMessage << "end projector" << std::endl;

  std::cout << GridLogMessage << "ReStart projector" << std::endl;
  SU3Adjoint::projector(hb, Gauss);
  std::cout << GridLogMessage << "end projector" << std::endl;
  SU3::LatticeAlgebraVector diff = ha -hb;
  std::cout << GridLogMessage << "Difference: " << norm2(diff) << std::endl;


  // Testing HMC representation classes
  AdjointRep<Nc> AdjRep(grid);

  // AdjointRepresentation has the predefined number of colours Nc
  Representations<FundamentalRepresentation, AdjointRepresentation> RepresentationTypes(grid);  

  LatticeGaugeField U(grid), V(grid);
  SU<Nc>::HotConfiguration<LatticeGaugeField>(gridRNG, U);
  SU<Nc>::HotConfiguration<LatticeGaugeField>(gridRNG, V);

  // Test group structure
  // (U_f * V_f)_r = U_r * V_r
  LatticeGaugeField UV(grid);
  UV = zero;
  for (int mu = 0; mu < Nd; mu++) {
    SU<Nc>::LatticeMatrix Umu = peekLorentz(U,mu);
    SU<Nc>::LatticeMatrix Vmu = peekLorentz(V,mu);
    pokeLorentz(UV,Umu*Vmu, mu);
  }

  AdjRep.update_representation(UV);
  typename AdjointRep<Nc>::LatticeField UVr = AdjRep.U;  // (U_f * V_f)_r


  AdjRep.update_representation(U);
  typename AdjointRep<Nc>::LatticeField Ur = AdjRep.U;  // U_r

  AdjRep.update_representation(V);
  typename AdjointRep<Nc>::LatticeField Vr = AdjRep.U;  // V_r

  typename AdjointRep<Nc>::LatticeField UrVr(grid);
  UrVr = zero;
  for (int mu = 0; mu < Nd; mu++) {
    typename AdjointRep<Nc>::LatticeMatrix Urmu = peekLorentz(Ur,mu);
    typename AdjointRep<Nc>::LatticeMatrix Vrmu = peekLorentz(Vr,mu);
    pokeLorentz(UrVr,Urmu*Vrmu, mu);
  }

  typename AdjointRep<Nc>::LatticeField Diff_check = UVr - UrVr;
  std::cout << GridLogMessage << "Group structure SU("<<Nc<<") check difference : " << norm2(Diff_check) << std::endl;

  // Check correspondence of algebra and group transformations
  // Create a random vector
  SU<Nc>::LatticeAlgebraVector h_adj(grid);
  typename AdjointRep<Nc>::LatticeMatrix Ar(grid);
  random(gridRNG,h_adj);
  h_adj = real(h_adj);
  SU_Adjoint<Nc>::AdjointLieAlgebraMatrix(h_adj,Ar);

  // Re-extract h_adj
  SU<Nc>::LatticeAlgebraVector h_adj2(grid);
  SU_Adjoint<Nc>::projectOnAlgebra(h_adj2, Ar);
  SU<Nc>::LatticeAlgebraVector h_diff = h_adj - h_adj2;
  std::cout << GridLogMessage << "Projections structure check vector difference : " << norm2(h_diff) << std::endl;

  // Exponentiate
  typename AdjointRep<Nc>::LatticeMatrix Uadj(grid);
  Uadj  = expMat(Ar, 1.0, 16);

  typename AdjointRep<Nc>::LatticeMatrix uno(grid);
  uno = 1.0;
  // Check matrix Uadj, must be real orthogonal
  typename AdjointRep<Nc>::LatticeMatrix Ucheck = Uadj - conjugate(Uadj);
  std::cout << GridLogMessage << "Reality check: " << norm2(Ucheck)
            << std::endl;

  Ucheck = Uadj * adj(Uadj) - uno;
  std::cout << GridLogMessage << "orthogonality check 1: " << norm2(Ucheck)
            << std::endl;
  Ucheck = adj(Uadj) * Uadj - uno;
  std::cout << GridLogMessage << "orthogonality check 2: " << norm2(Ucheck)
            << std::endl;
      

  // Construct the fundamental matrix in the group
  SU<Nc>::LatticeMatrix Af(grid);
  SU<Nc>::FundamentalLieAlgebraMatrix(h_adj,Af);
  SU<Nc>::LatticeMatrix Ufund(grid);
  Ufund  = expMat(Af, 1.0, 16);
  // Check unitarity
  SU<Nc>::LatticeMatrix uno_f(grid);
  uno_f = 1.0;
  SU<Nc>::LatticeMatrix UnitCheck(grid);
  UnitCheck = Ufund * adj(Ufund) - uno_f;
  std::cout << GridLogMessage << "unitarity check 1: " << norm2(UnitCheck)
            << std::endl;
  UnitCheck = adj(Ufund) * Ufund - uno_f;
  std::cout << GridLogMessage << "unitarity check 2: " << norm2(UnitCheck)
            << std::endl;


  // Tranform to the adjoint representation
  U = zero; // fill this with only one direction
  pokeLorentz(U,Ufund,0); // the representation transf acts on full gauge fields

  AdjRep.update_representation(U);
  Ur = AdjRep.U;  // U_r  
  typename AdjointRep<Nc>::LatticeMatrix Ur0 = peekLorentz(Ur,0); // this should be the same as Uadj

  typename AdjointRep<Nc>::LatticeMatrix Diff_check_mat = Ur0 - Uadj;
  std::cout << GridLogMessage << "Projections structure check group difference : " << norm2(Diff_check_mat) << std::endl;

  Grid_finalize();
}
