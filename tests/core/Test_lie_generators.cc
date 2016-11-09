/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_lie_generators.cc

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
#include <Grid/Grid.h>

#include <Grid/qcd/utils/CovariantCshift.h>

#include <Grid/qcd/utils/SUn.h>
#include <Grid/qcd/utils/SUnAdjoint.h>
#include <Grid/qcd/utils/SUnTwoIndex.h>

#include <Grid/qcd/representations/adjoint.h>
#include <Grid/qcd/representations/two_index.h>
#include <Grid/qcd/utils/WilsonLoops.h>

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
  Representations<FundamentalRepresentation, AdjointRepresentation, TwoIndexSymmetricRepresentation> RepresentationTypes(grid);  

  
  LatticeGaugeField U(grid), V(grid);
  SU<Nc>::HotConfiguration<LatticeGaugeField>(gridRNG, U);
  SU<Nc>::HotConfiguration<LatticeGaugeField>(gridRNG, V);

  // Adjoint representation
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
  std::cout << GridLogMessage << "Group structure SU("<<Nc<<") check difference (Adjoint representation) : " << norm2(Diff_check) << std::endl;

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
  std::cout << GridLogMessage << "Projections structure check vector difference (Adjoint representation) : " << norm2(h_diff) << std::endl;

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




  // TwoIndexRep tests

  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;


  
  std::cout << GridLogMessage << "* eS^{ij} base for SU(2)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "Dimension of Two Index Symmetric representation: "<< SU2TwoIndexSymm::Dimension << std::endl;
  SU2TwoIndexSymm::printBase();
      std::cout << GridLogMessage << "*********************************************"
            << std::endl;
        std::cout << GridLogMessage << "Generators of Two Index Symmetric representation: "<< SU2TwoIndexSymm::Dimension << std::endl;
  SU2TwoIndexSymm::printGenerators();
        std::cout << GridLogMessage << "Test of Two Index Symmetric Generators: "<< SU2TwoIndexSymm::Dimension << std::endl;
  SU2TwoIndexSymm::testGenerators();
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;


  
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* eAS^{ij} base for SU(2)" << std::endl;
  
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "Dimension of Two Index anti-Symmetric representation: "<< SU2TwoIndexAntiSymm::Dimension << std::endl;
  SU2TwoIndexAntiSymm::printBase();
      std::cout << GridLogMessage << "*********************************************"
            << std::endl;
        std::cout << GridLogMessage << "Dimension of Two Index anti-Symmetric representation: "<< SU2TwoIndexAntiSymm::Dimension << std::endl;
  SU2TwoIndexAntiSymm::printGenerators();
  std::cout << GridLogMessage << "Test of Two Index anti-Symmetric Generators: "<< SU2TwoIndexAntiSymm::Dimension << std::endl;
  SU2TwoIndexAntiSymm::testGenerators();
  
  
  
  std::cout << GridLogMessage << "*********************************************"
      << std::endl;
  std::cout << GridLogMessage << "Test for the Two Index Symmetric projectors"
      << std::endl;
  // Projectors 
  SU3TwoIndexSymm::LatticeTwoIndexMatrix Gauss2(grid);
  random(gridRNG,Gauss2);
  
  std::cout << GridLogMessage << "Start projectOnAlgebra" << std::endl;
  SU3TwoIndexSymm::projectOnAlgebra(ha, Gauss2);
  std::cout << GridLogMessage << "end projectOnAlgebra" << std::endl;
  std::cout << GridLogMessage << "Start projector" << std::endl;
  SU3TwoIndexSymm::projector(hb, Gauss2);
  std::cout << GridLogMessage << "end projector" << std::endl;
  
  std::cout << GridLogMessage << "ReStart projector" << std::endl;
  SU3TwoIndexSymm::projector(hb, Gauss2);
  std::cout << GridLogMessage << "end projector" << std::endl;
  SU3::LatticeAlgebraVector diff2 = ha - hb;
  std::cout << GridLogMessage << "Difference: " << norm2(diff) << std::endl;
  std::cout << GridLogMessage << "*********************************************"
      << std::endl;

  
    std::cout << GridLogMessage << "*********************************************"
      << std::endl;
  std::cout << GridLogMessage << "Test for the Two index anti-Symmetric projectors"
      << std::endl;
  // Projectors
  SU3TwoIndexAntiSymm::LatticeTwoIndexMatrix Gauss2a(grid);
  random(gridRNG,Gauss2a);
  
  std::cout << GridLogMessage << "Start projectOnAlgebra" << std::endl;
  SU3TwoIndexAntiSymm::projectOnAlgebra(ha, Gauss2a);
  std::cout << GridLogMessage << "end projectOnAlgebra" << std::endl;
  std::cout << GridLogMessage << "Start projector" << std::endl;
  SU3TwoIndexAntiSymm::projector(hb, Gauss2a);
  std::cout << GridLogMessage << "end projector" << std::endl;
  
  std::cout << GridLogMessage << "ReStart projector" << std::endl;
  SU3TwoIndexAntiSymm::projector(hb, Gauss2a);
  std::cout << GridLogMessage << "end projector" << std::endl;
  SU3::LatticeAlgebraVector diff2a = ha - hb;
  std::cout << GridLogMessage << "Difference: " << norm2(diff2a) << std::endl;
  std::cout << GridLogMessage << "*********************************************"
      << std::endl;
    
  
  std::cout << GridLogMessage << "Two index Symmetric: Checking Group Structure"
      << std::endl;
  // Testing HMC representation classes
  TwoIndexRep< Nc, Symmetric > TIndexRep(grid);

  // Test group structure
  // (U_f * V_f)_r = U_r * V_r
  LatticeGaugeField U2(grid), V2(grid);
  SU<Nc>::HotConfiguration<LatticeGaugeField>(gridRNG, U2);
  SU<Nc>::HotConfiguration<LatticeGaugeField>(gridRNG, V2);
  
  LatticeGaugeField UV2(grid);
  UV2 = zero;
  for (int mu = 0; mu < Nd; mu++) {
    SU<Nc>::LatticeMatrix Umu2 = peekLorentz(U2,mu);
    SU<Nc>::LatticeMatrix Vmu2 = peekLorentz(V2,mu);
    pokeLorentz(UV2,Umu2*Vmu2, mu);
  }
  
  TIndexRep.update_representation(UV2);
  typename TwoIndexRep< Nc, Symmetric >::LatticeField UVr2 = TIndexRep.U;  // (U_f * V_f)_r
  
  TIndexRep.update_representation(U2);
  typename TwoIndexRep< Nc, Symmetric >::LatticeField Ur2 = TIndexRep.U;  // U_r
  
  TIndexRep.update_representation(V2);
  typename TwoIndexRep< Nc, Symmetric >::LatticeField Vr2 = TIndexRep.U;  // V_r
  
  typename TwoIndexRep< Nc, Symmetric >::LatticeField Ur2Vr2(grid);
  Ur2Vr2 = zero;
  for (int mu = 0; mu < Nd; mu++) {
    typename TwoIndexRep< Nc, Symmetric >::LatticeMatrix Urmu2 = peekLorentz(Ur2,mu);
    typename TwoIndexRep< Nc, Symmetric >::LatticeMatrix Vrmu2 = peekLorentz(Vr2,mu);
    pokeLorentz(Ur2Vr2,Urmu2*Vrmu2, mu);
  }
  
  typename TwoIndexRep< Nc, Symmetric >::LatticeField Diff_check2 = UVr2 - Ur2Vr2;
  std::cout << GridLogMessage << "Group structure SU("<<Nc<<") check difference (Two Index Symmetric): " << norm2(Diff_check2) << std::endl;

  
  // Check correspondence of algebra and group transformations
  // Create a random vector
  SU<Nc>::LatticeAlgebraVector h_sym(grid);
  typename TwoIndexRep< Nc, Symmetric>::LatticeMatrix Ar_sym(grid);
  random(gridRNG,h_sym);
  h_sym = real(h_sym);
  SU_TwoIndex<Nc,Symmetric>::TwoIndexLieAlgebraMatrix(h_sym,Ar_sym);
  
  // Re-extract h_sym
  SU<Nc>::LatticeAlgebraVector h_sym2(grid);
  SU_TwoIndex< Nc, Symmetric>::projectOnAlgebra(h_sym2, Ar_sym);
  SU<Nc>::LatticeAlgebraVector h_diff_sym = h_sym - h_sym2;
  std::cout << GridLogMessage << "Projections structure check vector difference (Two Index Symmetric): " << norm2(h_diff_sym) << std::endl;

  
  // Exponentiate
  typename TwoIndexRep< Nc, Symmetric>::LatticeMatrix U2iS(grid);
  U2iS  = expMat(Ar_sym, 1.0, 16);
  
  typename TwoIndexRep< Nc, Symmetric>::LatticeMatrix uno2iS(grid);
  uno2iS = 1.0;
  // Check matrix U2iS, must be real orthogonal
  typename TwoIndexRep< Nc, Symmetric>::LatticeMatrix Ucheck2iS = U2iS - conjugate(U2iS);
  std::cout << GridLogMessage << "Reality check: " << norm2(Ucheck2iS)
      << std::endl;
  
  Ucheck2iS = U2iS * adj(U2iS) - uno2iS;
  std::cout << GridLogMessage << "orthogonality check 1: " << norm2(Ucheck2iS)
      << std::endl;
  Ucheck2iS = adj(U2iS) * U2iS - uno2iS;
  std::cout << GridLogMessage << "orthogonality check 2: " << norm2(Ucheck2iS)
      << std::endl;
  
  
  
  // Construct the fundamental matrix in the group
  SU<Nc>::LatticeMatrix Af_sym(grid);
  SU<Nc>::FundamentalLieAlgebraMatrix(h_sym,Af_sym);
  SU<Nc>::LatticeMatrix Ufund2(grid);
  Ufund2  = expMat(Af_sym, 1.0, 16);
  SU<Nc>::LatticeMatrix UnitCheck2(grid);
  UnitCheck2 = Ufund2 * adj(Ufund2) - uno_f;
  std::cout << GridLogMessage << "unitarity check 1: " << norm2(UnitCheck2)
      << std::endl;
  UnitCheck2 = adj(Ufund2) * Ufund2 - uno_f;
  std::cout << GridLogMessage << "unitarity check 2: " << norm2(UnitCheck2)
      << std::endl;
  

  // Tranform to the 2Index Sym representation
  U = zero; // fill this with only one direction
  pokeLorentz(U,Ufund2,0); // the representation transf acts on full gauge fields
  
  TIndexRep.update_representation(U);
  Ur2 = TIndexRep.U;  // U_r  
  typename TwoIndexRep< Nc, Symmetric>::LatticeMatrix Ur02 = peekLorentz(Ur2,0); // this should be the same as U2iS
  
  typename TwoIndexRep< Nc, Symmetric>::LatticeMatrix Diff_check_mat2 = Ur02 - U2iS;
  std::cout << GridLogMessage << "Projections structure check group difference (Two Index Symmetric): " << norm2(Diff_check_mat2) << std::endl;
  



  if (TwoIndexRep<Nc, AntiSymmetric >::Dimension != 1){

  std::cout << GridLogMessage << "*********************************************"
      << std::endl;
    
  
  std::cout << GridLogMessage << "Two Index anti-Symmetric: Check Group Structure"
      << std::endl;
  // Testing HMC representation classes
  TwoIndexRep< Nc, AntiSymmetric > TIndexRepA(grid);


  // Test group structure
  // (U_f * V_f)_r = U_r * V_r
  LatticeGaugeField U2A(grid), V2A(grid);
  SU<Nc>::HotConfiguration<LatticeGaugeField>(gridRNG, U2A);
  SU<Nc>::HotConfiguration<LatticeGaugeField>(gridRNG, V2A);
  
  LatticeGaugeField UV2A(grid);
  UV2A = zero;
  for (int mu = 0; mu < Nd; mu++) {
    SU<Nc>::LatticeMatrix Umu2A = peekLorentz(U2,mu);
    SU<Nc>::LatticeMatrix Vmu2A = peekLorentz(V2,mu);
    pokeLorentz(UV2A,Umu2A*Vmu2A, mu);
  }
  
  TIndexRep.update_representation(UV2A);
  typename TwoIndexRep< Nc, AntiSymmetric >::LatticeField UVr2A = TIndexRepA.U;  // (U_f * V_f)_r
  
  TIndexRep.update_representation(U2A);
  typename TwoIndexRep< Nc, AntiSymmetric >::LatticeField Ur2A = TIndexRepA.U;  // U_r
  
  TIndexRep.update_representation(V2A);
  typename TwoIndexRep< Nc, AntiSymmetric >::LatticeField Vr2A = TIndexRepA.U;  // V_r
  
  typename TwoIndexRep< Nc, AntiSymmetric >::LatticeField Ur2Vr2A(grid);
  Ur2Vr2A = zero;
  for (int mu = 0; mu < Nd; mu++) {
    typename TwoIndexRep< Nc, AntiSymmetric >::LatticeMatrix Urmu2A = peekLorentz(Ur2A,mu);
    typename TwoIndexRep< Nc, AntiSymmetric >::LatticeMatrix Vrmu2A = peekLorentz(Vr2A,mu);
    pokeLorentz(Ur2Vr2A,Urmu2A*Vrmu2A, mu);
  }
  
  typename TwoIndexRep< Nc, AntiSymmetric >::LatticeField Diff_check2A = UVr2A - Ur2Vr2A;
  std::cout << GridLogMessage << "Group structure SU("<<Nc<<") check difference (Two Index anti-Symmetric): " << norm2(Diff_check2A) << std::endl;

  
  // Check correspondence of algebra and group transformations
  // Create a random vector
  SU<Nc>::LatticeAlgebraVector h_Asym(grid);
  typename TwoIndexRep< Nc, AntiSymmetric>::LatticeMatrix Ar_Asym(grid);
  random(gridRNG,h_Asym);
  h_Asym = real(h_Asym);
  SU_TwoIndex< Nc, AntiSymmetric>::TwoIndexLieAlgebraMatrix(h_Asym,Ar_Asym);
  
  // Re-extract h_sym
  SU<Nc>::LatticeAlgebraVector h_Asym2(grid);
  SU_TwoIndex< Nc, AntiSymmetric>::projectOnAlgebra(h_Asym2, Ar_Asym);
  SU<Nc>::LatticeAlgebraVector h_diff_Asym = h_Asym - h_Asym2;
  std::cout << GridLogMessage << "Projections structure check vector difference (Two Index anti-Symmetric): " << norm2(h_diff_Asym) << std::endl;

  
  // Exponentiate
  typename TwoIndexRep< Nc, AntiSymmetric>::LatticeMatrix U2iAS(grid);
  U2iAS  = expMat(Ar_Asym, 1.0, 16);
  
  typename TwoIndexRep< Nc, AntiSymmetric>::LatticeMatrix uno2iAS(grid);
  uno2iAS = 1.0;
  // Check matrix U2iS, must be real orthogonal
  typename TwoIndexRep< Nc, AntiSymmetric>::LatticeMatrix Ucheck2iAS = U2iAS - conjugate(U2iAS);
  std::cout << GridLogMessage << "Reality check: " << norm2(Ucheck2iAS)
      << std::endl;
  
  Ucheck2iAS = U2iAS * adj(U2iAS) - uno2iAS;
  std::cout << GridLogMessage << "orthogonality check 1: " << norm2(Ucheck2iAS)
      << std::endl;
  Ucheck2iAS = adj(U2iAS) * U2iAS - uno2iAS;
  std::cout << GridLogMessage << "orthogonality check 2: " << norm2(Ucheck2iAS)
      << std::endl;
  
  
  
  // Construct the fundamental matrix in the group
  SU<Nc>::LatticeMatrix Af_Asym(grid);
  SU<Nc>::FundamentalLieAlgebraMatrix(h_Asym,Af_Asym);
  SU<Nc>::LatticeMatrix Ufund2A(grid);
  Ufund2A  = expMat(Af_Asym, 1.0, 16);
  SU<Nc>::LatticeMatrix UnitCheck2A(grid);
  UnitCheck2A = Ufund2A * adj(Ufund2A) - uno_f;
  std::cout << GridLogMessage << "unitarity check 1: " << norm2(UnitCheck2A)
      << std::endl;
  UnitCheck2A = adj(Ufund2A) * Ufund2A - uno_f;
  std::cout << GridLogMessage << "unitarity check 2: " << norm2(UnitCheck2A)
      << std::endl;
  

  // Tranform to the 2Index Sym representation
  U = zero; // fill this with only one direction
  pokeLorentz(U,Ufund2A,0); // the representation transf acts on full gauge fields
  
  TIndexRepA.update_representation(U);
  Ur2A = TIndexRepA.U;  // U_r  
  typename TwoIndexRep< Nc, AntiSymmetric>::LatticeMatrix Ur02A = peekLorentz(Ur2A,0); // this should be the same as U2iS
  
  typename TwoIndexRep< Nc, AntiSymmetric>::LatticeMatrix Diff_check_mat2A = Ur02A - U2iAS;
  std::cout << GridLogMessage << "Projections structure check group difference (Two Index anti-Symmetric): " << norm2(Diff_check_mat2A) << std::endl;
  
} else  {
  std::cout << GridLogMessage << "Skipping Two Index anti-Symmetric tests "
                                 "because representation is trivial (dim = 1)"
            << std::endl;
}









  
  Grid_finalize();
}
