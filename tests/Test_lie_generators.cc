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
  AdjointRep<3> AdjRep(grid);

  // AdjointRepresentation has the predefined number of colours Nc
  Representations<FundamentalRepresentation, AdjointRepresentation> RepresentationTypes(grid);  

  Grid_finalize();
}
