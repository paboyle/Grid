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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid.h>

#include <qcd/utils/CovariantCshift.h>
#include <qcd/utils/WilsonLoops.h>
#include <qcd/utils/SUn.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt({4,4,4,8});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(latt, 
							GridDefaultSimd(Nd,vComplex::Nsimd()),
							GridDefaultMpi());
  
  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  std::cout<<GridLogMessage<<"* Generators for SU(2)"<<std::endl;
  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  SU2::printGenerators();
  SU2::testGenerators();

  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  std::cout<<GridLogMessage<<"* Generators for SU(3)"<<std::endl;
  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  SU3::printGenerators();
  SU3::testGenerators();

  //  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  //  std::cout<<GridLogMessage<<"* Generators for SU(4)"<<std::endl;
  //  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  //  SU4::printGenerators();
  //  SU4::testGenerators();

  //  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  //  std::cout<<GridLogMessage<<"* Generators for SU(5)"<<std::endl;
  //  std::cout<<GridLogMessage<<"*********************************************"<<std::endl;
  //  SU5::printGenerators();
  //  SU5::testGenerators();


  Grid_finalize();
}


