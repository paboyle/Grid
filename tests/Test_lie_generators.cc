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
							GridDefaultSimd(Nd,vComplexF::Nsimd()),
							GridDefaultMpi());
  
  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  std::cout<<"*********************************************"<<std::endl;
  std::cout<<"* Generators for SU(2)"<<std::endl;
  std::cout<<"*********************************************"<<std::endl;
  SU2::printGenerators();
  SU2::testGenerators();

  std::cout<<"*********************************************"<<std::endl;
  std::cout<<"* Generators for SU(3)"<<std::endl;
  std::cout<<"*********************************************"<<std::endl;
  SU3::printGenerators();
  SU3::testGenerators();

  //  std::cout<<"*********************************************"<<std::endl;
  //  std::cout<<"* Generators for SU(4)"<<std::endl;
  //  std::cout<<"*********************************************"<<std::endl;
  //  SU4::printGenerators();
  //  SU4::testGenerators();

  //  std::cout<<"*********************************************"<<std::endl;
  //  std::cout<<"* Generators for SU(5)"<<std::endl;
  //  std::cout<<"*********************************************"<<std::endl;
  //  SU5::printGenerators();
  //  SU5::testGenerators();


  Grid_finalize();
}


