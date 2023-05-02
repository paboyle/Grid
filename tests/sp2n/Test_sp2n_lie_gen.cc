#include <Grid/Grid.h>

using namespace Grid;

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);
    
  //std::vector<int> latt({4, 4, 4, 8});
  //GridCartesian* grid = SpaceTimeGrid::makeFourDimGrid(
  //latt, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  //GridRedBlackCartesian* rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);
    
  std::cout << GridLogMessage << "*********************************************"
              << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(2)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
              << std::endl;
      
  Sp2::printGenerators();
  Sp2::testGenerators();
  
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(4)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
    
  Sp4::printGenerators();
  Sp4::testGenerators();
    
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(6)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  
  Sp6::printGenerators();
  Sp6::testGenerators();
    
  std::cout << GridLogMessage << "*********************************************"
              << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(8)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
              << std::endl;
    
  Sp8::printGenerators();
  Sp8::testGenerators(); 
    
  Grid_finalize();
}
