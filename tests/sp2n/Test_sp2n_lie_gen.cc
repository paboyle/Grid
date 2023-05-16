#include <Grid/Grid.h>

using namespace Grid;

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);
    
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
    
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(4) TwoIndex AntiSymmetric" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;

  Sp4TwoIndexAntiSymm::printGenerators();
  Sp4TwoIndexAntiSymm::testGenerators();

  Grid_finalize();
}
