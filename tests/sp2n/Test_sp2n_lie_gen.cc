#include <Grid/Grid.h>

using namespace Grid;

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  std::cout << GridLogMessage << "*********************************************"
              << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(2) (print and test)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
              << std::endl;
      
  Sp2::printGenerators();
  Sp2::testGenerators();
  
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(4) (print and test)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
    
  Sp4::printGenerators();
  Sp4::testGenerators();
    
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(6) (test)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;

  Sp6::testGenerators();

  std::cout << GridLogMessage << "*********************************************"
              << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(8) (test)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
              << std::endl;
    
  Sp8::testGenerators();
    
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(4) TwoIndexAS (test)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;

  Sp_TwoIndex<4, AntiSymmetric>:::testGenerators();
    
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(6) TwoIndexAS (test)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;

  Sp_TwoIndex<6, AntiSymmetric>::testGenerators();

  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(8) TwoIndexAS (test)" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;

  Sp_TwoIndex<8, AntiSymmetric>::testGenerators();

  Grid_finalize();
}
