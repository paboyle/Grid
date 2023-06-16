#include <Grid/Grid.h>

using namespace Grid;

template <int ncolour>
void run_checks(bool print_generators = 0) {
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for Sp(" << ncolour << ")"
            << "Fundamental" << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;

  if (print_generators) {
    Sp<ncolour>::printGenerators();
  }
  Sp<ncolour>::testGenerators();

  if (Sp_TwoIndex<ncolour, Symmetric>::Dimension > 1) {
    std::cout << GridLogMessage
              << "*********************************************" << std::endl;
    std::cout << GridLogMessage << "* Generators for Sp(" << ncolour << ")"
              << "TwoIndex Symmetric: " << std::endl;
    std::cout << GridLogMessage
              << "*********************************************" << std::endl;
    if (print_generators) {
      Sp_TwoIndex<ncolour, Symmetric>::printGenerators();
    }
    Sp_TwoIndex<ncolour, Symmetric>::testGenerators();
  }

  if (Sp_TwoIndex<ncolour, AntiSymmetric>::Dimension > 1) {
    std::cout << GridLogMessage
              << "*********************************************" << std::endl;
    std::cout << GridLogMessage << "* Generators for Sp(" << ncolour << ")"
              << "TwoIndex AntiSymmetric: " << std::endl;
    std::cout << GridLogMessage
              << "*********************************************" << std::endl;
    if (print_generators) {
      Sp_TwoIndex<ncolour, AntiSymmetric>::printGenerators();
    }
    Sp_TwoIndex<ncolour, AntiSymmetric>::testGenerators();
  }
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  run_checks<2>(1);  //  check and print Nc=2
  run_checks<4>(1);  //  check and print Nc=4
  run_checks<6>();   //  check Nc=6
  run_checks<8>();   //  check Nc=8

  Grid_finalize();
}
