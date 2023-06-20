#include <Grid/Grid.h>

#include <iostream>

using namespace Grid;

template <int ngroup>
std::ostream& operator<<(std::ostream& o, Sp<ngroup> g) {
  return o << "Sp(" << ngroup << ") Fundamental";
}

template <int ngroup, TwoIndexSymmetry S>
std::ostream& operator<<(std::ostream& o, Sp_TwoIndex<ngroup, S> g) {
  return o << "Sp(" << ngroup << ") TwoIndex "
           << (S == Symmetric ? "Symmetric" : "AntiSymmetric");
}

template <class Group>
void run_check_on(bool print_generators = false) {
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;
  std::cout << GridLogMessage << "* Generators for " << Group() << std::endl;
  std::cout << GridLogMessage << "*********************************************"
            << std::endl;

  if (print_generators) {
    Group::printGenerators();
  }
  Group::testGenerators();
}

template <int ngroup>
void run_checks() {
  run_check_on<Sp<ngroup>>();
  run_check_on<Sp_TwoIndex<ngroup, Symmetric>>();
  run_check_on<Sp_TwoIndex<ngroup, AntiSymmetric>>();
}

template <>
void run_checks<2>() {
  // Print generators because they are small enough to be actually helpful.
  run_check_on<Sp<2>>(true);
  run_check_on<Sp_TwoIndex<2, Symmetric>>(true);
  // The AntiSymmetric representation is 0 dimensional. This makes problems in
  // device code.
}

template <>
void run_checks<4>() {
  // Print generators because they are small enough to be actually helpful.
  run_check_on<Sp<4>>(true);
  run_check_on<Sp_TwoIndex<4, Symmetric>>(true);
  run_check_on<Sp_TwoIndex<4, AntiSymmetric>>(true);
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  run_checks<2>();
  run_checks<4>();
  run_checks<6>();
  run_checks<8>();

  Grid_finalize();
}
