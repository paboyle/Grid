#include <Grid/Hadrons/Modules/MSolver/A2AVectors.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TA2AVectors<FIMPL, FermionEigenPack<FIMPL>>;
template class Grid::Hadrons::MSolver::TA2AVectors<ZFIMPL, FermionEigenPack<ZFIMPL>>;
