#include <Grid/Hadrons/Modules/MSolver/A2AVectors.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TA2AVectors<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class Grid::Hadrons::MSolver::TA2AVectors<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>;
