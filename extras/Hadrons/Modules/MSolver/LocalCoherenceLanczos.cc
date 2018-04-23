#include <Grid/Hadrons/Modules/MSolver/LocalCoherenceLanczos.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<ZFIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;

