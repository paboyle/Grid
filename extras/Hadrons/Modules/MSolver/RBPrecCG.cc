#include <Grid/Hadrons/Modules/MSolver/RBPrecCG.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TRBPrecCG<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class Grid::Hadrons::MSolver::TRBPrecCG<ZFIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;

