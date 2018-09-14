#include <Hadrons/Modules/MSolver/MixedPrecisionRBPrecCG.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TMixedPrecisionRBPrecCG<FIMPLF, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class Grid::Hadrons::MSolver::TMixedPrecisionRBPrecCG<ZFIMPLF, ZFIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>;
