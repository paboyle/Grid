#include <Hadrons/Modules/MAction/ImprovedStaggered.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MAction;

template class Grid::Hadrons::MAction::TImprovedStaggered<StagIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MAction::TImprovedStaggered<StagIMPLF>;
#endif
