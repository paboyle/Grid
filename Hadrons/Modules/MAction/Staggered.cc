#include <Hadrons/Modules/MAction/Staggered.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MAction;

template class Grid::Hadrons::MAction::TStaggered<STAGIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MAction::TStaggered<STAGIMPLF>;
#endif
