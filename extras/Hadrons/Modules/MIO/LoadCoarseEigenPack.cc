#include <Grid/Hadrons/Modules/MIO/LoadCoarseEigenPack.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>>;

