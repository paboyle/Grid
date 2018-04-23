#include <Grid/Hadrons/Modules/MIO/LoadBinary.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

template class Grid::Hadrons::MIO::TLoadBinary<GIMPL>;
template class Grid::Hadrons::MIO::TLoadBinary<ScalarNxNAdjImplR<2>>;
template class Grid::Hadrons::MIO::TLoadBinary<ScalarNxNAdjImplR<3>>;
template class Grid::Hadrons::MIO::TLoadBinary<ScalarNxNAdjImplR<4>>;
template class Grid::Hadrons::MIO::TLoadBinary<ScalarNxNAdjImplR<5>>;
template class Grid::Hadrons::MIO::TLoadBinary<ScalarNxNAdjImplR<6>>;

