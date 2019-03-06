/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MNoise/FullVolumeSpinColorDiagonal.cc

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
Author: Vera Guelpers <vmg1n14@soton.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Hadrons/Modules/MNoise/FullVolumeSpinColorDiagonal.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MNoise;

template class Grid::Hadrons::MNoise::TFullVolumeSpinColorDiagonal<FIMPL>;
template class Grid::Hadrons::MNoise::TFullVolumeSpinColorDiagonal<ZFIMPL>;
