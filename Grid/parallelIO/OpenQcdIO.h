/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/parallelIO/OpenQcdIO.h

Copyright (C) 2015 - 2020

Author: Daniel Richtmann <daniel.richtmann@ur.de>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#pragma once

NAMESPACE_BEGIN(Grid);

struct OpenQcdHeader : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(OpenQcdHeader,
                                  int,    Nt,
                                  int,    Nx,
                                  int,    Ny,
                                  int,    Nz,
                                  double, plaq);
};

class OpenQcdIO : public BinaryIO {
public:
  static constexpr double normalisationFactor = Nc; // normalisation difference: grid 18, openqcd 6

  static inline int readHeader(std::string file, GridBase* grid, FieldMetaData& field) {
    OpenQcdHeader header;

    {
      std::ifstream fin(file, std::ios::in | std::ios::binary);
      fin.read(reinterpret_cast<char*>(&header), sizeof(OpenQcdHeader));
      assert(!fin.fail());
      field.data_start = fin.tellg();
      fin.close();
    }

    header.plaq /= normalisationFactor;

    // sanity check (should trigger on endian issues)
    assert(0 < header.Nt && header.Nt <= 1024);
    assert(0 < header.Nx && header.Nx <= 1024);
    assert(0 < header.Ny && header.Ny <= 1024);
    assert(0 < header.Nz && header.Nz <= 1024);

    field.dimension[0] = header.Nx;
    field.dimension[1] = header.Ny;
    field.dimension[2] = header.Nz;
    field.dimension[3] = header.Nt;

    assert(grid->_ndimension == Nd);
    for(int d = 0; d < Nd; d++)
      assert(grid->_fdimensions[d] == field.dimension[d]);

    field.plaquette = header.plaq;

    return field.data_start;
  }

  template<class vsimd>
  static inline void readConfiguration(Lattice<iLorentzColourMatrix<vsimd>>& Umu,
                                       FieldMetaData&                        header,
                                       std::string                           file) {
    auto grid = dynamic_cast<GridCartesian*>(Umu.Grid());
    assert(grid != nullptr);
    assert((grid->_ndimension == Nd) && (Nd == 4));

    uint64_t offset = readHeader(file, Umu.Grid(), header);
    FieldMetaData clone(header);

    // NOTE: This version is suboptimal because we read in the full file on every rank
    std::vector<ColourMatrix> data(grid->gSites() * 4);
    {
      auto fin = std::fstream(file, std::ios::in | std::ios::binary);
      fin.seekg(offset);
      fin.read((char *)data.data(), data.size() * sizeof(ColourMatrix));
      fin.close();
    }

    // global lattice size
    Coordinate fdim = grid->FullDimensions();

    // coordinate of this process
    Coordinate pcoor;
    grid->ProcessorCoorFromRank(CartesianCommunicator::RankWorld(), pcoor);

    // loop over local indices
    thread_for(idx, grid->lSites(), {
      // convert local index to global coordinate
      Coordinate lcoor, gcoor;
      grid->LocalIndexToLocalCoor(idx, lcoor);
      grid->ProcessorCoorLocalCoorToGlobalCoor(pcoor, lcoor, gcoor);

      // openQCD stores links attached to odd sites
      bool neg = (gcoor[Xdir] + gcoor[Ydir] + gcoor[Zdir] + gcoor[Tdir]) % 2 != 1;

      LorentzColourMatrix site_data;
      for (int mu = 0; mu < 4; ++mu) {
        // determine the site at which it is stored
        Coordinate c = gcoor;
        if (neg)
          c[mu] = (c[mu] + 1) % grid->FullDimensions()[mu];

        // site-index in the OpenQCD format (which uses t,x,y,z order)
        int openqcd_idx = (c[Tdir] * fdim[Xdir] * fdim[Ydir] * fdim[Zdir]
                        +  c[Xdir] * fdim[Ydir] * fdim[Zdir]
                        +  c[Ydir] * fdim[Zdir]
                        +  c[Zdir])/2;
        int openqcd_mu = (mu + 1) % 4;

        // pick the colour-matrix out
        site_data(mu) =
            data[8 * openqcd_idx + 2 * openqcd_mu + (neg ? 1 : 0)]();
      }

      pokeLocalSite(site_data, Umu, lcoor);
    });

    GaugeStatistics(Umu, clone);

    std::cout << GridLogMessage << "OpenQcd Configuration " << file << " plaquette "
              << std::setprecision(15)
              << clone.plaquette << " header " << header.plaquette
              << " difference " << fabs(clone.plaquette - header.plaquette)
              << std::endl;

    if(fabs(clone.plaquette - header.plaquette) >= 1.0e-5) std::cout << " Plaquette mismatch " << std::endl;
    assert(fabs(clone.plaquette - header.plaquette) < 1.0e-5);

    std::cout << GridLogMessage << "OpenQcd Configuration " << file << " and plaquette agree" << std::endl;
  }
};

NAMESPACE_END(Grid);
