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

    std::cout << GridLogDebug << "header: " << header << std::endl;
    std::cout << GridLogDebug << "grid dimensions: " << grid->_fdimensions << std::endl;
    std::cout << GridLogDebug << "file dimensions: " << field.dimension << std::endl;

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
    typedef Lattice<iDoubleStoredColourMatrix<vsimd>> DoubleStoredGaugeField;

    assert(Ns == 4 and Nd == 4 and Nc == 3);

    auto grid = dynamic_cast<GridCartesian*>(Umu.Grid());
    assert(grid != nullptr); assert(grid->_ndimension == Nd);

    uint64_t offset = readHeader(file, Umu.Grid(), header);

    FieldMetaData clone(header);

    std::string format("IEEE64"); // they always store little endian double precsision
    uint32_t    nersc_csum, scidac_csuma, scidac_csumb;

    GridCartesian*         grid_openqcd = createOpenQcdGrid(grid);
    GridRedBlackCartesian* grid_rb      = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

    typedef DoubleStoredColourMatrixD                                              fobj;
    typedef typename DoubleStoredGaugeField::vector_object::scalar_object          sobj;
    typedef typename DoubleStoredGaugeField::vector_object::Realified::scalar_type word;

    word w = 0;

    std::vector<fobj> iodata(grid_openqcd->lSites()); // Munge, checksum, byte order in here
    std::vector<sobj> scalardata(grid->lSites());

    IOobject(w, grid_openqcd, iodata, file, offset, format, BINARYIO_READ | BINARYIO_LEXICOGRAPHIC,
             nersc_csum, scidac_csuma, scidac_csumb);

    GridStopWatch timer;
    timer.Start();

    DoubleStoredGaugeField Umu_ds(grid);

    auto munge = GaugeDoubleStoredMunger<DoubleStoredColourMatrixD, DoubleStoredColourMatrix>();

    Coordinate ldim = grid->LocalDimensions();
    thread_for(idx_g, grid->lSites(), {
        Coordinate coor;
        grid->LocalIndexToLocalCoor(idx_g, coor);

        bool isOdd = grid_rb->CheckerBoard(coor) == Odd;

        if(!isOdd) continue;

        int idx_o = (coor[Tdir] * ldim[Xdir] * ldim[Ydir] * ldim[Zdir]
                  +  coor[Xdir] * ldim[Ydir] * ldim[Zdir]
                  +  coor[Ydir] * ldim[Zdir]
                  +  coor[Zdir])/2;

        munge(iodata[idx_o], scalardata[idx_g]);
    });

    grid->Barrier(); timer.Stop();
    std::cout << Grid::GridLogMessage << "OpenQcdIO::readConfiguration: munge overhead " << timer.Elapsed() << std::endl;

    timer.Reset(); timer.Start();

    vectorizeFromLexOrdArray(scalardata, Umu_ds);

    grid->Barrier(); timer.Stop();
    std::cout << Grid::GridLogMessage << "OpenQcdIO::readConfiguration: vectorize overhead " << timer.Elapsed() << std::endl;

    timer.Reset(); timer.Start();

    undoDoubleStore(Umu, Umu_ds);

    grid->Barrier(); timer.Stop();
    std::cout << Grid::GridLogMessage << "OpenQcdIO::readConfiguration: redistribute overhead " << timer.Elapsed() << std::endl;

    PeriodicGaugeStatistics Stats; Stats(Umu, clone);

    RealD plaq_diff = fabs(clone.plaquette - header.plaquette);

    // clang-format off
    std::cout << GridLogMessage << "OpenQcd Configuration " << file
              << " plaquette " << clone.plaquette
              << " header " << header.plaquette
              << " difference " << plaq_diff
              << std::endl;
    // clang-format on

    RealD precTol = (getPrecision<vsimd>::value == 1) ? 2e-7 : 2e-15;
    RealD tol     = precTol * std::sqrt(grid->_Nprocessors); // taken from RQCD chroma code

    if(plaq_diff >= tol)
      std::cout << " Plaquette mismatch (diff = " << plaq_diff << ", tol = " << tol << ")" << std::endl;
    assert(plaq_diff < tol);

    std::cout << GridLogMessage << "OpenQcd Configuration " << file << " and plaquette agree" << std::endl;
  }

  template<class vsimd>
  static inline void writeConfiguration(Lattice<iLorentzColourMatrix<vsimd>>& Umu,
                                        std::string                           file) {
    std::cout << GridLogError << "Writing to openQCD file format is not implemented" << std::endl;
    exit(EXIT_FAILURE);
  }

private:
  static inline GridCartesian* createOpenQcdGrid(GridCartesian* grid) {
    // exploit GridCartesian to be able to still use IOobject
    Coordinate gdim  = grid->GlobalDimensions();
    Coordinate ldim  = grid->LocalDimensions();
    Coordinate pcoor = grid->ThisProcessorCoor();

    // openqcd does rb on the z direction
    gdim[Zdir] /= 2;
    ldim[Zdir] /= 2;

    // and has the order T X Y Z (from slowest to fastest)
    std::swap(gdim[Xdir], gdim[Zdir]);
    std::swap(ldim[Xdir], ldim[Zdir]);
    std::swap(pcoor[Xdir], pcoor[Zdir]);

    GridCartesian* ret   = SpaceTimeGrid::makeFourDimGrid(gdim, grid->_simd_layout, grid->ProcessorGrid());
    ret->_ldimensions    = ldim;
    ret->_processor_coor = pcoor;
    return ret;
  }

  template<class vsimd>
  static inline void undoDoubleStore(Lattice<iLorentzColourMatrix<vsimd>>&            Umu,
                                     Lattice<iDoubleStoredColourMatrix<vsimd>> const& Umu_ds) {
    conformable(Umu.Grid(), Umu_ds.Grid());
    Lattice<iColourMatrix<vsimd>> U(Umu.Grid());

    // they store T+, T-, X+, X-, Y+, Y-, Z+, Z-
    for(int mu_g = 0; mu_g < Nd; ++mu_g) {
      int mu_o = (mu_g + 1) % Nd;
      U        = PeekIndex<LorentzIndex>(Umu_ds, 2 * mu_o)
               + Cshift(PeekIndex<LorentzIndex>(Umu_ds, 2 * mu_o + 1), mu_g, +1);
      PokeIndex<LorentzIndex>(Umu, U, mu_g);
    }
  }
};

NAMESPACE_END(Grid);
