/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/parallelIO/OpenQcdIOChromaReference.h

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

#include <ios>
#include <iostream>
#include <limits>
#include <iomanip>
#include <mpi.h>
#include <ostream>
#include <string>

#define CHECK {std::cerr << __FILE__ << " @l " << __LINE__ << ": CHECK" << grid->ThisRank() << std::endl;}
#define CHECK_VAR(a)   { std::cerr << __FILE__ << "@l" << __LINE__ << " on "<< grid->ThisRank() << ": " << __func__ << " " << #a << "=" << (a) << std::endl; }
// #undef CHECK
// #define CHECK

NAMESPACE_BEGIN(Grid);

class ParRdr {
private:
  bool const swap;

  MPI_Status status;
  MPI_File   fp;

  int err;

  MPI_Datatype oddSiteType;
  MPI_Datatype fileViewType;

  GridBase* grid;

public:
  ParRdr(MPI_Comm comm, std::string const& filename, GridBase* gridPtr)
    : swap(false)
    , grid(gridPtr) {
    err = MPI_File_open(comm, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
    assert(err == MPI_SUCCESS);
  }

  virtual ~ParRdr() { MPI_File_close(&fp); }

  inline void errInfo(int const err, std::string const& func) {
    static char estring[MPI_MAX_ERROR_STRING];
    int         eclass = -1, len = 0;
    MPI_Error_class(err, &eclass);
    MPI_Error_string(err, estring, &len);
    std::cerr << func << " - Error " << eclass << ": " << estring << std::endl;
  }

  int readHeader(FieldMetaData& field) {
    assert((grid->_ndimension == Nd) && (Nd == 4));
    assert(Nc == 3);

    OpenQcdHeader header;

    readBlock(reinterpret_cast<char*>(&header), 0, sizeof(OpenQcdHeader), MPI_CHAR);

    header.plaq /= 3.; // TODO change this into normalizationfactor

    // sanity check (should trigger on endian issues) TODO remove?
    assert(0 < header.Nt && header.Nt <= 1024);
    assert(0 < header.Nx && header.Nx <= 1024);
    assert(0 < header.Ny && header.Ny <= 1024);
    assert(0 < header.Nz && header.Nz <= 1024);

    field.dimension[0] = header.Nx;
    field.dimension[1] = header.Ny;
    field.dimension[2] = header.Nz;
    field.dimension[3] = header.Nt;

    for(int d = 0; d < Nd; d++)
      assert(grid->FullDimensions()[d] == field.dimension[d]);

    field.plaquette = header.plaq;

    field.data_start = sizeof(OpenQcdHeader);

    return field.data_start;
  }

  void readBlock(void* const dest, uint64_t const pos, uint64_t const nbytes, MPI_Datatype const datatype) {
    err = MPI_File_read_at_all(fp, pos, dest, nbytes, datatype, &status);
    errInfo(err, "MPI_File_read_at_all");
    // CHECK_VAR(err)

    int read = -1;
    MPI_Get_count(&status, datatype, &read);
    // CHECK_VAR(read)
    assert(nbytes == (uint64_t)read);
    assert(err == MPI_SUCCESS);
  }

  void createTypes() {
    constexpr int elem_size = Nd * 2 * 2 * Nc * Nc * sizeof(double); // 2_complex 2_fwdbwd

    err = MPI_Type_contiguous(elem_size, MPI_BYTE, &oddSiteType); assert(err == MPI_SUCCESS);
    err = MPI_Type_commit(&oddSiteType); assert(err == MPI_SUCCESS);

    Coordinate const L = grid->GlobalDimensions();
    Coordinate const l = grid->LocalDimensions();
    Coordinate const i = grid->ThisProcessorCoor();

    Coordinate sizes({L[2] / 2, L[1], L[0], L[3]});
    Coordinate subsizes({l[2] / 2, l[1], l[0], l[3]});
    Coordinate starts({i[2] * l[2] / 2, i[1] * l[1], i[0] * l[0], i[3] * l[3]});

    err = MPI_Type_create_subarray(grid->_ndimension, &sizes[0], &subsizes[0], &starts[0], MPI_ORDER_FORTRAN, oddSiteType, &fileViewType); assert(err == MPI_SUCCESS);
    err = MPI_Type_commit(&fileViewType); assert(err == MPI_SUCCESS);
  }

  void freeTypes() {
    err = MPI_Type_free(&fileViewType); assert(err == MPI_SUCCESS);
    err = MPI_Type_free(&oddSiteType); assert(err == MPI_SUCCESS);
  }

  bool readGauge(std::vector<ColourMatrixD>& domain_buff, FieldMetaData& meta) {
    auto hdr_offset = readHeader(meta);
    CHECK
    createTypes();
    err = MPI_File_set_view(fp, hdr_offset, oddSiteType, fileViewType, "native", MPI_INFO_NULL); errInfo(err, "MPI_File_set_view0"); assert(err == MPI_SUCCESS);
    CHECK
    int const domainSites = grid->lSites();
    domain_buff.resize(Nd * domainSites); // 2_fwdbwd * 4_Nd * domainSites / 2_onlyodd

    // the actual READ
    constexpr uint64_t cm_size   = 2 * Nc * Nc * sizeof(double);    // 2_complex
    constexpr uint64_t os_size   = Nd * 2 * cm_size;                // 2_fwdbwd
    constexpr uint64_t max_elems = std::numeric_limits<int>::max(); // int adressable elems: floor is fine
    uint64_t const     n_os      = domainSites / 2;

    for(uint64_t os_idx = 0; os_idx < n_os;) {
      uint64_t const read_os = os_idx + max_elems <= n_os ? max_elems : n_os - os_idx;
      uint64_t const cm      = os_idx * Nd * 2;
      readBlock(&(domain_buff[cm]), os_idx, read_os, oddSiteType);
      os_idx += read_os;
    }

    CHECK
    err = MPI_File_set_view(fp, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
  errInfo(err, "MPI_File_set_view1");
    assert(err == MPI_SUCCESS);
    freeTypes();

    std::cout << GridLogMessage << "read sum: " << n_os * os_size << " bytes" << std::endl;
    return true;
  }
};

class OpenQcdIOChromaReference : public BinaryIO {
public:
  template<class vsimd>
  static inline void readConfiguration(Lattice<iLorentzColourMatrix<vsimd>>& Umu,
                                       Grid::FieldMetaData&                  header,
                                       std::string                           file) {
    typedef Lattice<iDoubleStoredColourMatrix<vsimd>> DoubledGaugeField;

    assert(Ns == 4 and Nd == 4 and Nc == 3);

    auto grid = Umu.Grid();

    typedef ColourMatrixD fobj;

    std::vector<fobj> iodata(
      Nd * grid->lSites()); // actual size = 2*Nd*lsites but have only lsites/2 sites in file

    {
      ParRdr rdr(MPI_COMM_WORLD, file, grid);
      rdr.readGauge(iodata, header);
    } // equivalent to using binaryio

    std::vector<iDoubleStoredColourMatrix<typename vsimd::scalar_type>> Umu_ds_scalar(grid->lSites());

    copyToLatticeObject(Umu_ds_scalar, iodata, grid); // equivalent to munging

    DoubledGaugeField Umu_ds(grid);

    vectorizeFromLexOrdArray(Umu_ds_scalar, Umu_ds);

    redistribute(Umu, Umu_ds); // equivalent to undoDoublestore

    FieldMetaData clone(header);

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

private:
  template<class vsimd>
  static inline void redistribute(Lattice<iLorentzColourMatrix<vsimd>>&            Umu,
                                  Lattice<iDoubleStoredColourMatrix<vsimd>> const& Umu_ds) {
    Grid::conformable(Umu.Grid(), Umu_ds.Grid());
    Lattice<iColourMatrix<vsimd>> U(Umu.Grid());

    U = PeekIndex<LorentzIndex>(Umu_ds, 2) + Cshift(PeekIndex<LorentzIndex>(Umu_ds, 3), 0, +1); PokeIndex<LorentzIndex>(Umu, U, 0);
    U = PeekIndex<LorentzIndex>(Umu_ds, 4) + Cshift(PeekIndex<LorentzIndex>(Umu_ds, 5), 1, +1); PokeIndex<LorentzIndex>(Umu, U, 1);
    U = PeekIndex<LorentzIndex>(Umu_ds, 6) + Cshift(PeekIndex<LorentzIndex>(Umu_ds, 7), 2, +1); PokeIndex<LorentzIndex>(Umu, U, 2);
    U = PeekIndex<LorentzIndex>(Umu_ds, 0) + Cshift(PeekIndex<LorentzIndex>(Umu_ds, 1), 3, +1); PokeIndex<LorentzIndex>(Umu, U, 3);
  }

  static inline void copyToLatticeObject(std::vector<DoubleStoredColourMatrix>& u_fb,
                                         std::vector<ColourMatrixD> const&      node_buff,
                                         GridBase*                              grid) {
    assert(node_buff.size() == Nd * grid->lSites());

    Coordinate const& l = grid->LocalDimensions();

    Coordinate coord(Nd);
    int&       x = coord[0];
    int&       y = coord[1];
    int&       z = coord[2];
    int&       t = coord[3];

    int buff_idx = 0;
    for(t = 0; t < l[3]; ++t) // IMPORTANT: openQCD file ordering
      for(x = 0; x < l[0]; ++x)
        for(y = 0; y < l[1]; ++y)
          for(z = 0; z < l[2]; ++z) {
            if((t + z + y + x) % 2 == 0) continue;

            int local_idx;
            Lexicographic::IndexFromCoor(coord, local_idx, grid->LocalDimensions());
            for(int mu = 0; mu < 2 * Nd; ++mu)
              for(int c1 = 0; c1 < Nc; ++c1) {
                for(int c2 = 0; c2 < Nc; ++c2) {
                  u_fb[local_idx](mu)()(c1,c2) = node_buff[mu+buff_idx]()()(c1,c2);
                }
              }
            buff_idx += 2 * Nd;
          }

    assert(node_buff.size() == buff_idx);
  }
};

NAMESPACE_END(Grid);
