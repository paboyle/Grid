/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/parallelIO/IldgIO.h

Copyright (C) 2015

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
#ifndef GRID_ILDG_IO_H
#define GRID_ILDG_IO_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

#include <pwd.h>
#include <sys/utsname.h>
#include <unistd.h>

extern "C" {  // for linkage
#include "lime.h"
}

namespace Grid {
namespace QCD {

inline void ILDGGrid(GridBase *grid, ILDGField &header) {
  assert(grid->_ndimension == 4);  // emit error if not
  header.dimension.resize(4);
  header.boundary.resize(4);
  for (int d = 0; d < 4; d++) {
    header.dimension[d] = grid->_fdimensions[d];
    // Read boundary conditions from ... ?
    header.boundary[d] = std::string("periodic");
  }
}

inline void ILDGChecksum(uint32_t *buf, uint32_t buf_size_bytes,
                         uint32_t &csum) {
  BinaryIO::Uint32Checksum(buf, buf_size_bytes, csum);
}

//////////////////////////////////////////////////////////////////////
// Utilities ; these are QCD aware
//////////////////////////////////////////////////////////////////////
template <class GaugeField>
inline void ILDGStatistics(GaugeField &data, ILDGField &header) {
  // How to convert data precision etc...
  header.link_trace = Grid::QCD::WilsonLoops<PeriodicGimplR>::linkTrace(data);
  header.plaquette = Grid::QCD::WilsonLoops<PeriodicGimplR>::avgPlaquette(data);
  // header.polyakov =
}

// Forcing QCD here
template <class fobj, class sobj>
struct ILDGMunger {
  void operator()(fobj &in, sobj &out, uint32_t &csum) {
    for (int mu = 0; mu < 4; mu++) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          out(mu)()(i, j) = in(mu)()(i, j);
        }
      }
    }
    ILDGChecksum((uint32_t *)&in, sizeof(in), csum);
  };
};

template <class fobj, class sobj>
struct ILDGSimpleUnmunger {
  void operator()(sobj &in, fobj &out, uint32_t &csum) {
    for (int mu = 0; mu < 4; mu++) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          out(mu)()(i, j) = in(mu)()(i, j);
        }
      }
    }
    ILDGChecksum((uint32_t *)&out, sizeof(out), csum);
  };
};

////////////////////////////////////////////////////////////////////////////////
// Write and read from fstream; compute header offset for payload
////////////////////////////////////////////////////////////////////////////////
class ILDGIO : public BinaryIO {
  FILE *outFile;
  LimeWriter *LimeW;
  LimeRecordHeader *LimeHeader;

 public:
  ILDGIO(std::string file) {
    outFile = fopen(file.c_str(), "w");
    // check if opened correctly

    LimeW = limeCreateWriter(outFile);
  }

  ~ILDGIO() { fclose(outFile); }

  unsigned int writeHeader(ILDGField &header) {
    // write header in LIME
    n_uint64_t nbytes;
    int MB_flag = 1, ME_flag = 0;

    char message[] = "ildg-format";
    nbytes = strlen(message);
    LimeHeader = limeCreateHeader(MB_flag, ME_flag, message, nbytes);
    limeWriteRecordHeader(LimeHeader, LimeW);
    limeDestroyHeader(LimeHeader);
    // save the xml header here
    // use the xml_writer to c++ streams in pugixml
    // and convert to char message
    // limeWriteRecordData(message, &nbytes, LimeW);
    limeWriterCloseRecord(LimeW);

    return 0;
  }

  unsigned int readHeader(std::string file, GridBase *grid, ILDGField &field) {
    return 0;
  }

  template <class vsimd>
  int readConfiguration(Lattice<iLorentzColourMatrix<vsimd> > &Umu,
                        ILDGField &header, std::string file) {
    typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;

    return 0;
  }

  template <class vsimd>
  int writeConfiguration(Lattice<iLorentzColourMatrix<vsimd> > &Umu,
                         ILDGField &header, std::string file) {
    typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;
    typedef iLorentzColourMatrix<vsimd> vobj;
    typedef typename vobj::scalar_object sobj;
    typedef LorentzColourMatrixD fobj;

    ILDGSimpleUnmunger<fobj, sobj> munge;
    unsigned int offset = writeHeader(header);

    BinaryIO::Uint32Checksum<vobj, fobj>(Umu, munge, header.checksum);

    // Write record header 
    LimeRecordHeader *h;
    std::cout << GridLogDebug << "ILDG Creating Header" << std::endl;
    char message[] = "ildg-binary-data";
    h = limeCreateHeader(1, 1, message, strlen(message));

    std::cout << GridLogDebug << "ILDG Writing Header" << std::endl;
    int status = limeWriteRecordHeader(h, LimeW);

    if (status < 0) {
      std::cerr << "ILDG Header error\n";
      return 1;
    }

    limeDestroyHeader(h);

    ILDGtype ILDGt(true, LimeW);
    uint32_t csum = BinaryIO::writeObjectParallel<vobj, fobj>(
        Umu, file, munge, offset, header.floating_point, ILDGt);

    limeWriterCloseRecord(LimeW);

    return 0;
  }

  // format for RNG?
};
}
}
#endif
