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
struct ILDGUnmunger {
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
enum ILDGstate {ILDGread, ILDGwrite};

class ILDGIO : public BinaryIO {
  FILE *File;
  LimeWriter *LimeW;
  LimeRecordHeader *LimeHeader;
  LimeReader *LimeR;
  std::string filename;


 public:
  ILDGIO(std::string file, ILDGstate RW) {
      filename = file;
    if (RW == ILDGwrite){
      File = fopen(file.c_str(), "w");
      // check if opened correctly

      LimeW = limeCreateWriter(File);
    } else {
      File = fopen(file.c_str(), "r");
      // check if opened correctly

      LimeR = limeCreateReader(File);
    }
  }

  ~ILDGIO() { fclose(File); }

  int createHeader(std::string message, int MB, int ME, size_t PayloadSize, LimeWriter* L){
    LimeRecordHeader *h;
    h = limeCreateHeader(MB, ME, const_cast<char *>(message.c_str()), PayloadSize);
    int status = limeWriteRecordHeader(h, L);
    if (status < 0) {
      std::cerr << "ILDG Header error\n";
      return status;
    }
    limeDestroyHeader(h);
    return LIME_SUCCESS;
  }

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
    limeWriteRecordData(message, &nbytes, LimeW);
    limeWriterCloseRecord(LimeW);

    return 0;
  }

  unsigned int readHeader(ILDGField &header) {
    return 0;
  }

  template <class vsimd>
  uint32_t readConfiguration(Lattice<iLorentzColourMatrix<vsimd> > &Umu) {
    typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;
    typedef LorentzColourMatrixD sobjd;
    typedef LorentzColourMatrixF sobjf;
    typedef iLorentzColourMatrix<vsimd> itype;
    typedef LorentzColourMatrix sobj;
    GridBase *grid = Umu._grid;

    ILDGField header;
    readHeader(header);

    // now just the conf, ignore the header
    std::string format = std::string("IEEE64BIG");
    do {limeReaderNextRecord(LimeR);}
    while (strncmp(limeReaderType(LimeR), "ildg-binary-data",16));

    n_uint64_t nbytes = limeReaderBytes(LimeR);//size of this record (configuration)


    ILDGtype ILDGt(true, LimeR);
    // this is special for double prec data, just for the moment
    uint32_t csum = BinaryIO::readObjectParallel< itype, sobjd >(
       Umu, filename, ILDGMunger<sobjd, sobj>(), 0, format, ILDGt);

    // Check configuration 
    // todo

    return csum;
  }

  template <class vsimd>
  uint32_t writeConfiguration(Lattice<iLorentzColourMatrix<vsimd> > &Umu, std::string format) {
    typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;
    typedef iLorentzColourMatrix<vsimd> vobj;
    typedef typename vobj::scalar_object sobj;
    typedef LorentzColourMatrixD fobj;

    ILDGField header;
    // fill the header
    header.floating_point = format;

    ILDGUnmunger<fobj, sobj> munge;
    unsigned int offset = writeHeader(header);

    BinaryIO::Uint32Checksum<vobj, fobj>(Umu, munge, header.checksum);

    // Write data record header
    n_uint64_t PayloadSize = sizeof(fobj) * Umu._grid->_gsites;
    createHeader("ildg-binary-data", 0, 1, PayloadSize, LimeW);

    ILDGtype ILDGt(true, LimeW);
    uint32_t csum = BinaryIO::writeObjectParallel<vobj, fobj>(
       Umu, filename, munge, 0, header.floating_point, ILDGt);

    limeWriterCloseRecord(LimeW);

    // Last record
    // the logical file name LNF
    // look into documentation on how to generate this string
    std::string LNF = "empty"; 


    PayloadSize = sizeof(LNF);
    createHeader("ildg-binary-lfn", 1 , 1, PayloadSize, LimeW);
    limeWriteRecordData(const_cast<char*>(LNF.c_str()), &PayloadSize, LimeW);

    limeWriterCloseRecord(LimeW);

    return csum;
  }

  // format for RNG?
};
}
}
#endif
