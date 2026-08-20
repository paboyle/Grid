/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_io.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

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

//
// Interchange between vectorised and lexicographic (Nsimd()==1) lattices.
//
//   - transfer() in both directions via an unvectorise/vectorise pair
//   - RNG seeded identically on both layouts produces identical fields
//   - binary I/O written from one layout and read into the other
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-10;

////////////////////////////////////////////////////////////////////////
// Layout interchange.  Both lattices share the same scalar_object, so a
// single lexicographic array mediates; works in either direction.
////////////////////////////////////////////////////////////////////////
template<class vobjOut,class vobjIn>
void transfer(Lattice<vobjOut> &out,const Lattice<vobjIn> &in)
{
  typedef typename vobjIn::scalar_object sobj;
  static_assert(std::is_same<sobj,typename vobjOut::scalar_object>::value,
		"transfer: lattices must share scalar_object");

  GRID_ASSERT(out.Grid()->gSites() == in.Grid()->gSites());

  std::vector<sobj> buf;
  unvectorizeToLexOrdArray(buf,in);
  vectorizeFromLexOrdArray(buf,out);
}

////////////////////////////////////////////////////////////////////////
// Binary write/read of a single lattice object
////////////////////////////////////////////////////////////////////////
template<class Field>
void writeField(const Field &f,std::string file)
{
  typedef typename Field::vector_object vobj;
  typedef typename Field::scalar_object sobj;

  BinarySimpleMunger<sobj,sobj> munge;
  std::string format = getFormatString<vobj>();
  uint64_t offset = 0;
  uint32_t nersc_csum,scidac_csuma,scidac_csumb;

  BinaryIO::writeLatticeObject<vobj,sobj>(const_cast<Field &>(f),file,munge,offset,format,
					  nersc_csum,scidac_csuma,scidac_csumb);
}

template<class Field>
void readField(Field &f,std::string file)
{
  typedef typename Field::vector_object vobj;
  typedef typename Field::scalar_object sobj;

  BinarySimpleMunger<sobj,sobj> munge;
  std::string format = getFormatString<vobj>();
  uint64_t offset = 0;
  uint32_t nersc_csum,scidac_csuma,scidac_csumb;

  BinaryIO::readLatticeObject<vobj,sobj>(f,file,munge,offset,format,
					 nersc_csum,scidac_csuma,scidac_csumb);
}

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  Coordinate latt = GridDefaultLatt();
  Coordinate mpi  = GridDefaultMpi();
  Coordinate vsimd = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate ssimd({1,1,1,1});

  GridCartesian vGrid(latt,vsimd,mpi);
  GridCartesian sGrid(latt,ssimd,mpi);

  std::cout << GridLogMessage << "vectorised   grid Nsimd = " << vGrid.Nsimd() << std::endl;
  std::cout << GridLogMessage << "lexicographic grid Nsimd = " << sGrid.Nsimd() << std::endl;
  GRID_ASSERT( sGrid.Nsimd() == 1 );

  typedef LatticeComplexD    vField;
  typedef lexLatticeComplexD sField;

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG vRNG(&vGrid); vRNG.SeedFixedIntegers(seeds);
  GridParallelRNG sRNG(&sGrid); sRNG.SeedFixedIntegers(seeds);

  vField v(&vGrid);  random(vRNG,v);
  sField s(&sGrid);  random(sRNG,s);

  RealD nv = norm2(v);
  RealD ns = norm2(s);
  std::cout << GridLogMessage << "norm2 vectorised   = " << nv << std::endl;
  std::cout << GridLogMessage << "norm2 lexicographic = " << ns << std::endl;
  GRID_ASSERT( nv > tol );

  ////////////////////////////////////////////////////////
  // A: transfer round trip v -> s -> v
  ////////////////////////////////////////////////////////
  {
    sField st(&sGrid);
    vField vt(&vGrid);

    transfer(st,v);
    transfer(vt,st);

    vField d(&vGrid); d = vt - v;
    RealD n = norm2(d);
    std::cout << GridLogMessage << "A: round trip v->s->v  norm2(diff) = " << n << std::endl;
    GRID_ASSERT( n < tol );

    // and the transferred copy must carry the same norm
    std::cout << GridLogMessage << "A: norm2 transferred = " << norm2(st) << std::endl;
    GRID_ASSERT( fabs(norm2(st) - nv) < tol*nv );
  }

  ////////////////////////////////////////////////////////
  // B: identically seeded RNGs agree across layouts
  ////////////////////////////////////////////////////////
  {
    vField vs(&vGrid);
    transfer(vs,s);

    vField d(&vGrid); d = vs - v;
    RealD n = norm2(d);
    std::cout << GridLogMessage << "B: same seed, both layouts  norm2(diff) = " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  ////////////////////////////////////////////////////////
  // C: write vectorised, read lexicographic
  ////////////////////////////////////////////////////////
  {
    sField sref(&sGrid);  transfer(sref,v);   // what the file should contain
    sField sin(&sGrid);   sin = Zero();

    writeField(v,"nonsimd_io_v.bin");
    readField(sin,"nonsimd_io_v.bin");

    sField d(&sGrid); d = sin - sref;
    RealD n = norm2(d);
    std::cout << GridLogMessage << "C: write simd / read lexicographic  norm2(diff) = " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  ////////////////////////////////////////////////////////
  // D: write lexicographic, read vectorised
  ////////////////////////////////////////////////////////
  {
    sField sout(&sGrid);  transfer(sout,v);
    vField vin(&vGrid);   vin = Zero();

    writeField(sout,"nonsimd_io_s.bin");
    readField(vin,"nonsimd_io_s.bin");

    vField d(&vGrid); d = vin - v;
    RealD n = norm2(d);
    std::cout << GridLogMessage << "D: write lexicographic / read simd  norm2(diff) = " << n << std::endl;
    GRID_ASSERT( n < tol );
  }

  std::cout << GridLogMessage << "Test_io: ALL PASS" << std::endl;

  Grid_finalize();
}
