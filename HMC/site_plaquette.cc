/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: 

Copyright (C) 2017

Author: Peter Boyle

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
#include <Grid/Grid.h>
#include <string>

NAMESPACE_BEGIN(Grid);
template <class T> void writeFile(T& out, std::string const fname){
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(out.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(out,record,0,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
  WR.close();
#endif
}
NAMESPACE_END(Grid);
int main(int argc, char **argv) {
  using namespace Grid;
  
  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size, simd_layout, mpi_layout);

  LatticeGaugeField Umu(&Grid);
  std::vector<LatticeColourMatrix> U(4,&Grid);
  LatticeComplexD plaq(&Grid);

  FieldMetaData header;

  double vol = Umu.Grid()->gSites();
  double faces = (1.0 * Nd * (Nd - 1)) / 2.0;
  double Ncdiv = 1.0/Nc;
  
  std::string file1(argv[1]);
  std::string file2(argv[2]);
  std::cout << "Reading "<<file1<<std::endl;
  NerscIO::readConfiguration(Umu,header,file1);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  SU3WilsonLoops::sitePlaquette(plaq,U);

  plaq = plaq *(Ncdiv/faces);
  
  std::cout << "Writing "<<file2<<std::endl;
  writeFile(plaq,file2);
  
  Grid_finalize();
}  // main


