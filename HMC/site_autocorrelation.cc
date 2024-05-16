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

template <class T> void readFile(T& out, std::string const fname){
  Grid::emptyUserRecord record;
  Grid::ScidacReader RD;
  RD.open(fname);
  RD.readScidacFieldRecord(out,record);
  RD.close();
}


int main(int argc, char **argv) {
  using namespace Grid;
  
  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size, simd_layout, mpi_layout);

  LatticeComplexD plaq1(&Grid), plaq2(&Grid);

  FieldMetaData header;

  double vol = plaq1.Grid()->gSites();
  
  std::string file1(argv[1]);
  std::cout << "Reading "<<file1<<std::endl;
  readFile(plaq1,file1);
  std::string file2(argv[2]);
  std::cout << "Reading "<<file2<<std::endl;
  readFile(plaq2,file2);
  
  auto p1bar = TensorRemove(sum(plaq1));
  auto p2bar = TensorRemove(sum(plaq2));

  p1bar = p1bar / vol;
  p2bar = p2bar / vol;

  std::cout<< GridLogMessage << "p1bar = "<<p1bar<<std::endl;
  std::cout<< GridLogMessage << "p2bar = "<<p2bar<<std::endl;

  auto corr_site = plaq1 * plaq2 - p1bar * p2bar;
  auto corr_bar  = TensorRemove(sum(corr_site))/vol;

  auto cov1_site = plaq1 * plaq1 - p1bar * p1bar;
  auto cov1_bar  = TensorRemove(sum(cov1_site))/vol;

  auto cov2_site = plaq2 * plaq2 - p2bar * p2bar;
  auto cov2_bar  = TensorRemove(sum(cov2_site))/vol;

  std::cout<< GridLogMessage << "cov_bar = "<<corr_bar<<std::endl;

  std::cout<< GridLogMessage << "corr_bar = "<<corr_bar/sqrt(cov1_bar*cov2_bar)<<std::endl;
  
  Grid_finalize();
}  // main


