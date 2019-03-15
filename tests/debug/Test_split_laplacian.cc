    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_mrhs_cg.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  typedef LatticeComplex ComplexField; 

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  int nd   = latt_size.size();
  int ndm1 = nd-1;

  std::vector<int> simd_layout = GridDefaultSimd(nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> mpi_split (mpi_layout.size(),1);

  std::cout << " Full " << GridCmdVectorIntToString(latt_size)  << " subgrid"         <<std::endl;
  std::cout << " Full " << GridCmdVectorIntToString(mpi_layout) << " sub communicator"<<std::endl;
  std::cout << " Full " << GridCmdVectorIntToString(simd_layout)<< " simd layout "    <<std::endl;

  GridCartesian         * GridN = new GridCartesian(latt_size,
						    simd_layout,
						    mpi_layout);

  std::vector<int> latt_m  = latt_size;   latt_m[nd-1] = 1;
  std::vector<int> mpi_m   = mpi_layout;  mpi_m [nd-1] = 1;
  std::vector<int> simd_m  = GridDefaultSimd(ndm1,vComplex::Nsimd()); simd_m.push_back(1);


  std::cout << " Requesting " << GridCmdVectorIntToString(latt_m)<< " subgrid"         <<std::endl;
  std::cout << " Requesting " << GridCmdVectorIntToString(mpi_m) << " sub communicator"<<std::endl;
  std::cout << " Requesting " << GridCmdVectorIntToString(simd_m)<< " simd layout "    <<std::endl;
  GridCartesian         * Grid_m = new GridCartesian(latt_m,
						     simd_m,
						     mpi_m,
						     *GridN); 

  Complex C(1.0);
  Complex tmp;

  ComplexField Full(GridN); Full = C;
  ComplexField Full_cpy(GridN);
  ComplexField Split(Grid_m);Split= C;

  std::cout << GridLogMessage<< " Full  volume "<< norm2(Full) <<std::endl;
  std::cout << GridLogMessage<< " Split volume "<< norm2(Split) <<std::endl;

  tmp=C;
  GridN->GlobalSum(tmp);
  std::cout << GridLogMessage<< " Full  nodes "<< tmp <<std::endl;

  tmp=C;
  Grid_m->GlobalSum(tmp);
  std::cout << GridLogMessage<< " Split nodes "<< tmp <<std::endl;
  GridN->Barrier();

  auto local_latt = GridN->LocalDimensions();

  Full_cpy = zero;
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          RNG(GridN);  RNG.SeedFixedIntegers(seeds);

  random(RNG,Full);
  for(int t=0;t<local_latt[nd-1];t++){
    ExtractSliceLocal(Split,Full,0,t,Tp);
    InsertSliceLocal (Split,Full_cpy,0,t,Tp);
  }
  Full_cpy = Full_cpy - Full;
  std::cout << " NormFull " << norm2(Full)<<std::endl;
  std::cout << " NormDiff " << norm2(Full_cpy)<<std::endl;
  Grid_finalize();
}
