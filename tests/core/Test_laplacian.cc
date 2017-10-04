    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/Test_laplacian.cc

    Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));


  std::vector<int> point({0,0,0,0});

  LatticeFermion src   (&Grid); //random(pRNG,src);
  SpinColourVectorD Sp;
  for (unsigned int s = 0; s < Ns; ++s)
      for (unsigned int c = 0; c < Nc; ++c)
        Sp()(s)(c) = 1;

  src = zero;
  pokeSite(Sp,src,point);

  LatticeFermion result(&Grid); result=zero;
  LatticeFermion tmp(&Grid); tmp=zero;

  // Gauge configuration
  LatticeGaugeField Umu(&Grid); SU3::HotConfiguration(pRNG,Umu);

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing the laplacian operator on a point source            "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  Laplacian<WilsonImplR> LaplaceOperator(src._grid);
  LaplaceOperator.ImportGauge(Umu);
  LaplaceOperator.M(src, result);

  std::cout << "Source vector" << std::endl;
  std::cout << src << std::endl;

  std::cout << "Result vector" << std::endl;
  std::cout << result << std::endl;

  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing the laplacian smearing operator on a point source   "<<std::endl;
  std::cout<<GridLogMessage<<"=============================================================="<<std::endl;

  LatticeFermion smeared  (&Grid); smeared = src;
  for (int smr = 0; smr < 10; ++smr)
  {
      LaplaceOperator.M(smeared, tmp);
      smeared += 0.1*tmp;
  }

  std::cout << "Smeared vector" << std::endl;
  std::cout << smeared << std::endl;

  // Norm of vector
  LatticeComplex smr_norm(&Grid);
  smr_norm = localNorm2(smeared);
  std::cout << "Smeared vector norm" << std::endl;
  std::cout << smr_norm << std::endl;



  Grid_finalize();
}
