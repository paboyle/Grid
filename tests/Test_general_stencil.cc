    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/Test_stencil.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/stencil/GeneralLocalStencil.h>

using namespace std;
using namespace Grid;

int main(int argc, char ** argv) 
{
  Grid_init(&argc, &argv);

  //  typedef LatticeColourMatrix Field;
  typedef LatticeComplex Field;
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object sobj;

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();

  double volume = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

  GridCartesian Fine(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian rbFine(&Fine);
  GridParallelRNG       fRNG(&Fine);

  //  fRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9});
  std::vector<int> seeds({1,2,3,4});
  fRNG.SeedFixedIntegers(seeds);

  Field Foo(&Fine);
  Field Bar(&Fine);
  Field Check(&Fine);
  Field Diff(&Fine);
  LatticeComplex lex(&Fine);

  lex = Zero();
  random(fRNG,Foo);
  gaussian(fRNG,Bar);

  for (int i=0;i<simd_layout.size();i++){
    std::cout <<" simd layout "<<i<<" = "<<simd_layout[i]<<std::endl;
  }
  Integer stride =1000;
  {
    LatticeComplex coor(&Fine);

    for(int d=0;d<Nd;d++){
      LatticeCoordinate(coor,d);
      lex = lex + coor*stride;
      stride=stride/10;
    }
    Foo=lex;
  }

  typedef GeneralLocalStencil        GeneralStencil;

    for(int dir1=0;dir1<4;dir1++){
    for(int dir2=0;dir2<4;dir2++){
      if ( dir2!= dir1){
      for(int disp1=0;disp1<Fine._fdimensions[dir1];disp1++){
      for(int disp2=0;disp2<Fine._fdimensions[dir2];disp2++){

	std::cout<< std::fixed <<GridLogMessage << "Using stencil to shift dim "<<dir1<< " by "<<disp1<<std::endl;
	std::cout<< std::fixed <<GridLogMessage << "                   and dim "<<dir2<< " by "<<disp2<<std::endl;
	// start to test the Cartesian npoint stencil infrastructure
	int npoint=1;

	Coordinate shift(Nd); 
	for(int d=0;d<Nd;d++) shift[d]=0;
	shift[dir1]=disp1;
	shift[dir2]=disp2;
	std::vector<Coordinate> shifts(npoint,shift);
	GeneralLocalStencil gStencil(&Fine,shifts);

	Bar = Cshift(Foo,dir1,disp1);
	Bar = Cshift(Bar,dir2,disp2);

	// Implement a stencil code that should agree with cshift!
	for(int i=0;i<Check.Grid()->oSites();i++){
	  auto SE = gStencil.GetEntry(0,i);
	  autoView(check, Check, CpuWrite);
	  autoView(  foo, Foo, CpuRead);
	  // Encapsulate in a general wrapper
	  check[i] = foo[SE->_offset];                   auto tmp=check[i];
	  if (SE->_permute & 0x1 ) { permute(check[i],tmp,0); tmp=check[i];}
	  if (SE->_permute & 0x2 ) { permute(check[i],tmp,1); tmp=check[i];}
	  if (SE->_permute & 0x4 ) { permute(check[i],tmp,2); tmp=check[i];}
	  if (SE->_permute & 0x8 ) { permute(check[i],tmp,3); tmp=check[i];}
	}

	Real nrmC = norm2(Check);
	Real nrmB = norm2(Bar);
	Diff = Check-Bar;
	Real nrm  = norm2(Diff);
	std::cout<<GridLogMessage<<"N2diff ="<<nrm<<" "<<nrmC<<" " <<nrmB<<std::endl;

	Coordinate coor(4);
	for(coor[3]=0;coor[3]<latt_size[3]/mpi_layout[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2]/mpi_layout[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1]/mpi_layout[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0]/mpi_layout[0];coor[0]++){

	  RealD diff;
	  sobj check,bar;
	  peekSite(check,Check,coor);
	  peekSite(bar,Bar,coor);

	  sobj ddiff;
	  ddiff = check -bar;
	  diff =norm2(ddiff);
	  if ( diff > 0){
	    std::cout <<"Coor (" << coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3]
		      <<") " <<check<<" vs "<<bar<<std::endl;
	  }


	}}}}

	if (nrm > 1.0e-4) {
	  autoView( check , Check, CpuRead);
	  autoView(   bar ,   Bar, CpuRead);
	  for(int i=0;i<check.size();i++){
	    std::cout << i<<" Check "<<check[i]<< "\n"<<i<<" Bar "<<bar[i]<<std::endl;
	  }
	}
	if (nrm > 1.0e-4) exit(-1);

      }}
      }
    }}
 Grid_finalize();
}
