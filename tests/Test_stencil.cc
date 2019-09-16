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

using namespace std;
using namespace Grid;
 ;

int main(int argc, char ** argv) {
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

  typedef CartesianStencil<vobj,vobj,int> Stencil;
    for(int dir=0;dir<4;dir++){
      for(int disp=0;disp<Fine._fdimensions[dir];disp++){

	std::cout<< std::fixed <<GridLogMessage << "Using stencil to shift dim "<<dir<< " by "<<disp<<std::endl;
	// start to test the Cartesian npoint stencil infrastructure
	int npoint=1;
	std::vector<int> directions(npoint,dir);
	std::vector<int> displacements(npoint,disp);

	Stencil myStencil(&Fine,npoint,0,directions,displacements,0);
	Coordinate ocoor(4);
	for(int o=0;o<Fine.oSites();o++){
	  Fine.oCoorFromOindex(ocoor,o);
	  ocoor[dir]=(ocoor[dir]+disp)%Fine._rdimensions[dir];
	}

	SimpleCompressor<vobj> compress;
	myStencil.HaloExchange(Foo,compress);

	Bar = Cshift(Foo,dir,disp);

	// Implement a stencil code that should agree with cshift!
	for(int i=0;i<Check.Grid()->oSites();i++){
	  
	  int permute_type;
	  StencilEntry *SE;
	  SE = myStencil.GetEntry(permute_type,0,i);
	  
	  auto check = Check.View();
	  auto foo   = Foo.View();
	  if ( SE->_is_local && SE->_permute )
	    permute(check[i],foo[SE->_offset],permute_type);
	  else if (SE->_is_local)
	    check[i] = foo[SE->_offset];
	  else { 
	    check[i] = myStencil.CommBuf()[SE->_offset];
	    //	    std::cout << " receive "<<i<<" " << Check[i]<<std::endl;
	    //	    std::cout << " Foo     "<<i<<" " <<   Foo[i]<<std::endl;
	  }
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
	  auto check = Check.View();
	  auto bar   = Bar.View();
	  for(int i=0;i<check.size();i++){
	    std::cout << i<<" Check "<<check[i]<< "\n"<<i<<" Bar "<<bar[i]<<std::endl;
	  }
	}
	if (nrm > 1.0e-4) exit(-1);

      }
    }

    std::cout<<GridLogMessage<<"Testing RedBlack\n ";


  Field EFoo(&rbFine);
  Field OFoo(&rbFine);
  Field ECheck(&rbFine);
  Field OCheck(&rbFine);
  pickCheckerboard(Even,EFoo,Foo);
  pickCheckerboard(Odd ,OFoo,Foo);

    for(int dir=0;dir<4;dir++){
      for(int disp=0;disp<rbFine._fdimensions[dir];disp++){

	std::cout<<GridLogMessage << "Using stencil to shift rb dim "<<dir<< " by "<<disp<<std::endl;
	// start to test the Cartesian npoint stencil infrastructure
	int npoint=1;
	std::vector<int> directions(npoint,dir);
	std::vector<int> displacements(npoint,disp);

	Stencil EStencil(&rbFine,npoint,Even,directions,displacements,0);
	Stencil OStencil(&rbFine,npoint,Odd,directions,displacements,0);

	Coordinate ocoor(4);
	for(int o=0;o<Fine.oSites();o++){
	  Fine.oCoorFromOindex(ocoor,o);
	  ocoor[dir]=(ocoor[dir]+disp)%Fine._rdimensions[dir];
	}

	SimpleCompressor<vobj> compress;

	Bar = Cshift(Foo,dir,disp);

	if ( disp & 0x1 ) {
	  ECheck.Checkerboard() = Even;
	  OCheck.Checkerboard() = Odd;
	} else {
	  ECheck.Checkerboard() = Odd;
	  OCheck.Checkerboard() = Even;
	}

	// Implement a stencil code that should agree with that darn cshift!
	EStencil.HaloExchange(EFoo,compress);
	for(int i=0;i<OCheck.Grid()->oSites();i++){
	  int permute_type;
	  StencilEntry *SE;
	  SE = EStencil.GetEntry(permute_type,0,i);
	  //	  std::cout << "Even source "<< i<<" -> " <<SE->_offset << " "<< SE->_is_local<<std::endl;

	  auto ocheck = OCheck.View();
	  auto efoo   = EFoo.View();
	  if ( SE->_is_local && SE->_permute )
	    permute(ocheck[i],efoo[SE->_offset],permute_type);
	  else if (SE->_is_local)
	    ocheck[i] = efoo[SE->_offset];
	  else
	    ocheck[i] = EStencil.CommBuf()[SE->_offset];
	}
	OStencil.HaloExchange(OFoo,compress);
	for(int i=0;i<ECheck.Grid()->oSites();i++){
	  int permute_type;
	  StencilEntry *SE;
	  SE = OStencil.GetEntry(permute_type,0,i);
	  //	  std::cout << "ODD source "<< i<<" -> " <<SE->_offset << " "<< SE->_is_local<<std::endl;

	  auto echeck = ECheck.View();
	  auto ofoo   = OFoo.View();
	  if ( SE->_is_local && SE->_permute )
	    permute(echeck[i],ofoo[SE->_offset],permute_type);
	  else if (SE->_is_local)
	    echeck[i] = ofoo[SE->_offset];
	  else
	    echeck[i] = OStencil.CommBuf()[SE->_offset];
	}

	setCheckerboard(Check,ECheck);
	setCheckerboard(Check,OCheck);

	Real nrmC = norm2(Check);
	Real nrmB = norm2(Bar);
	Diff = Check-Bar;
	Real nrm  = norm2(Diff);
	std::cout<<GridLogMessage<<"RB N2diff ="<<nrm<<" "<<nrmC<<" " <<nrmB<<std::endl;

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
	    std::cout <<"Coor (" << coor[0]<<","<<coor[1]<<","<<coor[2]<<","<<coor[3] <<") "
		      <<"shift "<<disp<<" dir "<< dir
		      << "  stencil impl " <<check<<" vs cshift impl "<<bar<<std::endl;
	  }

	}}}}

	if (nrm > 1.0e-4) exit(-1);

      }
    }

 Grid_finalize();
}
