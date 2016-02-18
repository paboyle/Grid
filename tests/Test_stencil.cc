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
#include "Grid.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  //  typedef LatticeColourMatrix Field;
  typedef LatticeComplex Field;
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object sobj;

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  double volume = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
    
  GridCartesian Fine(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian rbFine(latt_size,simd_layout,mpi_layout);
  GridParallelRNG       fRNG(&Fine);

  //  fRNG.SeedRandomDevice();
  std::vector<int> seeds({1,2,3,4});
  fRNG.SeedFixedIntegers(seeds);
  
  Field Foo(&Fine);
  Field Bar(&Fine);
  Field Check(&Fine);
  Field Diff(&Fine);
  LatticeComplex lex(&Fine);

  lex = zero;  
  random(fRNG,Foo);
  gaussian(fRNG,Bar);

  /*
  Integer stride =1000;
  {
    double nrm;
    LatticeComplex coor(&Fine);

    for(int d=0;d<Nd;d++){
      LatticeCoordinate(coor,d);
      lex = lex + coor*stride;
      stride=stride/10;
    }
    Foo=lex;
  }
  */

  typedef CartesianStencil<vobj,vobj> Stencil;
    for(int dir=0;dir<4;dir++){
      for(int disp=0;disp<Fine._fdimensions[dir];disp++){

	std::cout<< std::fixed <<GridLogMessage << "Using stencil to shift dim "<<dir<< " by "<<disp<<std::endl;
	// start to test the Cartesian npoint stencil infrastructure
	int npoint=1;
	std::vector<int> directions(npoint,dir);
	std::vector<int> displacements(npoint,disp);

	Stencil myStencil(&Fine,npoint,0,directions,displacements);

	std::vector<int> ocoor(4);
	for(int o=0;o<Fine.oSites();o++){
	  Fine.oCoorFromOindex(ocoor,o);
	  ocoor[dir]=(ocoor[dir]+disp)%Fine._rdimensions[dir];
	}
	
	SimpleCompressor<vobj> compress;
	myStencil.HaloExchange(Foo,compress);

	Bar = Cshift(Foo,dir,disp);

	// Implement a stencil code that should agree with cshift!
	for(int i=0;i<Check._grid->oSites();i++){
	  
	  int permute_type;
	  StencilEntry *SE;
	  SE = myStencil.GetEntry(permute_type,0,i);
	  
	  if ( SE->_is_local && SE->_permute )
	    permute(Check._odata[i],Foo._odata[SE->_offset],permute_type);
	  else if (SE->_is_local)
	    Check._odata[i] = Foo._odata[SE->_offset];
	  else 
	    Check._odata[i] = myStencil.comm_buf[SE->_offset];
	}

	Real nrmC = norm2(Check);
	Real nrmB = norm2(Bar);
	Diff = Check-Bar;
	Real nrm  = norm2(Diff);
	std::cout<<GridLogMessage<<"N2diff ="<<nrm<<" "<<nrmC<<" " <<nrmB<<std::endl;

	std::vector<int> coor(4);
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

	Stencil EStencil(&rbFine,npoint,Even,directions,displacements);
	Stencil OStencil(&rbFine,npoint,Odd,directions,displacements);

	std::vector<int> ocoor(4);
	for(int o=0;o<Fine.oSites();o++){
	  Fine.oCoorFromOindex(ocoor,o);
	  ocoor[dir]=(ocoor[dir]+disp)%Fine._rdimensions[dir];
	}
	
	SimpleCompressor<vobj> compress;

	EStencil.HaloExchange(EFoo,compress);
	OStencil.HaloExchange(OFoo,compress);
	
	Bar = Cshift(Foo,dir,disp);

	if ( disp & 0x1 ) {
	  ECheck.checkerboard = Even;
	  OCheck.checkerboard = Odd;
	} else { 
	  ECheck.checkerboard = Odd;
	  OCheck.checkerboard = Even;
	}

	// Implement a stencil code that should agree with that darn cshift!
	for(int i=0;i<OCheck._grid->oSites();i++){
	  int permute_type;
	  StencilEntry *SE;
	  SE = EStencil.GetEntry(permute_type,0,i);
	  //	  std::cout << "Even source "<< i<<" -> " <<SE->_offset << " "<< SE->_is_local<<std::endl;

	  if ( SE->_is_local && SE->_permute )
	    permute(OCheck._odata[i],EFoo._odata[SE->_offset],permute_type);
	  else if (SE->_is_local)
	    OCheck._odata[i] = EFoo._odata[SE->_offset];
	  else 
	    OCheck._odata[i] = EStencil.comm_buf[SE->_offset];
	}
	for(int i=0;i<ECheck._grid->oSites();i++){
	  int permute_type;
	  StencilEntry *SE;
	  SE = OStencil.GetEntry(permute_type,0,i);
	  //	  std::cout << "ODD source "<< i<<" -> " <<SE->_offset << " "<< SE->_is_local<<std::endl;
	  
	  if ( SE->_is_local && SE->_permute )
	    permute(ECheck._odata[i],OFoo._odata[SE->_offset],permute_type);
	  else if (SE->_is_local)
	    ECheck._odata[i] = OFoo._odata[SE->_offset];
	  else 
	    ECheck._odata[i] = OStencil.comm_buf[SE->_offset];
	}
	
	setCheckerboard(Check,ECheck);
	setCheckerboard(Check,OCheck);
	
	Real nrmC = norm2(Check);
	Real nrmB = norm2(Bar);
	Diff = Check-Bar;
	Real nrm  = norm2(Diff);
	std::cout<<GridLogMessage<<"RB N2diff ="<<nrm<<" "<<nrmC<<" " <<nrmB<<std::endl;

	std::vector<int> coor(4);
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


      }
    }

 Grid_finalize();
}
