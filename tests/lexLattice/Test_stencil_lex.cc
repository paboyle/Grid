/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_stencil_lex.cc

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
// CartesianStencil on a lexicographic (Nsimd()==1) lattice, checked against
// Cshift for every direction and displacement of either sign.
//
// On a unit simd grid no stencil entry may request a permute; the Grid_simd1
// permute is an assert, so a request would abort rather than pass silently.
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-10;

int main(int argc, char ** argv)
{
  Grid_init(&argc, &argv);

  typedef lexLatticeComplexD Field;
  typedef Field::vector_object vobj;

  Coordinate latt = GridDefaultLatt();
  Coordinate simd({1,1,1,1});
  Coordinate mpi  = GridDefaultMpi();

  GridCartesian Fine(latt,simd,mpi);

  std::cout << GridLogMessage << "Lexicographic grid, Nsimd = " << Fine.Nsimd() << std::endl;
  GRID_ASSERT( Fine.Nsimd() == 1 );

  GridParallelRNG fRNG(&Fine);
  fRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  Field Foo(&Fine);   random(fRNG,Foo);
  Field Bar(&Fine);
  Field Check(&Fine);
  Field Diff(&Fine);

  typedef CartesianStencil<vobj,vobj,SimpleStencilParams> Stencil;
  SimpleStencilParams p;

  for(int dir=0;dir<Nd;dir++){

    int L = Fine._fdimensions[dir];

    for(int disp=-(L-1);disp<=L-1;disp++){

      int npoint=1;
      std::vector<int> directions(npoint,dir);
      std::vector<int> displacements(npoint,disp);

      Stencil myStencil(&Fine,npoint,0,directions,displacements,p);

      SimpleCompressor<vobj> compress;
      myStencil.HaloExchange(Foo,compress);

      Bar = Cshift(Foo,dir,disp);

      {
	autoView( check , Check, AcceleratorWrite);
	autoView( foo   , Foo,   AcceleratorRead);
	autoView( st_v  , myStencil, AcceleratorRead);
	auto CBp=myStencil.CommBuf();
	accelerator_for(i,Check.Grid()->oSites(), 1, {

	    int permute_type;
	    StencilEntry *SE;
	    SE = st_v.GetEntry(permute_type,0,i);

	    if ( SE->_is_local && SE->_permute )
	      permute(check[i],foo[SE->_offset],permute_type);
	    else if (SE->_is_local)
	      check[i] = foo[SE->_offset];
	    else
	      check[i] = CBp[SE->_offset];
	});
      }

      Diff = Check-Bar;
      RealD nrm = norm2(Diff);

      if ( nrm > tol ) {
	std::cout << GridLogMessage << "FAIL dir " << dir << " disp " << disp
		  << " norm2(stencil-cshift) = " << nrm << std::endl;
      }
      GRID_ASSERT( nrm < tol );
    }
    std::cout << GridLogMessage << "Stencil == Cshift, dir " << dir
	      << ", all displacements " << -(L-1) << " .. " << L-1 << std::endl;
  }

  std::cout << GridLogMessage << "Test_stencil_lex: ALL PASS" << std::endl;

  Grid_finalize();
}
