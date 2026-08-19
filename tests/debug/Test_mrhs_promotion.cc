/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_mrhs_promotion.cc

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
// MrhsPromotedOperator: a D dimensional operator presented as D+1 with Nrhs.
//
// Each output slice must equal the D dimensional operator applied to the
// corresponding input slice, bit for bit -- the promotion moves data, it does
// not change arithmetic.  Slices are independent random fields so a mistake in
// the slice indexing cannot pass.
//
#include <Grid/Grid.h>

using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Nrhs = 4;

  Coordinate latt  = GridDefaultLatt();
  Coordinate simd  = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi   = GridDefaultMpi();

  GridCartesian         *UGrid  = new GridCartesian(latt,simd,mpi);
  GridRedBlackCartesian *UrbGrid= SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  // D+1 grid with Nrhs in dimension 0
  GridCartesian         *RGrid  = SpaceTimeGrid::makeFiveDimGrid(Nrhs,UGrid);

  std::cout << GridLogMessage << "Nrhs = " << Nrhs
	    << "  low grid " << UGrid->_ndimension << "d"
	    << "  high grid " << RGrid->_ndimension << "d" << std::endl;

  GridParallelRNG pRNG(UGrid); pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  LatticeGaugeFieldD Umu(UGrid); SU<Nc>::HotConfiguration(pRNG,Umu);

  RealD mass = 0.1;
  WilsonFermionD Dw(Umu,*UGrid,*UrbGrid,mass);

  MdagMLinearOperator<WilsonFermionD,LatticeFermionD> HermOp(Dw);

  MrhsPromotedOperator<LatticeFermionD> MrhsOp(HermOp,UGrid,Nrhs);

  //////////////////////////////////////////////////
  // Independent random source on each slice
  //////////////////////////////////////////////////
  LatticeFermionD hi_in (RGrid);
  LatticeFermionD hi_out(RGrid);
  std::vector<LatticeFermionD> lo_in (Nrhs,UGrid);
  std::vector<LatticeFermionD> lo_ref(Nrhs,UGrid);

  for(int r=0;r<Nrhs;r++){
    random(pRNG,lo_in[r]);
    InsertSliceFast(lo_in[r],hi_in,r,0);
  }

  //////////////////////////////////////////////////
  // Promoted application, and the reference
  //////////////////////////////////////////////////
  MrhsOp.Op(hi_in,hi_out);

  for(int r=0;r<Nrhs;r++){
    HermOp.Op(lo_in[r],lo_ref[r]);
  }

  RealD worst = 0.0;
  for(int r=0;r<Nrhs;r++){
    LatticeFermionD slice(UGrid);
    LatticeFermionD err(UGrid);
    ExtractSliceFast(slice,hi_out,r,0);
    err = slice - lo_ref[r];
    RealD n = norm2(err);
    std::cout << GridLogMessage << "rhs["<<r<<"] |promoted - reference|^2 = " << n
	      << "   |reference|^2 = " << norm2(lo_ref[r]) << std::endl;
    worst = std::max(worst,n);
    GRID_ASSERT( n == 0.0 );
  }
  std::cout << GridLogMessage << "worst slice difference " << worst << " (bit exact required)" << std::endl;

  //////////////////////////////////////////////////
  // Slices must be independent: perturbing one must
  // change only that output slice
  //////////////////////////////////////////////////
  {
    LatticeFermionD pert(UGrid);
    random(pRNG,pert);
    InsertSliceFast(pert,hi_in,1,0);
    MrhsOp.Op(hi_in,hi_out);

    for(int r=0;r<Nrhs;r++){
      LatticeFermionD slice(UGrid);
      LatticeFermionD err(UGrid);
      ExtractSliceFast(slice,hi_out,r,0);
      err = slice - lo_ref[r];
      RealD n = norm2(err);
      if ( r == 1 ) {
	std::cout << GridLogMessage << "perturbed rhs[1] changed by " << n << std::endl;
	GRID_ASSERT( n > 0.0 );
      } else {
	std::cout << GridLogMessage << "untouched rhs["<<r<<"] moved by " << n << std::endl;
	GRID_ASSERT( n == 0.0 );
      }
    }
  }

  //////////////////////////////////////////////////
  // AdjOp routes through the same slicing
  //////////////////////////////////////////////////
  {
    for(int r=0;r<Nrhs;r++){
      random(pRNG,lo_in[r]);
      InsertSliceFast(lo_in[r],hi_in,r,0);
      HermOp.AdjOp(lo_in[r],lo_ref[r]);
    }
    MrhsOp.AdjOp(hi_in,hi_out);
    for(int r=0;r<Nrhs;r++){
      LatticeFermionD slice(UGrid);
      LatticeFermionD err(UGrid);
      ExtractSliceFast(slice,hi_out,r,0);
      err = slice - lo_ref[r];
      GRID_ASSERT( norm2(err) == 0.0 );
    }
    std::cout << GridLogMessage << "AdjOp slices bit exact" << std::endl;
  }

  std::cout << GridLogMessage << "Test_mrhs_promotion: ALL PASS" << std::endl;

  Grid_finalize();
}
