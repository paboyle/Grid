/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_innerproduct_norm_lex.cc

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
// Lexicographic (Nsimd()==1) lattice types: linear combination, expression
// templates, norm2, innerProduct, innerProductNorm and RNG fill.
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-5;

template<class Field>
void BasicChecks(GridCartesian *grid,GridParallelRNG &pRNG,std::string precision)
{
  RealD vol = grid->gSites();

  Field a(grid);  a = 1.0;
  Field b(grid);  b = 2.0;
  Field c(grid);  c = 3.0;
  Field d(grid);

  ///////////////////////////////////////////////////
  // Linear combination through expression templates
  ///////////////////////////////////////////////////
  d = a + b - c;
  RealD n_zero = norm2(d);

  RealD n_a = norm2(a);
  RealD n_b = norm2(b);

  std::cout << GridLogMessage << precision << ": norm2(1+2-3) = " << n_zero << std::endl;
  std::cout << GridLogMessage << precision << ": norm2(1) = " << n_a << " expect " << vol << std::endl;
  std::cout << GridLogMessage << precision << ": norm2(2) = " << n_b << " expect " << 4.0*vol << std::endl;

  GRID_ASSERT( n_zero          < tol );
  GRID_ASSERT( fabs(n_a -     vol) < tol*vol );
  GRID_ASSERT( fabs(n_b - 4.0*vol) < tol*vol );

  ///////////////////////////////////////////////////
  // innerProduct of constant fields
  ///////////////////////////////////////////////////
  ComplexD ip = innerProduct(a,b);
  std::cout << GridLogMessage << precision << ": innerProduct(1,2) = " << ip
            << " expect (" << 2.0*vol << ",0)" << std::endl;
  GRID_ASSERT( fabs(real(ip) - 2.0*vol) < tol*vol );
  GRID_ASSERT( fabs(imag(ip))           < tol*vol );

  ///////////////////////////////////////////////////
  // RNG fill; norm2(r+r) == 4 norm2(r)
  ///////////////////////////////////////////////////
  Field r(grid); random(pRNG,r);
  Field s(grid); random(pRNG,s);

  RealD n_r  = norm2(r);
  Field rr(grid); rr = r + r;
  RealD n_rr = norm2(rr);

  std::cout << GridLogMessage << precision << ": norm2(r) = " << n_r
            << " norm2(r+r) = " << n_rr << " ratio " << n_rr/n_r << std::endl;
  GRID_ASSERT( n_r > tol );                          // RNG actually filled it
  GRID_ASSERT( fabs(n_rr - 4.0*n_r) < tol*n_rr );

  ///////////////////////////////////////////////////
  // Fused innerProductNorm against separate calls
  ///////////////////////////////////////////////////
  ComplexD ip_ref = innerProduct(r,s);
  RealD    n2_ref = norm2(r);

  ComplexD ip_res;
  RealD    n2_res;
  innerProductNorm(ip_res,n2_res,r,s);

  std::cout << GridLogMessage << precision << ": innerProductNorm ip diff "
            << abs(ip_ref-ip_res) << " norm2 diff " << fabs(n2_ref-n2_res) << std::endl;
  GRID_ASSERT( abs(ip_ref-ip_res)   < tol*abs(ip_ref) + tol );
  GRID_ASSERT( fabs(n2_ref-n2_res)  < tol*n2_ref );

  std::cout << GridLogMessage << precision << ": all checks passed" << std::endl;
}

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  Coordinate latt   = GridDefaultLatt();
  Coordinate simd({1,1,1,1});          // lexicographic: no SIMD
  Coordinate mpi    = GridDefaultMpi();

  GridCartesian grid(latt,simd,mpi);

  std::cout << GridLogMessage << "Lexicographic grid, Nsimd = " << grid.Nsimd() << std::endl;
  GRID_ASSERT( grid.Nsimd() == 1 );

  GridParallelRNG pRNG(&grid);
  pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  

  BasicChecks<lexLatticeComplexD>(&grid,pRNG,"Double");
  BasicChecks<lexLatticeComplexF>(&grid,pRNG,"Single");

  std::cout << GridLogMessage << "Test_innerproduct_norm_lex: ALL PASS" << std::endl;

  Grid_finalize();
}
