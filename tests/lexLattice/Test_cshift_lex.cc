/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/lexLattice/Test_cshift_lex.cc

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
// Lexicographic (Nsimd()==1) lattice types: LatticeCoordinate, Cshift in all
// directions and shifts, and predicated where().
//
// Cshift is checked by an identity rather than a reference loop: with
//    d = Cshift(coor_mu,mu,shift) - coor_mu
// every site has d = shift (no wrap) or d = shift - L (wrapped), so
//    (d - shift)*(d - shift + L) == 0
// everywhere.
//
#include <Grid/Grid.h>

using namespace Grid;

const RealD tol = 1.0e-5;

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  Coordinate latt   = GridDefaultLatt();
  Coordinate simd({1,1,1,1});          // lexicographic: no SIMD
  Coordinate mpi    = GridDefaultMpi();

  GridCartesian grid(latt,simd,mpi);

  std::cout << GridLogMessage << "Lexicographic grid, Nsimd = " << grid.Nsimd() << std::endl;
  GRID_ASSERT( grid.Nsimd() == 1 );

  

  //////////////////////////////////////////////////
  // Cshift against the wrap identity
  //////////////////////////////////////////////////
  for(int mu=0;mu<Nd;mu++){

    lexLatticeComplexD coor(&grid);
    LatticeCoordinate(coor,mu);

    RealD L = latt[mu];

    for(int shift=0;shift<latt[mu];shift++){

      lexLatticeComplexD shifted(&grid);
      shifted = Cshift(coor,mu,shift);

      lexLatticeComplexD d(&grid);
      d = shifted - coor;

      lexLatticeComplexD p(&grid);
      p = (d - ComplexD(shift)) * (d - ComplexD(shift) + ComplexD(L));

      RealD n = norm2(p);
      if ( n > tol ) {
        std::cout << GridLogMessage << "FAIL mu " << mu << " shift " << shift
                  << " residual " << n << std::endl;
      }
      GRID_ASSERT( n < tol );
    }
    std::cout << GridLogMessage << "Cshift mu = " << mu << " all shifts pass" << std::endl;
  }

  //////////////////////////////////////////////////
  // Predicated where(): zero the second half in time
  //////////////////////////////////////////////////
  {
    int    Tdir = Nd-1;
    RealD  T    = latt[Tdir];
    RealD  vol  = grid.gSites();

    lexLatticeInteger  tcoor(&grid);
    LatticeCoordinate(tcoor,Tdir);

    lexLatticeComplexD f(&grid);   f  = 2.0;
    lexLatticeComplexD zz(&grid);  zz = Zero();
    lexLatticeComplexD g(&grid);

    g = where( tcoor < Integer(T/2) , f , zz );

    RealD n_g = norm2(g);
    RealD expect = 4.0 * vol / 2.0;

    std::cout << GridLogMessage << "where(): norm2 = " << n_g
              << " expect " << expect << std::endl;
    DumpSliceNorm("where() slice norm",g,Tdir);

    GRID_ASSERT( fabs(n_g - expect) < tol*expect );
  }

  //////////////////////////////////////////////////
  // Cshift against sliceSum of a slice-zeroed random field.
  // Exhaustive over direction, zeroed slice and shift:
  //    sliceSum(Cshift(f,mu,s))[t] == sliceSum(f)[(t+s)%L]
  //////////////////////////////////////////////////
  {
    GridParallelRNG pRNG(&grid);
    pRNG.SeedFixedIntegers(std::vector<int>({7,8,9,10}));

    typedef lexLatticeComplexD::vector_object::scalar_object sobj;

    lexLatticeComplexD rnd(&grid);  random(pRNG,rnd);
    lexLatticeComplexD zz(&grid);   zz = Zero();

    for(int mu=0;mu<Nd;mu++){

      int L = latt[mu];

      lexLatticeInteger coor(&grid);
      LatticeCoordinate(coor,mu);

      for(int t0=0;t0<L;t0++){

        lexLatticeComplexD f(&grid);
        f = where( coor == Integer(t0) , zz , rnd );

        std::vector<sobj> ref;
        sliceSum(f,ref,mu);

        for(int shift=0;shift<L;shift++){

          lexLatticeComplexD g(&grid);
          g = Cshift(f,mu,shift);

          std::vector<sobj> res;
          sliceSum(g,res,mu);

          for(int t=0;t<L;t++){
            ComplexD got    = TensorRemove(res[t]);
            ComplexD expect = TensorRemove(ref[(t+shift)%L]);
            if ( abs(got-expect) > tol*(abs(expect)+1.0) ) {
              std::cout << GridLogMessage << "FAIL mu " << mu << " t0 " << t0
                        << " shift " << shift << " t " << t
                        << " got " << got << " expect " << expect << std::endl;
            }
            GRID_ASSERT( abs(got-expect) < tol*(abs(expect)+1.0) );
          }
        }
      }
      std::cout << GridLogMessage << "sliceSum/Cshift consistency mu = " << mu
                << " all zeroed slices and shifts pass" << std::endl;
    }
  }

  std::cout << GridLogMessage << "Test_cshift_lex: ALL PASS" << std::endl;

  Grid_finalize();
}
