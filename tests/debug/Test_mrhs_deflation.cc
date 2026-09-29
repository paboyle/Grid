/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_mrhs_deflation.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

// MultiRHSDeflation: the D+1 multiRHS interface (one permutation pass each
// way) must reproduce the vector-of-D-fields interface bit for bit up to
// the summation order of the GEMMs.  Random "eigenvectors" and values: the
// deflation formula G = E (E^dag R)/lambda does not care that they are not
// eigenpairs of anything.

using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int nbasis = 8;
  const int nev    = 12;
  const int nrhs   = 5;

  typedef iVector<sTComplexD,nbasis> siteVector;
  typedef Lattice<siteVector>        CoarseVector;

  // Unvectorised D and D+1 coarse grids, as MGCoarseGrids builds them
  Coordinate latt({1,4,4,4,4});
  Coordinate simd({1,1,1,1,1});
  Coordinate mpi = GridDefaultMpi();
  Coordinate mpi5({1,mpi[0],mpi[1],mpi[2],mpi[3]});
  GridCartesian Coarse5d(latt,simd,mpi5);

  Coordinate latt6({nrhs,1,4,4,4,4});
  Coordinate simd6({1,1,1,1,1,1});
  Coordinate mpi6({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  GridCartesian CoarseMrhs(latt6,simd6,mpi6);

  GridParallelRNG RNG(&Coarse5d); RNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  std::vector<CoarseVector> evec(nev,&Coarse5d);
  std::vector<RealD>        eval(nev);
  for(int e=0;e<nev;e++){ random(RNG,evec[e]); eval[e] = 0.5 + 0.1*e; }

  std::vector<CoarseVector> src(nrhs,&Coarse5d), guess(nrhs,&Coarse5d);
  for(int r=0;r<nrhs;r++) random(RNG,src[r]);

  MultiRHSDeflation<CoarseVector> Deflator;
  Deflator.ImportEigenBasis(evec,eval);

  // Reference: the vector interface
  Deflator.DeflateSources(src,guess);

  // The D+1 interface on the same sources
  CoarseVector src_mrhs(&CoarseMrhs), guess_mrhs(&CoarseMrhs);
  for(int r=0;r<nrhs;r++) InsertSliceFast(src[r],src_mrhs,r,0);
  guess_mrhs = Zero();
  Deflator.DeflateSources(src_mrhs,guess_mrhs);

  RealD worst = 0.0;
  CoarseVector g(&Coarse5d), d(&Coarse5d);
  for(int r=0;r<nrhs;r++){
    ExtractSliceFast(g,guess_mrhs,r,0);
    d = g - guess[r];
    RealD rel = std::sqrt(norm2(d)/norm2(guess[r]));
    std::cout << GridLogMessage << "rhs " << r << " ||guess_mrhs - guess||/||guess|| = " << rel << std::endl;
    worst = std::max(worst,rel);
  }
  std::cout << GridLogMessage << "MultiRHSDeflation D+1 vs vector interface: worst " << worst << std::endl;
  GRID_ASSERT( worst < 1.0e-12 );
  std::cout << GridLogMessage << "PASS" << std::endl;

  Grid_finalize();
  return 0;
}
