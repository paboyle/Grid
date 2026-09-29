/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_coarse.cc

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
// MultiGeneralCoarsenedOperator against the existing mrhs coarse operator.
//
// The reference (the deprecated operator) is constructed on the D+1 grid; the
// operator under test on the D dimensional grid, with SetGrid() adopting the
// caller owned D+1 grid and building its padded cell and neighbour table from
// the D dimensional stencil, with the Nrhs factor multiplied in.
//
// Both are given identical matrix elements, so any difference in the apply is
// the restructured neighbour table. The same geometry object is passed to
// both: the reference adds one to skip for the rhs direction, the operator
// under test uses it as is on the D dimensional grid, so both describe the
// same stencil over the D dimensions.
//
#include <Grid/Grid.h>

using namespace Grid;

const int nbasis = 8;

typedef vSpinColourVector FineObj;
typedef sTComplexD        CComplexT;   // unvectorised coarse space

typedef DeprecatedMultiGeneralCoarsenedMatrix    <FineObj,CComplexT,nbasis> RefOperator;
typedef MultiGeneralCoarsenedOperator<FineObj,CComplexT,nbasis> TestOperator;

////////////////////////////////////////////////////////////////////////
// Identical random matrix elements into both operators
////////////////////////////////////////////////////////////////////////
template<class OpA,class OpB>
void SeedMatrixElements(OpA &A,OpB &B,int npoint,GridSerialRNG &sRNG)
{
  typedef typename OpA::calcMatrix calcMatrix;

  deviceVector<calcMatrix> bufA,bufB;
  for(int p=0;p<npoint;p++){

    // MatrixPointOut/In carry one stencil point in lSite order whichever
    // internal layout the operator uses.
    A.MatrixPointOut(p,bufA);
    B.MatrixPointOut(p,bufB);
    GRID_ASSERT(bufA.size() == bufB.size());

    int64_t sites = bufA.size();
    std::vector<calcMatrix> host(sites);

    ComplexD *w = (ComplexD *)&host[0];
    int64_t words = sites*sizeof(calcMatrix)/sizeof(ComplexD);
    for(int64_t i=0;i<words;i++){
      RealD re,im;
      random(sRNG,re);
      random(sRNG,im);
      w[i] = ComplexD(re-0.5,im-0.5);
    }

    acceleratorCopyToDevice(&host[0],&bufA[0],sites*sizeof(calcMatrix));
    A.MatrixPointIn(p,bufA);
    B.MatrixPointIn(p,bufA);
  }
}

template<class Op>
RealD MatrixChecksum(Op &O,int npoint)
{
  typedef typename Op::calcMatrix calcMatrix;
  RealD sum=0.0;
  deviceVector<calcMatrix> buf;
  for(int p=0;p<npoint;p++){
    O.MatrixPointOut(p,buf);
    int64_t sites = buf.size();
    std::vector<calcMatrix> host(sites);
    acceleratorCopyFromDevice(&buf[0],&host[0],sites*sizeof(calcMatrix));
    ComplexD *w = (ComplexD *)&host[0];
    int64_t words = sites*sizeof(calcMatrix)/sizeof(ComplexD);
    for(int64_t i=0;i<words;i++) sum += real(w[i])*real(w[i]) + imag(w[i])*imag(w[i]);
  }
  return sum;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int nrhs = 4;

  Coordinate clatt = GridDefaultLatt();
  Coordinate csimd(Nd,1);                 // the coarse space is unvectorised
  Coordinate cmpi  = GridDefaultMpi();

  ////////////////////////////////////////////////
  // D dimensional coarse grid, and D+1 for the reference
  ////////////////////////////////////////////////
  GridCartesian *CoarseD = new GridCartesian(clatt,csimd,cmpi);

  std::cout << GridLogMessage << "coarse D   grid "; for(int d=0;d<Nd;d++) std::cout<<clatt[d]<<" ";
  std::cout << "  Nsimd " << CoarseD->Nsimd() << std::endl;

  ////////////////////////////////////////////////
  // One geometry object for both, and one D+1 grid
  // owned here and shared by both operators: fields
  // conform only across a shared grid object.
  ////////////////////////////////////////////////
  NextToNearestStencilGeometry4D geom(CoarseD);

  Coordinate mlatt(1,nrhs), msimd(1,1), mmpi(1,1);
  for(int d=0;d<Nd;d++){
    mlatt.push_back(clatt[d]);
    msimd.push_back(csimd[d]);
    mmpi .push_back(cmpi[d]);
  }
  GridCartesian *CoarseMulti = new GridCartesian(mlatt,msimd,mmpi);

  TestOperator OpTest(geom,CoarseD);
  RefOperator OpRef(geom,CoarseMulti);
  OpTest.SetGrid(CoarseMulti);

  std::cout << GridLogMessage << "coarse D+1 grid nrhs " << nrhs
	    << "  Nsimd " << CoarseMulti->Nsimd() << std::endl;

  std::cout << GridLogMessage << "npoint ref " << OpRef.geom.npoint
	    << "   npoint test " << OpTest.geom.npoint << std::endl;
  GRID_ASSERT(OpRef.geom.npoint == OpTest.geom.npoint);

  int npoint = OpRef.geom.npoint;

  ////////////////////////////////////////////////
  // Identical matrix elements
  ////////////////////////////////////////////////
  GridSerialRNG sRNG; sRNG.SeedFixedIntegers(std::vector<int>({7,8,9,10}));
  SeedMatrixElements(OpRef,OpTest,npoint,sRNG);

  RealD ckRef = MatrixChecksum(OpRef,npoint);
  RealD ckTest = MatrixChecksum(OpTest,npoint);
  std::cout << GridLogMessage << "matrix element checksum ref " << ckRef
	    << "  test " << ckTest << std::endl;
  GRID_ASSERT( ckRef == ckTest );

  ////////////////////////////////////////////////
  // Same input, compare the applies
  ////////////////////////////////////////////////
  typedef RefOperator::CoarseVector CoarseVector;
  // RNG on the D dimensional grid fills any D+1 field: the rhs direction is
  // undistributed and divides cleanly. One RNG serves every Nrhs.
  GridParallelRNG pRNG(CoarseD); pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  CoarseVector in  (CoarseMulti); random(pRNG,in);
  CoarseVector out1(CoarseMulti);
  CoarseVector out2(CoarseMulti);
  CoarseVector err (CoarseMulti);

  OpRef.M(in,out1);
  OpTest.M(in,out2);

  err = out1 - out2;
  std::cout << GridLogMessage << "|ref out|^2 = " << norm2(out1)
	    << "   |test out|^2 = " << norm2(out2) << std::endl;
  std::cout << GridLogMessage << "|ref - test|^2 = " << norm2(err) << std::endl;
  GRID_ASSERT( norm2(out1) > 0.0 );
  GRID_ASSERT( norm2(err) == 0.0 );

  ////////////////////////////////////////////////
  // SetGrid is idempotent on pointer identity
  ////////////////////////////////////////////////
  OpTest.SetGrid(CoarseMulti);
  GRID_ASSERT( MatrixChecksum(OpTest,npoint) == ckTest );
  OpTest.M(in,out2);
  err = out1 - out2;
  GRID_ASSERT( norm2(err) == 0.0 );
  std::cout << GridLogMessage << "SetGrid idempotent on identity" << std::endl;

  ////////////////////////////////////////////////
  // Move to a different Nrhs and back. The matrix
  // elements are Nrhs independent and must survive
  // both the release and the rebuild.
  ////////////////////////////////////////////////
  // Nrhs 1 is the single RHS case through the multiRHS path, and each slice
  // of the Nrhs 4 apply must come back unchanged.
  OpTest.M(in,out2);
  for(int nr=2;nr>=1;nr--){
    Coordinate latt2(1,nr), simd2(1,1), mpi2(1,1);
    for(int d=0;d<Nd;d++){
      latt2.push_back(clatt[d]);
      simd2.push_back(csimd[d]);
      mpi2 .push_back(cmpi[d]);
    }
    GridCartesian *CoarseMulti2 = new GridCartesian(latt2,simd2,mpi2);

    OpTest.SetGrid(CoarseMulti2);
    GRID_ASSERT( OpTest.Nrhs() == nr );
    GRID_ASSERT( MatrixChecksum(OpTest,npoint) == ckTest );

    CoarseVector in2 (CoarseMulti2);
    CoarseVector out(CoarseMulti2);
    for(int r=0;r<nr;r++){
      CoarseVector slice(CoarseD);
      ExtractSliceFast(slice,in,r,0);
      InsertSliceFast(slice,in2,r,0);
    }
    OpTest.M(in2,out);

    RealD sdiff=0.0;
    for(int r=0;r<nr;r++){
      CoarseVector a(CoarseD),b(CoarseD),e(CoarseD);
      ExtractSliceFast(a,out ,r,0);
      ExtractSliceFast(b,out2,r,0);
      e = a-b;
      sdiff += norm2(e);
    }
    // Not bit exact: a different Nrhs is a different GEMM shape
    std::cout << GridLogMessage << "Nrhs " << nr << " slices agree with Nrhs "
	      << nrhs << " : |diff|^2/|out|^2 = " << sdiff/norm2(out) << std::endl;
    GRID_ASSERT( norm2(out) > 0.0 );
    GRID_ASSERT( sdiff/norm2(out) < 1.0e-20 );

    OpTest.ReleaseGrid();
    GRID_ASSERT( MatrixChecksum(OpTest,npoint) == ckTest );  // survives release

    delete CoarseMulti2;
  }

  OpTest.SetGrid(CoarseMulti);
  GRID_ASSERT( OpTest.Nrhs() == nrhs );

  OpTest.M(in,out2);
  err = out1 - out2;
  std::cout << GridLogMessage << "after Nrhs 4 -> 2 -> 1 -> release -> 4, |ref - test|^2 = "
	    << norm2(err) << std::endl;
  GRID_ASSERT( norm2(err) == 0.0 );

  std::cout << GridLogMessage << "Test_coarse: ALL PASS" << std::endl;

  Grid_finalize();
}
