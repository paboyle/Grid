/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_coarse_coarsen.cc

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
// Coarsen the same fine operator two ways and compare the matrix elements.
//
//   ref  : DeprecatedMultiGeneralCoarsenedMatrix::CoarsenOperator, vectorised
//          coarse space matched to the fine SIMD layout, single RHS fine
//          applications
//   test : MultiGeneralCoarsenedOperator::CoarsenOperator on the D+1 grid,
//          unvectorised sComplexD coarse space, the batch of phased basis
//          vectors carried in the rhs direction and applied through
//          MrhsPromotedOperator
//
// BLAS_A is written by GridtoBLAS in lSite order, which does not depend on
// the SIMD layout, so the two are directly comparable.
//
#include <Grid/Grid.h>

using namespace Grid;

const int nbasis = 8;
const int batch  = 9;

typedef vSpinColourVector FineObj;
typedef vTComplex         CComplexV;   // vectorised coarse space, reference
typedef sTComplexD        CComplexS;   // unvectorised coarse space, under test

typedef DeprecatedMultiGeneralCoarsenedMatrix    <FineObj,CComplexV,nbasis> RefOperator;
typedef MultiGeneralCoarsenedOperator<FineObj,CComplexS,nbasis> TestOperator;

template<class Op>
void ReadMatrix(Op &O,int npoint,std::vector<std::vector<typename Op::calcMatrix> > &host)
{
  typedef typename Op::calcMatrix calcMatrix;
  host.resize(npoint);
  // MatrixPointOut carries one stencil point in lSite order whichever internal
  // layout the operator uses.
  deviceVector<calcMatrix> buf;
  for(int p=0;p<npoint;p++){
    O.MatrixPointOut(p,buf);
    int64_t sites = buf.size();
    host[p].resize(sites);
    acceleratorCopyFromDevice(&buf[0],&host[p][0],sites*sizeof(calcMatrix));
  }
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate flatt = GridDefaultLatt();
  Coordinate fsimd = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate fmpi  = GridDefaultMpi();

  Coordinate block({2,2,2,2});
  Coordinate clatt(Nd);
  for(int d=0;d<Nd;d++){
    GRID_ASSERT(flatt[d]%block[d]==0);
    clatt[d] = flatt[d]/block[d];
  }

  GridCartesian         *FineGrid   = new GridCartesian(flatt,fsimd,fmpi);
  GridRedBlackCartesian *FrbGrid    = SpaceTimeGrid::makeFourDimRedBlackGrid(FineGrid);

  ////////////////////////////////////////////////
  // Reference coarse space: SIMD layout matched to the fine
  ////////////////////////////////////////////////
  Coordinate cvsimd = GridDefaultSimd(Nd,CComplexV::Nsimd());
  GridCartesian *CoarseV = new GridCartesian(clatt,cvsimd,fmpi);

  // The reference puts all of the SIMD in the rhs direction. CoarsenOperator
  // does not use the multiRHS grid; it only sizes BLAS_A, so one lane of rhs suffices.
  int nrhs_ref = CComplexV::Nsimd();
  Coordinate reflatt(1,nrhs_ref),refsimd(1,CComplexV::Nsimd()),refmpi(1,1);
  for(int d=0;d<Nd;d++){
    reflatt.push_back(clatt[d]);
    refsimd.push_back(1);
    refmpi .push_back(fmpi[d]);
  }
  GridCartesian *CoarseVMulti = new GridCartesian(reflatt,refsimd,refmpi);

  ////////////////////////////////////////////////
  // Coarse space under test: unvectorised
  ////////////////////////////////////////////////
  Coordinate cssimd(Nd,1);
  GridCartesian *CoarseS = new GridCartesian(clatt,cssimd,fmpi);

  Coordinate cmlatt(1,batch),cmsimd(1,1),cmmpi(1,1);
  for(int d=0;d<Nd;d++){
    cmlatt.push_back(clatt[d]);
    cmsimd.push_back(1);
    cmmpi .push_back(fmpi[d]);
  }
  GridCartesian *CoarseSMulti = new GridCartesian(cmlatt,cmsimd,cmmpi);

  // D+1 fine grid carrying the batch
  Coordinate fmlatt(1,batch), fmsimd(1,1), fmmpi(1,1);
  for(int d=0;d<Nd;d++){
    fmlatt.push_back(flatt[d]);
    fmsimd.push_back(fsimd[d]);
    fmmpi .push_back(fmpi[d]);
  }
  GridCartesian *FineGridMulti = new GridCartesian(fmlatt,fmsimd,fmmpi);

  std::cout << GridLogMessage << "fine "<<flatt<<"   coarse "<<clatt<<"   batch "<<batch<<std::endl;
  std::cout << GridLogMessage << "Nsimd   fine "<<FineGrid->Nsimd()
	    << "   coarse ref "<<CoarseV->Nsimd()
	    << "   coarse test "<<CoarseS->Nsimd()<<std::endl;

  ////////////////////////////////////////////////
  // Fine operator
  ////////////////////////////////////////////////
  GridParallelRNG pRNG(FineGrid); pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));
  LatticeGaugeFieldD Umu(FineGrid); SU<Nc>::HotConfiguration(pRNG,Umu);

  RealD mass = 0.1;
  WilsonFermionD Dw(Umu,*FineGrid,*FrbGrid,mass);
  MdagMLinearOperator<WilsonFermionD,LatticeFermionD> HermOp(Dw);

  ////////////////////////////////////////////////
  // One random subspace, two copies: CoarsenOperator orthogonalises in place
  ////////////////////////////////////////////////
  Aggregation<FineObj,CComplexV,nbasis> Subspace(CoarseV,FineGrid,0);
  std::vector<LatticeFermionD> subspace(nbasis,FineGrid);
  for(int i=0;i<nbasis;i++){
    random(pRNG,Subspace.subspace[i]);
    subspace[i] = Subspace.subspace[i];
  }

  NextToNearestStencilGeometry4D geomV(CoarseV);
  NextToNearestStencilGeometry4D geomS(CoarseS);

  ////////////////////////////////////////////////
  // Reference coarsening, matched layouts
  ////////////////////////////////////////////////
  RefOperator OpRef(geomV,CoarseVMulti);
  std::cout << GridLogMessage << "reference CoarsenOperator" << std::endl;
  OpRef.CoarsenOperator(HermOp,Subspace,CoarseV);

  ////////////////////////////////////////////////
  // Coarsening under test, D+1 fine applications, unvectorised coarse
  ////////////////////////////////////////////////
  TestOperator OpTest(geomS,CoarseS);
  OpTest.SetGrid(CoarseSMulti);

  MrhsPromotedOperator<LatticeFermionD> MrhsHermOp(HermOp,FineGrid,batch);

  std::cout << GridLogMessage << "test CoarsenOperator (D+1)" << std::endl;
  OpTest.CoarsenOperator(MrhsHermOp,FineGridMulti,subspace,CoarseS);

  ////////////////////////////////////////////////
  // Compare matrix elements
  ////////////////////////////////////////////////
  int npoint = OpRef.geom.npoint;
  GRID_ASSERT(npoint == OpTest.geom.npoint);
  typedef RefOperator::calcMatrix calcMatrix;
  std::vector<std::vector<calcMatrix> > A1,A2;
  ReadMatrix(OpRef,npoint,A1);
  ReadMatrix(OpTest,npoint,A2);

  RealD num=0.0, den=0.0;
  for(int p=0;p<npoint;p++){
    GRID_ASSERT(A1[p].size()==A2[p].size());
    ComplexD *w1 = (ComplexD *)&A1[p][0];
    ComplexD *w2 = (ComplexD *)&A2[p][0];
    int64_t words = A1[p].size()*sizeof(calcMatrix)/sizeof(ComplexD);
    for(int64_t i=0;i<words;i++){
      ComplexD d = w1[i]-w2[i];
      num += real(d)*real(d)+imag(d)*imag(d);
      den += real(w1[i])*real(w1[i])+imag(w1[i])*imag(w1[i]);
    }
  }
  std::cout << GridLogMessage << "|A_ref|^2 = " << den << std::endl;
  std::cout << GridLogMessage << "|A_ref - A_test|^2 / |A_ref|^2 = " << num/den << std::endl;
  GRID_ASSERT( den > 0.0 );
  GRID_ASSERT( num/den < 1.0e-20 );

  ////////////////////////////////////////////////
  // Same coarsening through the single RHS variant: no mrhs packing, the
  // batch assembled on the coarse side by the mixed blockProject. Block
  // Gram-Schmidt is idempotent so the subspace may be reused in place.
  ////////////////////////////////////////////////
  TestOperator OpTestSrhs(geomS,CoarseS);
  OpTestSrhs.SetGrid(CoarseSMulti);

  std::cout << GridLogMessage << "test CoarsenOperator (single RHS fine op)" << std::endl;
  OpTestSrhs.CoarsenOperator(HermOp,subspace,CoarseS,batch);

  std::vector<std::vector<calcMatrix> > A3;
  ReadMatrix(OpTestSrhs,npoint,A3);

  RealD nums=0.0;
  for(int p=0;p<npoint;p++){
    GRID_ASSERT(A1[p].size()==A3[p].size());
    ComplexD *w1 = (ComplexD *)&A1[p][0];
    ComplexD *w3 = (ComplexD *)&A3[p][0];
    int64_t words = A1[p].size()*sizeof(calcMatrix)/sizeof(ComplexD);
    for(int64_t i=0;i<words;i++){
      ComplexD d = w1[i]-w3[i];
      nums += real(d)*real(d)+imag(d)*imag(d);
    }
  }
  std::cout << GridLogMessage << "|A_ref - A_test_srhs|^2 / |A_ref|^2 = " << nums/den << std::endl;
  GRID_ASSERT( nums/den < 1.0e-20 );

  ////////////////////////////////////////////////
  // GetMatrix/SetMatrix round trip through the BLAS layout array. This is
  // how a bilingual DenseCoarseMatrix will retrieve the elements, and it
  // needs no SetGrid since the matrix elements are Nrhs independent.
  ////////////////////////////////////////////////
  {
    typedef TestOperator::CoarseMatrix CoarseMatrixS;
    std::vector<CoarseMatrixS> Aget(npoint,CoarseS);
    for(int p=0;p<npoint;p++) OpTest.GetMatrix(p,Aget);

    TestOperator OpTestRoundTrip(geomS,CoarseS);
    for(int p=0;p<npoint;p++) OpTestRoundTrip.SetMatrix(p,Aget);

    std::vector<std::vector<calcMatrix> > A4;
    ReadMatrix(OpTestRoundTrip,npoint,A4);

    RealD numrt=0.0;
    for(int p=0;p<npoint;p++){
      GRID_ASSERT(A2[p].size()==A4[p].size());
      ComplexD *w2 = (ComplexD *)&A2[p][0];
      ComplexD *w4 = (ComplexD *)&A4[p][0];
      int64_t words = A2[p].size()*sizeof(calcMatrix)/sizeof(ComplexD);
      for(int64_t i=0;i<words;i++){
        ComplexD d = w2[i]-w4[i];
        numrt += real(d)*real(d)+imag(d)*imag(d);
      }
    }
    std::cout << GridLogMessage << "GetMatrix/SetMatrix round trip |diff|^2 = " << numrt << std::endl;
    GRID_ASSERT( numrt == 0.0 );
  }

  ////////////////////////////////////////////////
  // and the operator under test applies the matrix it just built
  ////////////////////////////////////////////////
  typedef TestOperator::CoarseVector CoarseVectorS;
  GridParallelRNG cRNG(CoarseS); cRNG.SeedFixedIntegers(std::vector<int>({5,6,7,8}));
  CoarseVectorS in(CoarseSMulti);  random(cRNG,in);
  CoarseVectorS out(CoarseSMulti);

  OpTest.M(in,out);
  std::cout << GridLogMessage << "|in|^2 = " << norm2(in)
	    << "   |M_test in|^2 = " << norm2(out) << std::endl;
  GRID_ASSERT( norm2(out) > 0.0 );

  std::cout << GridLogMessage << "Test_coarse_coarsen: ALL PASS" << std::endl;

  Grid_finalize();
}
