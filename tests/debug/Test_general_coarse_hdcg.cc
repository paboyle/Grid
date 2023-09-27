    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_general_coarse_hdcg.cc

    Copyright (C) 2023

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
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/GeneralCoarsenedMatrix.h>
#include <Grid/algorithms/iterative/AdefGeneric.h>

using namespace std;
using namespace Grid;

template<class Field> class TestSolver : public LinearFunction<Field> {
public:
  TestSolver() {};
  void operator() (const Field &in, Field &out){    out = Zero();  }     
};


RealD InverseApproximation(RealD x){
  return 1.0/x;
}

// Want Op in CoarsenOp to call MatPcDagMatPc
template<class Field>
class HermOpAdaptor : public LinearOperatorBase<Field>
{
  LinearOperatorBase<Field> & wrapped;
public:
  HermOpAdaptor(LinearOperatorBase<Field> &wrapme) : wrapped(wrapme)  {};
  void Op     (const Field &in, Field &out)   { wrapped.HermOp(in,out);  }
  void HermOp(const Field &in, Field &out)    { wrapped.HermOp(in,out); }
  void AdjOp     (const Field &in, Field &out){ wrapped.HermOp(in,out);  }
  void OpDiag (const Field &in, Field &out)                  {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out)  {    assert(0);  };
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
};
template<class Field,class Matrix> class ChebyshevSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _SmootherOperator;
  Chebyshev<Field> Cheby;
  ChebyshevSmoother(RealD _lo,RealD _hi,int _ord, FineOperator &SmootherOperator) :
    _SmootherOperator(SmootherOperator),
    Cheby(_lo,_hi,_ord,InverseApproximation)
  {
    std::cout << GridLogMessage<<" Chebyshev smoother order "<<_ord<<" ["<<_lo<<","<<_hi<<"]"<<std::endl;
  };
  void operator() (const Field &in, Field &out) 
  {
    Field tmp(in.Grid());
    tmp = in;
    Cheby(_SmootherOperator,tmp,out);
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=16;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Construct a coarsened grid with 4^4 cell
  Coordinate clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/4;
  }
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt,
							    GridDefaultSimd(Nd,vComplex::Nsimd()),
							    GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  ///////////////////////// RNGs /////////////////////////////////
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});

  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  ///////////////////////// Configuration /////////////////////////////////
  LatticeGaugeField Umu(UGrid);

  FieldMetaData header;
  std::string file("ckpoint_lat.4000");
  NerscIO::readConfiguration(Umu,header,file);

  //////////////////////// Fermion action //////////////////////////////////
  RealD mass=0.01;
  RealD M5=1.8;
  RealD b=1.5;
  RealD c=0.5;
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);

  SchurDiagMooeeOperator<MobiusFermionD, LatticeFermion> HermOpEO(Ddwf);

  typedef HermOpAdaptor<LatticeFermionD> HermFineMatrix;
  HermFineMatrix FineHermOp(HermOpEO);
  
  LatticeFermion result(FrbGrid); result=Zero();

  LatticeFermion    src(FrbGrid); random(RNG5,src);

  // Run power method on FineHermOp
  PowerMethod<LatticeFermion>       PM;   PM(HermOpEO,src);

 
  ////////////////////////////////////////////////////////////
  ///////////// Coarse basis and Little Dirac Operator ///////
  ////////////////////////////////////////////////////////////
  const int nbasis = 40;
  const int cb = 0 ;
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNextToNextToNearestStencilGeometry5D geom(Coarse5d);
  
  // Warning: This routine calls PVdagM.Op, not PVdagM.HermOp
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FrbGrid,cb);
  Aggregates.CreateSubspaceChebyshev(RNG5,
				     HermOpEO,
				     nbasis,
				     //				     100.0,
				     //	0.1, // Low pass is pretty high still -- 311 iters
				     //				     250.0,
				     //				     0.01, // subspace too low filter power wrong
				     //				     250.0,
				     //				     0.2, // slower
				     95.0,
				     //				     0.05, // nbasis 12 - 311 -- wrong coarse inv
				     //				     0.05, // nbasis 12 - 154 -- right filt
				     //				     0.1, // nbasis 12 - 169 oops
				     //				     0.05, // nbasis 16 -- 127 iters
				     //				     0.03, // nbasis 16 -- 13-
				     //				     0.1,  // nbasis 16 -- 142; sloppy solve
				     0.1,  // nbasis 24 
				     300);
  ////////////////////////////////////////////////////////////
  // Need to check about red-black grid coarsening
  ////////////////////////////////////////////////////////////
  LittleDiracOperator LittleDiracOp(geom,FrbGrid,Coarse5d);
  LittleDiracOp.CoarsenOperatorColoured(FineHermOp,Aggregates);

  // Try projecting to one hop only
  LittleDiracOperator LittleDiracOpProj(LittleDiracOp);
  LittleDiracOpProj.ProjectNearestNeighbour(0.5);

  typedef HermitianLinearOperator<LittleDiracOperator,CoarseVector> HermMatrix;
  HermMatrix CoarseOp (LittleDiracOp);

  //////////////////////////////////////////
  // Build a coarse lanczos
  //////////////////////////////////////////
  Chebyshev<CoarseVector>      IRLCheby(0.5,60.0,71);  // 1 iter
  FunctionHermOp<CoarseVector> IRLOpCheby(IRLCheby,CoarseOp);
  PlainHermOp<CoarseVector>    IRLOp    (CoarseOp);
  int Nk=48;
  int Nm=64;
  int Nstop=Nk;
  ImplicitlyRestartedLanczos<CoarseVector> IRL(IRLOpCheby,IRLOp,Nstop,Nk,Nm,1.0e-5,20);

  int Nconv;
  std::vector<RealD>            eval(Nm);
  std::vector<CoarseVector>     evec(Nm,Coarse5d);
  CoarseVector c_src(Coarse5d); c_src=1.0;

  PowerMethod<CoarseVector>       cPM;   cPM(CoarseOp,c_src);

  IRL.calc(eval,evec,c_src,Nconv);
  DeflatedGuesser<CoarseVector> DeflCoarseGuesser(evec,eval);
  
  //////////////////////////////////////////
  // Build a coarse space solver
  //////////////////////////////////////////
  int maxit=20000;
  ConjugateGradient<CoarseVector>  CG(1.0e-8,maxit,false);
  ConjugateGradient<LatticeFermionD>  CGfine(1.0e-8,10000,false);
  ZeroGuesser<CoarseVector> CoarseZeroGuesser;

  //  HPDSolver<CoarseVector> HPDSolve(CoarseOp,CG,CoarseZeroGuesser);
  HPDSolver<CoarseVector> HPDSolve(CoarseOp,CG,DeflCoarseGuesser);

  //////////////////////////////////////////
  // Build a smoother
  //////////////////////////////////////////
  //  ChebyshevSmoother<LatticeFermionD,HermFineMatrix > Smoother(10.0,100.0,10,FineHermOp); //499
  //  ChebyshevSmoother<LatticeFermionD,HermFineMatrix > Smoother(3.0,100.0,10,FineHermOp);  //383
  //  ChebyshevSmoother<LatticeFermionD,HermFineMatrix > Smoother(1.0,100.0,10,FineHermOp);  //328
  //  std::vector<RealD> los({0.5,1.0,3.0}); // 147/142/146 nbasis 1
  //  std::vector<RealD> los({1.0,2.0}); // Nbasis 24: 88,86 iterations
  //  std::vector<RealD> los({2.0,4.0}); // Nbasis 32 == 52, iters
  //  std::vector<RealD> los({2.0,4.0}); // Nbasis 40 == 36,36 iters

  //
  // Turns approx 2700 iterations into 340 fine multiplies with Nbasis 40
  // Need to measure cost of coarse space.
  //
  // -- i) Reduce coarse residual   -- 0.04
  // -- ii) Lanczos on coarse space -- done
  // -- iii) Possible 1 hop project and/or preconditioning it - easy - PrecCG it and
  //         use a limited stencil. Reread BFM code to check on evecs / deflation strategy with prec
  //
  std::vector<RealD> los({3.0}); // Nbasis 40 == 36,36 iters

  //  std::vector<int> ords({7,8,10}); // Nbasis 40 == 40,38,36 iters (320,342,396 mults)
  std::vector<int> ords({7}); // Nbasis 40 == 40 iters (320 mults)  

  // Standard CG
  //      result=Zero();
  //      CGfine(HermOpEO, src, result);

  for(int l=0;l<los.size();l++){

    RealD lo = los[l];

    for(int o=0;o<ords.size();o++){

      ConjugateGradient<CoarseVector>  CGsloppy(4.0e-2,maxit,false);
      HPDSolver<CoarseVector> HPDSolveSloppy(CoarseOp,CGsloppy,DeflCoarseGuesser);
      
      //    ChebyshevSmoother<LatticeFermionD,HermFineMatrix > Smoother(lo,92,10,FineHermOp); // 36 best case
      ChebyshevSmoother<LatticeFermionD,HermFineMatrix > Smoother(lo,92,ords[o],FineHermOp);  // 311

      //////////////////////////////////////////
      // Build a HDCG solver
      //////////////////////////////////////////
      TwoLevelFlexiblePcg<LatticeFermion,CoarseVector,Subspace>
	HDCG(1.0e-8, 3000,
	     FineHermOp,
	     Smoother,
	     HPDSolveSloppy,
	     HPDSolve,
	     Aggregates);

      //    result=Zero();
      //    HDCG(src,result);
    
      result=Zero();
      HDCG.Inflexible(src,result);
    }
  }
  
  Grid_finalize();
  return 0;
}
