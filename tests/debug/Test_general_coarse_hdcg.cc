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
//#include <Grid/algorithms/GeneralCoarsenedMatrix.h>
#include <Grid/algorithms/iterative/AdefGeneric.h>

using namespace std;
using namespace Grid;

template<class Coarsened>
void SaveOperator(Coarsened &Operator,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(Operator.Grid()->IsBoss());
  assert(Operator._A.size()==Operator.geom.npoint);
  WR.open(file);
  for(int p=0;p<Operator._A.size();p++){
    auto tmp = Operator.Cell.Extract(Operator._A[p]);
    WR.writeScidacFieldRecord(tmp,record);
  }
  WR.close();
#endif
}
template<class Coarsened>
void LoadOperator(Coarsened &Operator,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  Grid::ScidacReader RD ;
  RD.open(file);
  assert(Operator._A.size()==Operator.geom.npoint);
  for(int p=0;p<Operator.geom.npoint;p++){
    conformable(Operator._A[p].Grid(),Operator.CoarseGrid());
    RD.readScidacFieldRecord(Operator._A[p],record);
  }    
  RD.close();
  Operator.ExchangeCoarseLinks();
#endif
}
template<class aggregation>
void SaveBasis(aggregation &Agg,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(Agg.FineGrid->IsBoss());
  WR.open(file);
  for(int b=0;b<Agg.subspace.size();b++){
    WR.writeScidacFieldRecord(Agg.subspace[b],record);
  }
  WR.close();
#endif
}
template<class aggregation>
void LoadBasis(aggregation &Agg, std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacReader RD ;
  RD.open(file);
  for(int b=0;b<Agg.subspace.size();b++){
    RD.readScidacFieldRecord(Agg.subspace[b],record);
  }    
  RD.close();
#endif
}


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

  const int Ls=24;
  const int nbasis = 40;
  const int cb = 0 ;
  RealD mass=0.00078;
  RealD M5=1.8;
  RealD b=1.5;
  RealD c=0.5;

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
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNextToNextToNearestStencilGeometry5D geom(Coarse5d);
  NearestStencilGeometry5D geom_nn(Coarse5d);
  
  // Warning: This routine calls PVdagM.Op, not PVdagM.HermOp
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FrbGrid,cb);

  ////////////////////////////////////////////////////////////
  // Need to check about red-black grid coarsening
  ////////////////////////////////////////////////////////////
  LittleDiracOperator LittleDiracOp(geom,FrbGrid,Coarse5d);

  bool load=false;
  if ( load ) {
    LoadBasis(Aggregates,"/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.scidac");
    LoadOperator(LittleDiracOp,"/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/LittleDiracOp.scidac");
  } else {
    Aggregates.CreateSubspaceChebyshev(RNG5,HermOpEO,nbasis,
				       95.0,0.1,
				       //				     400,200,200 -- 48 iters
				       //				     600,200,200 -- 38 iters, 162s
				       //				     600,200,100 -- 38 iters, 169s
				       //				     600,200,50  -- 88 iters. 370s 
				       800,
				       200,
				       100,
				       0.0);
    LittleDiracOp.CoarsenOperator(FineHermOp,Aggregates);
    SaveBasis(Aggregates,"/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.scidac");
    SaveOperator(LittleDiracOp,"/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/LittleDiracOp.scidac");
  }
  
  // Try projecting to one hop only
  LittleDiracOperator LittleDiracOpProj(geom_nn,FrbGrid,Coarse5d);
  LittleDiracOpProj.ProjectNearestNeighbour(0.01,LittleDiracOp); // smaller shift 0.02? n

  typedef HermitianLinearOperator<LittleDiracOperator,CoarseVector> HermMatrix;
  HermMatrix CoarseOp     (LittleDiracOp);
  HermMatrix CoarseOpProj (LittleDiracOpProj);
  
  //////////////////////////////////////////
  // Build a coarse lanczos
  //////////////////////////////////////////
  Chebyshev<CoarseVector>      IRLCheby(0.2,40.0,71);  // 1 iter
  FunctionHermOp<CoarseVector> IRLOpCheby(IRLCheby,CoarseOp);
  PlainHermOp<CoarseVector>    IRLOp    (CoarseOp);
  int Nk=48;
  int Nm=64;
  int Nstop=Nk;
  ImplicitlyRestartedLanczos<CoarseVector> IRL(IRLOpCheby,IRLOp,Nstop,Nk,Nm,1.0e-5,20);

  int Nconv;
  std::vector<RealD>            eval(Nm);
  std::vector<CoarseVector>     evec(Nm,Coarse5d);
  CoarseVector c_src(Coarse5d);
  //c_src=1.0;
  random(CRNG,c_src);

  CoarseVector c_res(Coarse5d); 
  CoarseVector c_ref(Coarse5d); 

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
  c_res=Zero();
  HPDSolve(c_src,c_res); c_ref = c_res;
  std::cout << GridLogMessage<<"src norm "<<norm2(c_src)<<std::endl;
  std::cout << GridLogMessage<<"ref norm "<<norm2(c_ref)<<std::endl;
  //////////////////////////////////////////////////////////////////////////
  // Deflated (with real op EV's) solve for the projected coarse op
  // Work towards ADEF1 in the coarse space
  //////////////////////////////////////////////////////////////////////////
  HPDSolver<CoarseVector> HPDSolveProj(CoarseOpProj,CG,DeflCoarseGuesser);
  c_res=Zero();
  HPDSolveProj(c_src,c_res);
  std::cout << GridLogMessage<<"src norm "<<norm2(c_src)<<std::endl;
  std::cout << GridLogMessage<<"res norm "<<norm2(c_res)<<std::endl;
  c_res = c_res - c_ref;
  std::cout << "Projected solver error "<<norm2(c_res)<<std::endl;

  //////////////////////////////////////////////////////////////////////
  // Coarse ADEF1 with deflation space
  //////////////////////////////////////////////////////////////////////
  ChebyshevSmoother<CoarseVector,HermMatrix >
    CoarseSmoother(1.0,37.,8,CoarseOpProj);  // just go to sloppy 0.1 convergence
    //  CoarseSmoother(0.1,37.,8,CoarseOpProj);  //
  //  CoarseSmoother(0.5,37.,6,CoarseOpProj);  //  8 iter 0.36s
  //    CoarseSmoother(0.5,37.,12,CoarseOpProj);  // 8 iter, 0.55s
  //    CoarseSmoother(0.5,37.,8,CoarseOpProj);// 7-9 iter
  //  CoarseSmoother(1.0,37.,8,CoarseOpProj); // 0.4 - 0.5s solve to 0.04, 7-9 iter
  //  ChebyshevSmoother<CoarseVector,HermMatrix > CoarseSmoother(0.5,36.,10,CoarseOpProj);  // 311

  ////////////////////////////////////////////////////////
  // CG, Cheby mode spacing 200,200
  // Unprojected Coarse CG solve to 1e-8 : 190 iters, 4.9s
  // Unprojected Coarse CG solve to 4e-2 :  33 iters, 0.8s
  // Projected Coarse CG solve to 1e-8 : 100 iters, 0.36s
  ////////////////////////////////////////////////////////
  // CoarseSmoother(1.0,48.,8,CoarseOpProj); 48 evecs 
  ////////////////////////////////////////////////////////
  // ADEF1 Coarse solve to 1e-8 : 44 iters, 2.34s  2.1x gain
  // ADEF1 Coarse solve to 4e-2 : 7 iters, 0.4s
  // HDCG 38 iters 162s
  //
  // CoarseSmoother(1.0,40.,8,CoarseOpProj); 48 evecs 
  // ADEF1 Coarse solve to 1e-8 : 37 iters, 2.0s  2.1x gain
  // ADEF1 Coarse solve to 4e-2 : 6 iters, 0.36s
  // HDCG 38 iters 169s

  TwoLevelADEF1defl<CoarseVector>
    cADEF1(1.0e-8, 500,
	   CoarseOp,
	   CoarseSmoother,
	   evec,eval);

  c_res=Zero();
  cADEF1(c_src,c_res);
  std::cout << GridLogMessage<<"src norm "<<norm2(c_src)<<std::endl;
  std::cout << GridLogMessage<<"cADEF1 res norm "<<norm2(c_res)<<std::endl;
  c_res = c_res - c_ref;
  std::cout << "cADEF1 solver error "<<norm2(c_res)<<std::endl;
  
  //  cADEF1.Tolerance = 4.0e-2;
  //  cADEF1.Tolerance = 1.0e-1;
  cADEF1.Tolerance = 5.0e-2;
  c_res=Zero();
  cADEF1(c_src,c_res);
  std::cout << GridLogMessage<<"src norm "<<norm2(c_src)<<std::endl;
  std::cout << GridLogMessage<<"cADEF1 res norm "<<norm2(c_res)<<std::endl;
  c_res = c_res - c_ref;
  std::cout << "cADEF1 solver error "<<norm2(c_res)<<std::endl;
  
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
      TwoLevelADEF2<LatticeFermion,CoarseVector,Subspace>
	HDCG(1.0e-8, 100,
	     FineHermOp,
	     Smoother,
	     HPDSolveSloppy,
	     HPDSolve,
	     Aggregates);

      TwoLevelADEF2<LatticeFermion,CoarseVector,Subspace>
	HDCGdefl(1.0e-8, 100,
		 FineHermOp,
		 Smoother,
		 cADEF1,
		 HPDSolve,
		 Aggregates);
      
      result=Zero();
      HDCGdefl(src,result);

      result=Zero();
      HDCG(src,result);

      
    }
  }

  // Standard CG
  result=Zero();
  CGfine(HermOpEO, src, result);
  
  Grid_finalize();
  return 0;
}
