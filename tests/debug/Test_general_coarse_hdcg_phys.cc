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
    WR.writeScidacFieldRecord(tmp,record,0,0);
    //    WR.writeScidacFieldRecord(tmp,record,0,BINARYIO_LEXICOGRAPHIC);
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
    //    RD.readScidacFieldRecord(Operator._A[p],record,BINARYIO_LEXICOGRAPHIC);
    RD.readScidacFieldRecord(Operator._A[p],record,0);
  }    
  RD.close();
  Operator.ExchangeCoarseLinks();
#endif
}
template<class Coarsened>
void ReLoadOperator(Coarsened &Operator,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  Grid::ScidacReader RD ;
  RD.open(file);
  assert(Operator._A.size()==Operator.geom.npoint);
  for(int p=0;p<Operator.geom.npoint;p++){
    auto tmp=Operator.Cell.Extract(Operator._A[p]);
    RD.readScidacFieldRecord(tmp,record,0);
    Operator._A[p] = Operator.Cell.ExchangePeriodic(tmp);
  }    
  RD.close();
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
    //WR.writeScidacFieldRecord(Agg.subspace[b],record,0,BINARYIO_LEXICOGRAPHIC);
    WR.writeScidacFieldRecord(Agg.subspace[b],record,0,0);
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
    //    RD.readScidacFieldRecord(Agg.subspace[b],record,BINARYIO_LEXICOGRAPHIC);
    RD.readScidacFieldRecord(Agg.subspace[b],record,0);
  }    
  RD.close();
#endif
}

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
template<class Field> class ChebyshevSmoother : public LinearFunction<Field>
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

template<class Field> class CGSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _SmootherOperator;
  int iters;
  CGSmoother(int _iters, FineOperator &SmootherOperator) :
    _SmootherOperator(SmootherOperator),
    iters(_iters)
  {
    std::cout << GridLogMessage<<" Mirs smoother order "<<iters<<std::endl;
  };
  void operator() (const Field &in, Field &out) 
  {
    ConjugateGradient<Field>  CG(0.0,iters,false); // non-converge is just fine in a smoother
    CG(_SmootherOperator,in,out);
  }
};


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=24;
  const int nbasis = 62;
  //  const int nbasis = 56;
  //  const int nbasis = 44;
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
  Coordinate Block({4,4,6,4});
  Coordinate clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/Block[d];
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
  std::string file("ckpoint_lat.1000");
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

  std::string subspace_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.rat.scidac.62");
  std::string refine_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Refine.phys48.rat.scidac.62");
  std::string ldop_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/LittleDiracOp.phys48.rat.scidac.62");
  bool load_agg=true;
  bool load_refine=true;
  bool load_mat=true;
  if ( load_agg ) {
    LoadBasis(Aggregates,subspace_file);
  } else {

    // NBASIS=40
    // Best so far: ord 2000 [0.01,95], 500,500  -- 466 iters
    // slurm-398626.out:Grid : Message : 141.295253 s : 500 filt [1] <n|MdagM|n> 0.000103622063


    //Grid : Message : 33.870465 s :  Chebyshev subspace pass-1 : ord 2000 [0.001,95]
    //Grid : Message : 33.870485 s :  Chebyshev subspace pass-2 : nbasis40 min 1000 step 1000 lo0
    //slurm-1482200.out : filt ~ 0.004 -- not as low mode projecting -- took 626 iters

    // To try: 2000 [0.1,95]  ,2000,500,500 -- slurm-1482213.out 586 iterations

    // To try: 2000 [0.01,95] ,2000,500,500 -- 469 (think I bumped 92 to 95) (??)
    // To try: 2000 [0.025,95],2000,500,500
    // To try: 2000 [0.005,95],2000,500,500

    // NBASIS=44 -- HDCG paper was 64 vectors; AMD compiler craps out at 48
    // To try: 2000 [0.01,95] ,2000,500,500 -- 419 lowest slurm-1482355.out
    // To try: 2000 [0.025,95] ,2000,500,500 -- 487 
    // To try: 2000 [0.005,95] ,2000,500,500
    /*
      Smoother [3,92] order 16
slurm-1482355.out:Grid : Message : 35.239686 s :  Chebyshev subspace pass-1 : ord 2000 [0.01,95]
slurm-1482355.out:Grid : Message : 35.239714 s :  Chebyshev subspace pass-2 : nbasis44 min 500 step 500 lo0
slurm-1482355.out:Grid : Message : 5561.305552 s : HDCG: Pcg converged in 419 iterations and 2616.202598 s

slurm-1482367.out:Grid : Message : 43.157235 s :  Chebyshev subspace pass-1 : ord 2000 [0.025,95]
slurm-1482367.out:Grid : Message : 43.157257 s :  Chebyshev subspace pass-2 : nbasis44 min 500 step 500 lo0
slurm-1482367.out:Grid : Message : 6169.469330 s : HDCG: Pcg converged in 487 iterations and 3131.185821 s
    */
		 /*
		   Aggregates.CreateSubspaceChebyshev(RNG5,HermOpEO,nbasis,
				       95.0,0.0075,
				       2500,
				       500,
				       500,
				       0.0);
		 */

		 /*
		   Aggregates.CreateSubspaceChebyshevPowerLaw(RNG5,HermOpEO,nbasis,
							      95.0,
							      2000);
		 */

    Aggregates.CreateSubspaceMultishift(RNG5,HermOpEO,
					0.0003,1.0e-5,2000); // Lo, tol, maxit
  /*
    Aggregates.CreateSubspaceChebyshev(RNG5,HermOpEO,nbasis,
				       95.0,0.05,
				       2000,
				       500,
				       500,
				       0.0);
 */
    /*
      Aggregates.CreateSubspaceChebyshev(RNG5,HermOpEO,nbasis,
				       95.0,0.01,
				       2000,
				       500,
				       500,
				       0.0);
    */
    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermOpEO,nbasis,95.,0.01,1500); -- running slurm-1484934.out nbasis 56

    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermOpEO,nbasis,95.,0.01,1500); <== last run
    SaveBasis(Aggregates,subspace_file);
  }

  int refine=1;
  if(refine){
    if ( load_refine ) {
      LoadBasis(Aggregates,refine_file);
    } else {
      // HDCG used Pcg to refine
      Aggregates.RefineSubspace(HermOpEO,0.001,1.0e-3,3000);
      SaveBasis(Aggregates,refine_file);
    }
  }

  Aggregates.Orthogonalise();
  if ( load_mat ) {
    LoadOperator(LittleDiracOp,ldop_file);
  } else {
    LittleDiracOp.CoarsenOperator(FineHermOp,Aggregates);
    SaveOperator(LittleDiracOp,ldop_file);
  }

  // I/O test:
  CoarseVector c_src(Coarse5d);   random(CRNG,c_src);
  CoarseVector c_res(Coarse5d); 
  CoarseVector c_ref(Coarse5d);

  // Try projecting to one hop only
  //  LittleDiracOp.ShiftMatrix(1.0e-4);
  LittleDiracOperator LittleDiracOpProj(geom_nn,FrbGrid,Coarse5d);
  LittleDiracOpProj.ProjectNearestNeighbour(0.01,LittleDiracOp); // smaller shift 0.02? n

  typedef HermitianLinearOperator<LittleDiracOperator,CoarseVector> HermMatrix;
  HermMatrix CoarseOp     (LittleDiracOp);
  HermMatrix CoarseOpProj (LittleDiracOpProj);
  
  //////////////////////////////////////////
  // Build a coarse lanczos
  //////////////////////////////////////////
  //  Chebyshev<CoarseVector>      IRLCheby(0.012,40.0,201);  //500 HDCG iters
  //  int Nk=512; // Didn't save much
  //  int Nm=640;
  //  int Nstop=400;

  //  Chebyshev<CoarseVector>      IRLCheby(0.005,40.0,201);  //319 HDCG iters @ 128//160 nk.
  //  int Nk=128;
  //  int Nm=160;
  Chebyshev<CoarseVector>      IRLCheby(0.005,40.0,201);  //319 HDCG iters @ 128//160 nk.
  int Nk=192;
  int Nm=256;
  int Nstop=Nk;
  
  //  Chebyshev<CoarseVector>      IRLCheby(0.010,45.0,201);  // 1 iter
  FunctionHermOp<CoarseVector> IRLOpCheby(IRLCheby,CoarseOp);
  PlainHermOp<CoarseVector>    IRLOp    (CoarseOp);
  
  ImplicitlyRestartedLanczos<CoarseVector> IRL(IRLOpCheby,IRLOp,Nstop,Nk,Nm,1e-5,10);

  int Nconv;
  std::vector<RealD>            eval(Nm);
  std::vector<CoarseVector>     evec(Nm,Coarse5d);

  PowerMethod<CoarseVector>       cPM;   cPM(CoarseOp,c_src);

  IRL.calc(eval,evec,c_src,Nconv);
  DeflatedGuesser<CoarseVector> DeflCoarseGuesser(evec,eval);

  //////////////////////////////////////////
  // Build a coarse space solver
  //////////////////////////////////////////
  int maxit=30000;
  ConjugateGradient<CoarseVector>  CG(1.0e-10,maxit,false);
  ConjugateGradient<LatticeFermionD>  CGfine(1.0e-8,30000,false);
  ZeroGuesser<CoarseVector> CoarseZeroGuesser;

  //  HPDSolver<CoarseVector> HPDSolve(CoarseOp,CG,CoarseZeroGuesser);
  HPDSolver<CoarseVector> HPDSolve(CoarseOp,CG,DeflCoarseGuesser);
  c_res=Zero();
  //  HPDSolve(c_src,c_res); c_ref = c_res;
  //  std::cout << GridLogMessage<<"src norm "<<norm2(c_src)<<std::endl;
  //  std::cout << GridLogMessage<<"ref norm "<<norm2(c_ref)<<std::endl;
  //////////////////////////////////////////////////////////////////////////
  // Deflated (with real op EV's) solve for the projected coarse op
  // Work towards ADEF1 in the coarse space
  //////////////////////////////////////////////////////////////////////////
  HPDSolver<CoarseVector> HPDSolveProj(CoarseOpProj,CG,DeflCoarseGuesser);
  c_res=Zero();
  //  HPDSolveProj(c_src,c_res);
  //  std::cout << GridLogMessage<<"src norm "<<norm2(c_src)<<std::endl;
  //  std::cout << GridLogMessage<<"res norm "<<norm2(c_res)<<std::endl;
  //  c_res = c_res - c_ref;
  //  std::cout << "Projected solver error "<<norm2(c_res)<<std::endl;

  //////////////////////////////////////////////////////////////////////
  // Coarse ADEF1 with deflation space
  //////////////////////////////////////////////////////////////////////
  ChebyshevSmoother<CoarseVector >  CoarseSmoother(1.0,37.,8,CoarseOpProj);  // just go to sloppy 0.1 convergence
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

  //  c_res=Zero();
  //  cADEF1(c_src,c_res);
  //  std::cout << GridLogMessage<<"src norm "<<norm2(c_src)<<std::endl;
  //  std::cout << GridLogMessage<<"cADEF1 res norm "<<norm2(c_res)<<std::endl;
  //  c_res = c_res - c_ref;
  //  std::cout << "cADEF1 solver error "<<norm2(c_res)<<std::endl;
  
  //  cADEF1.Tolerance = 4.0e-2;
  //  cADEF1.Tolerance = 1.0e-1;
  //  cADEF1.Tolerance = 5.0e-2;
  //  c_res=Zero();
  //  cADEF1(c_src,c_res);
  //  std::cout << GridLogMessage<<"src norm "<<norm2(c_src)<<std::endl;
  //  std::cout << GridLogMessage<<"cADEF1 res norm "<<norm2(c_res)<<std::endl;
  //  c_res = c_res - c_ref;
  //  std::cout << "cADEF1 solver error "<<norm2(c_res)<<std::endl;
  
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
  //
  //
  //
  
  std::vector<RealD> los({2.0,2.5}); // Nbasis 40 == 36,36 iters

  //  std::vector<int> ords({7,8,10}); // Nbasis 40 == 40,38,36 iters (320,342,396 mults)
  //  std::vector<int> ords({7}); // Nbasis 40 == 40 iters (320 mults)
  std::vector<int> ords({9}); // Nbasis 40 == 40 iters (320 mults)  

 /*
   Smoother opt @56 nbasis, 0.04 convergence, 192 evs
 ord lo

 16   0.1  no converge -- likely sign indefinite
 32   0.1  no converge -- likely sign indefinite(?)

 16   0.5  422
 32   0.5  302
 
 8   1.0  575
 12  1.0  449
 16  1.0  375
 32  1.0  302

 12  3.0  476
 16  3.0  319
 32  3.0  306

 Powerlaw setup 62 vecs
slurm-1494943.out:Grid : Message : 4874.186617 s : HDCG: Pcg converged in 171 iterations and 1706.548006 s 1.0 32
slurm-1494943.out:Grid : Message : 6490.121648 s : HDCG: Pcg converged in 194 iterations and 1616.219654 s 1.0 16

 Cheby setup: 56vecs
 -- CG smoother O(16): 487
 
Power law setup, 56 vecs -- lambda^-5
slurm-1494383.out:Grid : Message : 4377.173265 s : HDCG: Pcg converged in 204 iterations and 1153.548935 s 1.0 32

Power law setup, 56 vecs -- lambda^-3

slurm-1494242.out:Grid : Message : 4370.464814 s : HDCG: Pcg converged in 204 iterations and 1143.494776 s  1.0 32
slurm-1494242.out:Grid : Message : 5432.414820 s : HDCG: Pcg converged in 237 iterations and 1061.455882 s  1.0 16
slurm-1494242.out:Grid : Message : 6588.727977 s : HDCG: Pcg converged in 205 iterations and 1156.565210 s  0.5 32

 Power law setup, 56 vecs -- lambda^-4
 -- CG smoother    O(16): 290
 -- Cheby smoother O(16): 218 -- getting close to the deflation level I expect 169 from BFM paper @O(7) smoother and 64 nbasis

Grid : Message : 2790.797194 s : HDCG: Pcg converged in 190 iterations and 1049.563182 s 1.0 32
Grid : Message : 3766.374396 s : HDCG: Pcg converged in 218 iterations and 975.455668 s  1.0 16
Grid : Message : 4888.746190 s : HDCG: Pcg converged in 191 iterations and 1122.252055 s 0.5 32
Grid : Message : 5956.679661 s : HDCG: Pcg converged in 231 iterations and 1067.812850 s 0.5 16

Grid : Message : 2767.405829 s : HDCG: Pcg converged in 218 iterations and 967.214067 s -- 16
Grid : Message : 3816.165905 s : HDCG: Pcg converged in 251 iterations and 1048.636269 s -- 12
Grid : Message : 5121.206572 s : HDCG: Pcg converged in 318 iterations and 1304.916168 s -- 8

 
[paboyle@login2.crusher debug]$ grep -v Memory slurm-402426.out  | grep converged | grep HDCG -- [1.0,16] cheby
Grid : Message : 5185.521063 s : HDCG: Pcg converged in 377 iterations and 1595.843529 s

[paboyle@login2.crusher debug]$ grep HDCG  slurm-402184.out | grep onver
Grid : Message : 3760.438160 s : HDCG: Pcg converged in 422 iterations and 2129.243141 s
Grid : Message : 5660.588015 s : HDCG: Pcg converged in 308 iterations and 1900.026821 s

 
Grid : Message : 4238.206528 s : HDCG: Pcg converged in 575 iterations and 2657.430676 s
Grid : Message : 6345.880344 s : HDCG: Pcg converged in 449 iterations and 2108.505208 s

grep onverg slurm-401663.out | grep HDCG
Grid : Message : 3900.817781 s : HDCG: Pcg converged in 476 iterations and 1992.591311 s
Grid : Message : 5647.202699 s : HDCG: Pcg converged in 306 iterations and 1746.838660 s


[paboyle@login2.crusher debug]$ grep converged slurm-401775.out | grep HDCG
Grid : Message : 3583.177025 s : HDCG: Pcg converged in 375 iterations and 1800.896037 s
Grid : Message : 5348.342243 s : HDCG: Pcg converged in 302 iterations and 1765.045018 s

Conclusion: higher order smoother is doing better. Much better. Use a Krylov smoother instead Mirs as in BFM version.

 */
				      //
  for(int l=0;l<los.size();l++){

    RealD lo = los[l];

    for(int o=0;o<ords.size();o++){

      ConjugateGradient<CoarseVector>  CGsloppy(4.0e-2,maxit,false);
      HPDSolver<CoarseVector> HPDSolveSloppy(CoarseOp,CGsloppy,DeflCoarseGuesser);
      
      //    ChebyshevSmoother<LatticeFermionD,HermFineMatrix > Smoother(lo,92,10,FineHermOp); // 36 best case
      ChebyshevSmoother<LatticeFermionD > ChebySmooth(lo,95,ords[o],FineHermOp);  // 311

      /*
       * CG smooth 11 iter: 
       slurm-403825.out:Grid : Message : 4369.824339 s : HDCG: fPcg converged in 215 iterations 3.0
       slurm-403908.out:Grid : Message : 3955.897470 s : HDCG: fPcg converged in 236 iterations 1.0
       slurm-404273.out:Grid : Message : 3843.792191 s : HDCG: fPcg converged in 210 iterations 2.0
       * CG smooth 9 iter: 
      */
      //
      RealD MirsShift = lo;
      ShiftedHermOpLinearOperator<LatticeFermionD> ShiftedFineHermOp(HermOpEO,MirsShift);
      CGSmoother<LatticeFermionD> CGsmooth(ords[o],ShiftedFineHermOp) ;
  
      //////////////////////////////////////////
      // Build a HDCG solver
      //////////////////////////////////////////
      TwoLevelADEF2<LatticeFermion,CoarseVector,Subspace>
	HDCG(1.0e-8, 700,
	     FineHermOp,
	     //	     ChebySmooth,
	     CGsmooth,
	     HPDSolveSloppy,
	     HPDSolve,
	     Aggregates);

      /*
	TwoLevelADEF2<LatticeFermion,CoarseVector,Subspace>
	HDCGdefl(1.0e-8, 700,
		 FineHermOp,
		 Smoother,
		 cADEF1,
		 HPDSolve,
		 Aggregates);
      */
      
      //      result=Zero();
      //      HDCGdefl(src,result);

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
