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
/*
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
*/
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

  //////////////////////////////////////////
  // Build a coarse lanczos
  //////////////////////////////////////////
  typedef HermitianLinearOperator<LittleDiracOperator,CoarseVector> HermMatrix;
  HermMatrix CoarseOp     (LittleDiracOp);
  
  int Nk=192;
  int Nm=256;
  int Nstop=Nk;
  
  Chebyshev<CoarseVector>      IRLCheby(0.005,40.0,201);
  //  Chebyshev<CoarseVector>      IRLCheby(0.010,45.0,201);  // 1 iter
  FunctionHermOp<CoarseVector> IRLOpCheby(IRLCheby,CoarseOp);
  PlainHermOp<CoarseVector>    IRLOp    (CoarseOp);
  
  ImplicitlyRestartedLanczos<CoarseVector> IRL(IRLOpCheby,IRLOp,Nstop,Nk,Nm,1e-5,10);

  int Nconv;
  std::vector<RealD>            eval(Nm);
  std::vector<CoarseVector>     evec(Nm,Coarse5d);

  PowerMethod<CoarseVector>       cPM;   cPM(CoarseOp,c_src);

  IRL.calc(eval,evec,c_src,Nconv);

  //////////////////////////////////////////
  // Deflated guesser
  //////////////////////////////////////////
  DeflatedGuesser<CoarseVector> DeflCoarseGuesser(evec,eval);

  int maxit=30000;
  ConjugateGradient<CoarseVector>  CG(1.0e-10,maxit,false);
  ConjugateGradient<LatticeFermionD>  CGfine(1.0e-8,30000,false);

  //////////////////////////////////////////
  // HDCG
  //////////////////////////////////////////
  
  std::vector<RealD> los({2.0,2.5}); // Nbasis 40 == 36,36 iters
  std::vector<int> ords({9}); // Nbasis 40 == 40 iters (320 mults)  

  for(int l=0;l<los.size();l++){

    RealD lo = los[l];

    for(int o=0;o<ords.size();o++){

      //////////////////////////////////////////
      // Sloppy coarse solve
      //////////////////////////////////////////
      
      ConjugateGradient<CoarseVector>  CGsloppy(4.0e-2,maxit,false);
      HPDSolver<CoarseVector> HPDSolveSloppy(CoarseOp,CGsloppy,DeflCoarseGuesser);
      HPDSolver<CoarseVector> HPDSolve(CoarseOp,CG,DeflCoarseGuesser);

      //////////////////////////////////////////
      // IRS shifted smoother based on CG
      //////////////////////////////////////////
      RealD MirsShift = lo;
      ShiftedHermOpLinearOperator<LatticeFermionD> ShiftedFineHermOp(HermOpEO,MirsShift);
      CGSmoother<LatticeFermionD> CGsmooth(ords[o],ShiftedFineHermOp) ;
  
      //////////////////////////////////////////
      // Build a HDCG solver
      //////////////////////////////////////////
      TwoLevelADEF2<LatticeFermion,CoarseVector,Subspace>
	HDCG(1.0e-8, 700,
	     FineHermOp,
	     CGsmooth,
	     HPDSolveSloppy,
	     HPDSolve,
	     Aggregates);

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
