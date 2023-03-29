/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_hdcr.cc

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>

using namespace std;
using namespace Grid;
/* Params
 * Grid: 
 * block1(4)
 * block2(4)
 * 
 * Subspace
 * * Fine  : Subspace(nbasis,hi,lo,order,first,step) -- 32, 60,0.02,500,100,100
 * * Coarse: Subspace(nbasis,hi,lo,order,first,step) -- 32, 18,0.02,500,100,100

 * Smoother:
 * * Fine: Cheby(hi, lo, order)            --  60,0.5,10
 * * Coarse: Cheby(hi, lo, order)          --  12,0.1,4

 * Lanczos:
 * CoarseCoarse IRL( Nk, Nm, Nstop, poly(lo,hi,order))   24,36,24,0.002,4.0,61 
 */
RealD InverseApproximation(RealD x){
  return 1.0/x;
}

template<class Field,class Matrix> class ChebyshevSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field>                            FineOperator;
  Matrix         & _SmootherMatrix;
  FineOperator   & _SmootherOperator;
  
  Chebyshev<Field> Cheby;

  ChebyshevSmoother(RealD _lo,RealD _hi,int _ord, FineOperator &SmootherOperator,Matrix &SmootherMatrix) :
    _SmootherOperator(SmootherOperator),
    _SmootherMatrix(SmootherMatrix),
    Cheby(_lo,_hi,_ord,InverseApproximation)
  {};

  void operator() (const Field &in, Field &out) 
  {
    Field tmp(in.Grid());
    MdagMLinearOperator<Matrix,Field>   MdagMOp(_SmootherMatrix); 
    _SmootherOperator.AdjOp(in,tmp);
    Cheby(MdagMOp,tmp,out);         
  }
};
template<class Field,class Matrix> class MirsSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field>                            FineOperator;
  Matrix         & SmootherMatrix;
  FineOperator   & SmootherOperator;
  RealD tol;
  RealD shift;
  int   maxit;

  MirsSmoother(RealD _shift,RealD _tol,int _maxit,FineOperator &_SmootherOperator,Matrix &_SmootherMatrix) :
    shift(_shift),tol(_tol),maxit(_maxit),
    SmootherOperator(_SmootherOperator),
    SmootherMatrix(_SmootherMatrix)
  {};

  void operator() (const Field &in, Field &out) 
  {
    ZeroGuesser<Field> Guess;
    ConjugateGradient<Field>  CG(tol,maxit,false);
 
    Field src(in.Grid());

    ShiftedMdagMLinearOperator<SparseMatrixBase<Field>,Field> MdagMOp(SmootherMatrix,shift);
    SmootherOperator.AdjOp(in,src);
    Guess(src,out);
    CG(MdagMOp,src,out); 
  }
};

template<class Fobj,class CComplex,int nbasis, class Matrix, class Guesser, class CoarseSolver>
class MultiGridPreconditioner : public LinearFunction< Lattice<Fobj> > {
public:
  using LinearFunction<Lattice<Fobj> >::operator();

  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;
  typedef CoarsenedMatrix<Fobj,CComplex,nbasis> CoarseOperator;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseMatrix CoarseMatrix;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::FineField    FineField;
  typedef LinearOperatorBase<FineField>                            FineOperator;
  typedef LinearFunction    <FineField>                            FineSmoother;

  Aggregates     & _Aggregates;
  CoarseOperator & _CoarseOperator;
  Matrix         & _FineMatrix;
  FineOperator   & _FineOperator;
  Guesser        & _Guess;
  FineSmoother   & _Smoother;
  CoarseSolver   & _CoarseSolve;

  int    level;  void Level(int lv) {level = lv; };

#define GridLogLevel std::cout << GridLogMessage <<std::string(level,'\t')<< " Level "<<level <<" "

  MultiGridPreconditioner(Aggregates &Agg, CoarseOperator &Coarse, 
			  FineOperator &Fine,Matrix &FineMatrix,
			  FineSmoother &Smoother,
			  Guesser &Guess_,
			  CoarseSolver &CoarseSolve_)
    : _Aggregates(Agg),
      _CoarseOperator(Coarse),
      _FineOperator(Fine),
      _FineMatrix(FineMatrix),
      _Smoother(Smoother),
      _Guess(Guess_),
      _CoarseSolve(CoarseSolve_),
      level(1)  {  }

  virtual void operator()(const FineField &in, FineField & out) 
  {
    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid()); 
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    double t;
    // Fine Smoother
    t=-usecond();
    _Smoother(in,out);
    t+=usecond();
    GridLogLevel << "Smoother took "<< t/1000.0<< "ms" <<std::endl;

    // Update the residual
    _FineOperator.Op(out,vec1);  sub(vec1, in ,vec1);   

    // Fine to Coarse 
    t=-usecond();
    _Aggregates.ProjectToSubspace  (Csrc,vec1);
    t+=usecond();
    GridLogLevel << "Project to coarse took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse correction
    t=-usecond();
    _CoarseSolve(Csrc,Csol);
    t+=usecond();
    GridLogLevel << "Coarse solve took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse to Fine
    t=-usecond();
    _Aggregates.PromoteFromSubspace(Csol,vec1); 
    add(out,out,vec1);
    t+=usecond();
    GridLogLevel << "Promote to this level took "<< t/1000.0<< "ms" <<std::endl;

    // Residual
    _FineOperator.Op(out,vec1);  sub(vec1 ,in , vec1);  

    // Fine Smoother
    t=-usecond();
    _Smoother(vec1,vec2);
    t+=usecond();
    GridLogLevel << "Smoother took "<< t/1000.0<< "ms" <<std::endl;

    add( out,out,vec2);
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=24;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  ///////////////////////////////////////////////////
  // Construct a coarsened grid; utility for this?
  ///////////////////////////////////////////////////
  std::vector<int> block ({2,2,2,2});
  std::vector<int> blockc ({2,2,2,2});
  const int nbasis= 40;
  const int nbasisc= 40;
  auto clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/block[d];
  }
  auto cclatt = clatt;
  for(int d=0;d<clatt.size();d++){
    cclatt[d] = clatt[d]/blockc[d];
  }

  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);
  //  GridCartesian *CoarseCoarse4d =  SpaceTimeGrid::makeFourDimGrid(cclatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  //  GridCartesian *CoarseCoarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);
  LatticeFermion    src(FGrid); gaussian(RNG5,src);// src=src+g5*src;
  LatticeFermion result(FGrid); 
  LatticeGaugeField Umu(UGrid); 

  FieldMetaData header;
  //  std::string file("./ckpoint_lat.4000");
  std::string file("./ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building g5R5 hermitian DWF operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  RealD mass=0.00078;
  RealD M5=1.8;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>              Subspace;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>          CoarseOperator;
  typedef CoarseOperator::CoarseVector                                 CoarseVector;
  typedef CoarseOperator::siteVector siteVector;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling Aggregation class to build subspace" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  MdagMLinearOperator<DomainWallFermionR,LatticeFermion> HermDefOp(Ddwf);

  Subspace Aggregates(Coarse5d,FGrid,0);

  assert ( (nbasis & 0x1)==0);
  {
    int nb=nbasis/2;
    LatticeFermion A(FGrid);
    LatticeFermion B(FGrid);
    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.002,1000,800,100,0.0);
    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.02,1000,800,100,0.0); 
    Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.01,1000,100,100,0.0); // Slightly faster

    for(int n=0;n<nb;n++){
      std::cout << GridLogMessage << " G5R5 "<<n<<std::endl;
      G5R5(Aggregates.subspace[n+nb],Aggregates.subspace[n]);
      std::cout << GridLogMessage << " Projection "<<n<<std::endl;
      A = Aggregates.subspace[n];
      B = Aggregates.subspace[n+nb];
      std::cout << GridLogMessage << " Copy "<<n<<std::endl;
      Aggregates.subspace[n]   = A+B; // 1+G5 // eigen value of G5R5 is +1
      std::cout << GridLogMessage << " P+ "<<n<<std::endl;
      Aggregates.subspace[n+nb]= A-B; // 1-G5 // eigen value of G5R5 is -1
      std::cout << GridLogMessage << " P- "<<n<<std::endl;
    }
  }
  
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building coarse representation of Indef operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>    Level1Op;
  typedef CoarsenedMatrix<siteVector,iScalar<vTComplex>,nbasisc> Level2Op;

  Gamma5R5HermitianLinearOperator<DomainWallFermionR,LatticeFermion> HermIndefOp(Ddwf);

  
  GridRedBlackCartesian * Coarse4dRB = SpaceTimeGrid::makeFourDimRedBlackGrid(Coarse4d);
  std::cout << " Making 5D coarse RB grid " <<std::endl;
  GridRedBlackCartesian * Coarse5dRB = SpaceTimeGrid::makeFiveDimRedBlackGrid(1,Coarse4d);
  std::cout << " Made 5D coarse RB grid " <<std::endl;
  Level1Op LDOp(*Coarse5d,*Coarse5dRB,1); LDOp.CoarsenOperator(FGrid,HermIndefOp,Aggregates);


  //////////////////////////////////////////////////
  // Deflate the course space. Recursive multigrid?
  //////////////////////////////////////////////////
  typedef Aggregation<siteVector,iScalar<vTComplex>,nbasisc>                   CoarseSubspace;
  //  CoarseSubspace CoarseAggregates(CoarseCoarse5d,Coarse5d,0);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Build deflation space in coarse operator "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  MdagMLinearOperator<CoarseOperator,CoarseVector> PosdefLdop(LDOp);
  /*
  {
    int nb=nbasisc/2;
    CoarseAggregates.CreateSubspaceChebyshev(CRNG,PosdefLdop,nb,15.0,0.02,1000,800,100,0.0);
    for(int n=0;n<nb;n++){
      autoView( subspace   , CoarseAggregates.subspace[n],CpuWrite);
      autoView( subspace_g5, CoarseAggregates.subspace[n+nb],CpuWrite);
      for(int nn=0;nn<nb;nn++){
	for(int site=0;site<Coarse5d->oSites();site++){
	  subspace_g5[site](nn)   = subspace[site](nn);
	  subspace_g5[site](nn+nb)=-subspace[site](nn+nb);
	}
      }
    }
  }
  */
  typedef Level2Op::CoarseVector CoarseCoarseVector;
  /*
  Level2Op L2Op(*CoarseCoarse5d,1); // Hermitian matrix
  HermitianLinearOperator<Level1Op,CoarseVector> L1LinOp(LDOp);
  L2Op.CoarsenOperator(Coarse5d,L1LinOp,CoarseAggregates);


  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Running CoarseCoarse grid Lanczos "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  MdagMLinearOperator<Level2Op,CoarseCoarseVector> IRLHermOpL2(L2Op);
  CoarseCoarseVector cc_src(CoarseCoarse5d); cc_src=1.0;
  */
  /*
  Chebyshev<CoarseCoarseVector> IRLChebyL2(0.001,15.0,301);
  FunctionHermOp<CoarseCoarseVector> IRLOpChebyL2(IRLChebyL2,IRLHermOpL2);
  PlainHermOp<CoarseCoarseVector> IRLOpL2    (IRLHermOpL2);
  int cNk=24;
  int cNm=36;
  int cNstop=24;
  ImplicitlyRestartedLanczos<CoarseCoarseVector> IRLL2(IRLOpChebyL2,IRLOpL2,cNstop,cNk,cNm,1.0e-3,20);

  int cNconv;
  std::vector<RealD>          eval2(cNm);
  std::vector<CoarseCoarseVector>   evec2(cNm,CoarseCoarse5d);
  IRLL2.calc(eval2,evec2,cc_src,cNconv);

  ConjugateGradient<CoarseCoarseVector>  CoarseCoarseCG(0.1,1000);
  DeflatedGuesser<CoarseCoarseVector> DeflCoarseCoarseGuesser(evec2,eval2);
  NormalEquations<CoarseCoarseVector> DeflCoarseCoarseCGNE(L2Op,CoarseCoarseCG,DeflCoarseCoarseGuesser);
  */

  /*
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Running Coarse grid Lanczos "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  MdagMLinearOperator<Level1Op,CoarseVector> IRLHermOp(LDOp);
  //  Chebyshev<CoarseVector>      IRLCheby(0.001,15.0,301);
  Chebyshev<CoarseVector>      IRLCheby(0.03,12.0,101);
  FunctionHermOp<CoarseVector> IRLOpCheby(IRLCheby,IRLHermOp);
  PlainHermOp<CoarseVector>    IRLOp    (IRLHermOp);
  int Nk=64;
  int Nm=128;
  int Nstop=Nk;
  ImplicitlyRestartedLanczos<CoarseVector> IRL(IRLOpCheby,IRLOp,Nstop,Nk,Nm,1.0e-3,20);

  int Nconv;
  std::vector<RealD>            eval(Nm);
  std::vector<CoarseVector>     evec(Nm,Coarse5d);
  IRL.calc(eval,evec,c_src,Nconv);
  */
  CoarseVector c_src(Coarse5d); c_src=1.0;
  //  DeflatedGuesser<CoarseVector> DeflCoarseGuesser(evec,eval);
  //  NormalEquations<CoarseVector> DeflCoarseCGNE(LDOp,CoarseCG,DeflCoarseGuesser);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building 3 level Multigrid            "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  //  typedef MultiGridPreconditioner<vSpinColourVector,  vTComplex,nbasis, DomainWallFermionR,DeflatedGuesser<CoarseVector> , NormalEquations<CoarseVector> >   TwoLevelMG;
  typedef MultiGridPreconditioner<vSpinColourVector,  vTComplex,nbasis, DomainWallFermionR,ZeroGuesser<CoarseVector> , NormalEquations<CoarseVector> >   TwoLevelMG;
  typedef MultiGridPreconditioner<siteVector,iScalar<vTComplex>,nbasisc,Level1Op, DeflatedGuesser<CoarseCoarseVector>, NormalEquations<CoarseCoarseVector> > CoarseMG;
  typedef MultiGridPreconditioner<vSpinColourVector,  vTComplex,nbasis, DomainWallFermionR,ZeroGuesser<CoarseVector>, LinearFunction<CoarseVector> >     ThreeLevelMG;

  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(0.25,60.0,12,HermIndefOp,Ddwf);
  /*
  // MultiGrid preconditioner acting on the coarse space <-> coarsecoarse space
  ChebyshevSmoother<CoarseVector,  Level1Op >        CoarseSmoother(0.1,15.0,3,L1LinOp,LDOp);

  //  MirsSmoother<CoarseVector,  Level1Op >        CoarseCGSmoother(0.1,0.1,4,L1LinOp,LDOp);
  //  MirsSmoother<LatticeFermion,DomainWallFermionR> FineCGSmoother(0.0,0.01,8,HermIndefOp,Ddwf);

  CoarseMG Level2Precon (CoarseAggregates, L2Op,
			 L1LinOp,LDOp,
			 CoarseSmoother,
			 DeflCoarseCoarseGuesser,	
		 DeflCoarseCoarseCGNE);
  Level2Precon.Level(2);

  // PGCR Applying this solver to solve the coarse space problem
  PrecGeneralisedConjugateResidual<CoarseVector>  l2PGCR(0.1, 100, L1LinOp,Level2Precon,16,16);
  l2PGCR.Level(2);
  
  // Wrap the 2nd level solver in a MultiGrid preconditioner acting on the fine space
  ZeroGuesser<CoarseVector> CoarseZeroGuesser;
  ThreeLevelMG ThreeLevelPrecon(Aggregates, LDOp,
				HermIndefOp,Ddwf,
				FineSmoother,
				CoarseZeroGuesser,
				l2PGCR);
  ThreeLevelPrecon.Level(1);

  // Apply the fine-coarse-coarsecoarse 2 deep MG preconditioner in an outer PGCR on the fine fgrid
  PrecGeneralisedConjugateResidual<LatticeFermion> l1PGCR(1.0e-8,1000,HermIndefOp,ThreeLevelPrecon,16,16);
  l1PGCR.Level(1);
  */
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling 2 level Multigrid            "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  result=Zero();


    ZeroGuesser<CoarseVector> CoarseZeroGuesser;
    ConjugateGradient<CoarseVector>  CoarseCG(0.01,1000);
    NormalEquations<CoarseVector> CoarseCGNE(LDOp,CoarseCG,CoarseZeroGuesser);
    TwoLevelMG TwoLevelPrecon(Aggregates, LDOp,
			      HermIndefOp,Ddwf,
			      FineSmoother,
			      CoarseZeroGuesser,	
			      CoarseCGNE);
    TwoLevelPrecon.Level(1);
    PrecGeneralisedConjugateResidual<LatticeFermion> l1PGCR(1.0e-8,20,HermIndefOp,TwoLevelPrecon,16,16);
    l1PGCR.Level(1);
    l1PGCR(src,result);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling CG            "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  ConjugateGradient<LatticeFermion> pCG(1.0e-8,60000);
  result=Zero();
  //  pCG(HermDefOp,src,result);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling red black CG            "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  result=Zero();

    LatticeFermion    src_o(FrbGrid);
    LatticeFermion result_o(FrbGrid);
    pickCheckerboard(Odd,src_o,src);
    result_o=Zero();
    SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> HermOpEO(Ddwf);
    pCG(HermOpEO,src_o,result_o);
  
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Fine        PowerMethod           "<< std::endl;
  PowerMethod<LatticeFermion>       PM;   PM(HermDefOp,src);
  std::cout<<GridLogMessage << " Coarse       PowerMethod           "<< std::endl;
  PowerMethod<CoarseVector>        cPM;  cPM(PosdefLdop,c_src);
  //  std::cout<<GridLogMessage << " CoarseCoarse PowerMethod           "<< std::endl;
  //  PowerMethod<CoarseCoarseVector> ccPM; ccPM(IRLHermOpL2,cc_src);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Grid_finalize();
}
