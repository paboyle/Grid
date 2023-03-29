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
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>

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
template<class Field,class Matrix> class RedBlackSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field>                            FineOperator;
  Matrix         & SmootherMatrix;
  RealD tol;
  RealD shift;
  int   maxit;

  RedBlackSmoother(RealD _shift,RealD _tol,int _maxit,Matrix &_SmootherMatrix) :
    shift(_shift),tol(_tol),maxit(_maxit),
    SmootherMatrix(_SmootherMatrix)
  {};

  void operator() (const Field &in, Field &out) 
  {
    std::cout << " Red Black Smootheer "<<norm2(in)<<" " <<norm2(out)<<std::endl;
    ConjugateGradient<Field>  CG(tol,maxit,false);
    out =Zero();
    SchurRedBlackDiagMooeeSolve<Field> RBSolver(CG);
    RBSolver(SmootherMatrix,in,out); 
    std::cout << " Red Black Smootheer "<<norm2(in)<<" " <<norm2(out)<<std::endl;
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
  FineSmoother   & _Smoother1;
  FineSmoother   & _Smoother2;
  CoarseSolver   & _CoarseSolve;

  int    level;  void Level(int lv) {level = lv; };

#define GridLogLevel std::cout << GridLogMessage <<std::string(level,'\t')<< " Level "<<level <<" "

  MultiGridPreconditioner(Aggregates &Agg, CoarseOperator &Coarse, 
			  FineOperator &Fine,Matrix &FineMatrix,
			  FineSmoother &Smoother1,
			  FineSmoother &Smoother2,
			  Guesser &Guess_,
			  CoarseSolver &CoarseSolve_)
    : _Aggregates(Agg),
      _CoarseOperator(Coarse),
      _FineOperator(Fine),
      _FineMatrix(FineMatrix),
      _Smoother1(Smoother1),
      _Smoother2(Smoother2),
      _Guess(Guess_),
      _CoarseSolve(CoarseSolve_),
      level(1)  {  }

  MultiGridPreconditioner(Aggregates &Agg, CoarseOperator &Coarse, 
			  FineOperator &Fine,Matrix &FineMatrix,
			  FineSmoother &Smoother,
			  Guesser &Guess_,
			  CoarseSolver &CoarseSolve_)
    : _Aggregates(Agg),
      _CoarseOperator(Coarse),
      _FineOperator(Fine),
      _FineMatrix(FineMatrix),
      _Smoother1(Smoother),
      _Smoother2(Smoother),
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
    _Smoother1(in,out);
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
    _Smoother2(vec1,vec2);
    t+=usecond();
    GridLogLevel << "Smoother took "<< t/1000.0<< "ms" <<std::endl;

    add( out,out,vec2);
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=16;
  //  const int rLs=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  ///////////////////////////////////////////////////
  // Construct a coarsened grid; utility for this?
  ///////////////////////////////////////////////////
  std::vector<int> block ({2,2,2,2});
  const int nbasis= 32;

  auto clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/block[d];
  }

  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);
  GridRedBlackCartesian * Coarse4dRB = SpaceTimeGrid::makeFourDimRedBlackGrid(Coarse4d);
  GridRedBlackCartesian * Coarse5dRB = SpaceTimeGrid::makeFiveDimRedBlackGrid(1,Coarse4d);

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
  std::string file("./ckpoint_lat");
  NerscIO::readConfiguration(Umu,header,file);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building g5R5 hermitian DWF operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  RealD mass=0.001;
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
    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.05,500,200,100,0.0);// 18s
    //   rAggregates.CreateSubspaceChebyshev(RNG5,rHermDefOp,nb,60.0,0.05,500,200,150,0.0);// 15.7 23iter
    Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.05,500,200,150,0.0);//
    // pad out the rAggregates.

    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.05,500,500,150,0.0);// 19s 

    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.05,500,200,200,0.0); 15.2s
    //    Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.05,500,500,200,0.0); 16.3s

    for(int n=0;n<nb;n++){
      G5R5(Aggregates.subspace[n+nb],Aggregates.subspace[n]);
    }
    LatticeFermion A(FGrid);
    LatticeFermion B(FGrid);
    for(int n=0;n<nb;n++){
      A = Aggregates.subspace[n];
      B = Aggregates.subspace[n+nb];
      Aggregates.subspace[n]   = A+B; // 1+G5 // eigen value of G5R5 is +1
      Aggregates.subspace[n+nb]= A-B; // 1-G5 // eigen value of G5R5 is -1
    }
  }
  
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building coarse representation of Indef operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>    Level1Op;

  Gamma5R5HermitianLinearOperator<DomainWallFermionR,LatticeFermion> HermIndefOp(Ddwf);

  Level1Op LDOp(*Coarse5d,*Coarse5dRB,1); LDOp.CoarsenOperator(FGrid,HermIndefOp,Aggregates);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Running Coarse grid Lanczos "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  MdagMLinearOperator<Level1Op,CoarseVector> IRLHermOp(LDOp);
  Chebyshev<CoarseVector> IRLCheby(0.002,12.,151);
  FunctionHermOp<CoarseVector> IRLOpCheby(IRLCheby,IRLHermOp);
  PlainHermOp<CoarseVector> IRLOp    (IRLHermOp);
  int Nk=48;
  int Nm=64;
  int Nstop=48;
  int Nconv;
  ImplicitlyRestartedLanczos<CoarseVector> IRL(IRLOpCheby,IRLOp,Nstop,Nk,Nm,1.0e-3,20);

  std::vector<RealD>          eval(Nm);
  std::vector<CoarseVector>   evec(Nm,Coarse5d);
  CoarseVector c_src(Coarse5d);
  gaussian(CRNG,c_src);
  IRL.calc(eval,evec,c_src,Nconv);

  //  ConjugateGradient<CoarseVector>  CoarseCG(0.01,1000);
  
  ConjugateGradient<CoarseVector>  CoarseCG(0.01,2000);// 14.7s
  eval.resize(0);
  evec.resize(0,Coarse5d);
  DeflatedGuesser<CoarseVector> DeflCoarseGuesser(evec,eval);
  NormalEquations<CoarseVector> DeflCoarseCGNE(LDOp,CoarseCG,DeflCoarseGuesser);

  c_src=1.0;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Fine        PowerMethod           "<< std::endl;
  PowerMethod<LatticeFermion>       PM;   PM(HermDefOp,src);
  std::cout<<GridLogMessage << " Coarse       PowerMethod           "<< std::endl;
  MdagMLinearOperator<CoarseOperator,CoarseVector> PosdefLdop(LDOp);
  PowerMethod<CoarseVector>        cPM;  cPM(PosdefLdop,c_src);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building 2 level Multigrid            "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  typedef MultiGridPreconditioner<vSpinColourVector,  vTComplex,nbasis, DomainWallFermionR,DeflatedGuesser<CoarseVector> , NormalEquations<CoarseVector> >   TwoLevelMG;

  // MultiGrid preconditioner acting on the coarse space <-> coarsecoarse space
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(0.5,60.0,14,HermIndefOp,Ddwf); // 72 iter 63s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(0.1,60.0,20,HermIndefOp,Ddwf); // 66 iter 69s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(0.5,60.0,20,HermIndefOp,Ddwf); // 63 iter 65  s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(1.0,60.0,20,HermIndefOp,Ddwf); // 69, 70
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(1.0,60.0,14,HermIndefOp,Ddwf); // 77

  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(0.5,60.0,10,HermIndefOp,Ddwf); // 23 iter 15.9s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(0.5,60.0,14,HermIndefOp,Ddwf); // 20, 16.9s
  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(0.5,60.0,12,HermIndefOp,Ddwf); // 21, 15.6s

  //  MirsSmoother<LatticeFermion,DomainWallFermionR> FineCGSmoother(0.05,0.01,20,HermIndefOp,Ddwf);
  //  RedBlackSmoother<LatticeFermion,DomainWallFermionR> FineRBSmoother(0.00,0.001,100,Ddwf);

  // Wrap the 2nd level solver in a MultiGrid preconditioner acting on the fine space
  //  ZeroGuesser<CoarseVector> CoarseZeroGuesser;
  TwoLevelMG TwoLevelPrecon(Aggregates, LDOp,
			    HermIndefOp,Ddwf,
			    FineSmoother,
			    DeflCoarseGuesser,
			    DeflCoarseCGNE);
  TwoLevelPrecon.Level(1);

  // Apply the fine-coarse-coarsecoarse 2 deep MG preconditioner in an outer PGCR on the fine fgrid
  PrecGeneralisedConjugateResidual<LatticeFermion> l1PGCR(1.0e-8,1000,HermIndefOp,TwoLevelPrecon,16,16);
  l1PGCR.Level(1);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling 2 level Multigrid            "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  result=Zero();
  l1PGCR(src,result);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Fine CG prec DiagMooee "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;


  ConjugateGradient<LatticeFermion>          FineCG(1.0e-8,10000);
  SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> FineDiagMooee(Ddwf); //  M_ee - Meo Moo^-1 Moe 
  LatticeFermion f_src_e(FrbGrid); f_src_e=1.0;
  LatticeFermion f_res_e(FrbGrid); f_res_e=Zero();
  FineCG(FineDiagMooee,f_src_e,f_res_e);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Grid_finalize();
}
