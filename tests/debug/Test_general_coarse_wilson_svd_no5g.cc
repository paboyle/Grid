/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_padded_cell.cc

    Copyright (C) 2023

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>

using namespace std;
using namespace Grid;

template<class Fobj,class CComplex,int nbasis>
class MGPreconditioner : public LinearFunction< Lattice<Fobj> > {
public:
  using LinearFunction<Lattice<Fobj> >::operator();

  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::FineField    FineField;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseMatrix CoarseMatrix;
  typedef LinearOperatorBase<FineField>                            FineOperator;
  typedef LinearFunction    <FineField>                            FineSmoother;
  typedef LinearOperatorBase<CoarseVector>                         CoarseOperator;
  typedef LinearFunction    <CoarseVector>                         CoarseSolver;
  Aggregates     & _Aggregates;
  FineOperator   & _FineOperator;
  FineSmoother   & _PreSmoother;
  FineSmoother   & _PostSmoother;
  CoarseOperator & _CoarseOperator;
  CoarseSolver   & _CoarseSolve;

  int    level;  void Level(int lv) {level = lv; };

  MGPreconditioner(Aggregates &Agg,
		   FineOperator &Fine,
		   FineSmoother &PreSmoother,
		   FineSmoother &PostSmoother,
		   CoarseOperator &CoarseOperator_,
		   CoarseSolver &CoarseSolve_)
    : _Aggregates(Agg),
      _FineOperator(Fine),
      _PreSmoother(PreSmoother),
      _PostSmoother(PostSmoother),
      _CoarseOperator(CoarseOperator_),
      _CoarseSolve(CoarseSolve_),
      level(1)  {  }

  virtual void operator()(const FineField &in, FineField & out) 
  {
    GridBase *CoarseGrid = _Aggregates.CoarseGrid;
    //    auto CoarseGrid = _CoarseOperator.Grid();
    CoarseVector Csrc(CoarseGrid);
    CoarseVector Csol(CoarseGrid);
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    std::cout<<GridLogMessage << "Calling PreSmoother " <<std::endl;

    //    std::cout<<GridLogMessage << "Calling PreSmoother input residual "<<norm2(in) <<std::endl;
    double t;
    // Fine Smoother
    //    out = in;
    out = Zero();
    t=-usecond();
    _PreSmoother(in,out);
    t+=usecond();

    std::cout<<GridLogMessage << "PreSmoother took "<< t/1000.0<< "ms" <<std::endl;

    // Update the residual
    _FineOperator.Op(out,vec1);  sub(vec1, in ,vec1);   
    //    std::cout<<GridLogMessage <<"Residual-1 now " <<norm2(vec1)<<std::endl;

    // Fine to Coarse 
    t=-usecond();
    _Aggregates.ProjectToSubspace  (Csrc,vec1);
    t+=usecond();
    std::cout<<GridLogMessage << "Project to coarse took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse correction
    t=-usecond();
    Csol = Zero();
    _CoarseSolve(Csrc,Csol);
    //Csol=Zero();
    t+=usecond();
    std::cout<<GridLogMessage << "Coarse solve took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse to Fine
    t=-usecond();  
    //    _CoarseOperator.PromoteFromSubspace(_Aggregates,Csol,vec1);
    _Aggregates.PromoteFromSubspace(Csol,vec1); 
    add(out,out,vec1);
    t+=usecond();
    std::cout<<GridLogMessage << "Promote to this level took "<< t/1000.0<< "ms" <<std::endl;

    // Residual
    _FineOperator.Op(out,vec1);  sub(vec1 ,in , vec1);  
    //    std::cout<<GridLogMessage <<"Residual-2 now " <<norm2(vec1)<<std::endl;

    // Fine Smoother
    t=-usecond();
    //    vec2=vec1;
    vec2=Zero();
    _PostSmoother(vec1,vec2);
    t+=usecond();
    std::cout<<GridLogMessage << "PostSmoother took "<< t/1000.0<< "ms" <<std::endl;

    add( out,out,vec2);
    std::cout<<GridLogMessage << "Done " <<std::endl;
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=16;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = UGrid;
  GridRedBlackCartesian * FrbGrid = UrbGrid;

  // Construct a coarsened grid
  Coordinate clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/2;
    //    clatt[d] = clatt[d]/4;
  }
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> cseeds({5,6,7,8});
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse4d);CRNG.SeedFixedIntegers(cseeds);

  LatticeFermion    src(FGrid); random(RNG4,src);
  LatticeFermion result(FGrid); result=Zero();
  LatticeFermion    ref(FGrid); ref=Zero();
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid);

  FieldMetaData header;
  std::string file("ckpoint_lat");
  NerscIO::readConfiguration(Umu,header,file);
  
  RealD csw =0.0;
  RealD mass=-0.92;

  WilsonCloverFermionD Dw(Umu,*UGrid,*UrbGrid,mass,csw,csw);

  const int nbasis = 40;
  const int cb = 0 ;
  LatticeFermion prom(FGrid);

  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NearestStencilGeometry4D geom(Coarse4d);

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
 
  // Warning: This routine calls Linop.Op, not LinOpo.HermOp
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse4d,FGrid,cb);

  MdagMLinearOperator<WilsonCloverFermionD,LatticeFermion> MdagMOpDw(Dw);
  NonHermitianLinearOperator<WilsonCloverFermionD,LatticeFermion> LinOpDw(Dw);
  ShiftedNonHermitianLinearOperator<WilsonCloverFermionD,LatticeFermion> ShiftedLinOpDw(Dw,0.5);

  //  Aggregates.CreateSubspaceGCR(RNG4,
  //			       LinOpDw,
  //			       nbasis);
  Aggregates.CreateSubspace(RNG4,MdagMOpDw,nbasis);

  LittleDiracOperator LittleDiracOp(geom,FGrid,Coarse4d);
  LittleDiracOp.CoarsenOperator(LinOpDw,Aggregates);

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"Testing coarsened operator "<<std::endl;
  
  CoarseVector c_src (Coarse4d);
  CoarseVector c_res (Coarse4d);
  CoarseVector c_proj(Coarse4d);

  std::vector<LatticeFermion> subspace(nbasis,FGrid);
  subspace=Aggregates.subspace;

  Complex one(1.0);
  c_src = one;  // 1 in every element for vector 1.
  blockPromote(c_src,err,subspace);

  prom=Zero();
  for(int b=0;b<nbasis;b++){
    prom=prom+subspace[b];
  }
  err=err-prom; 
  std::cout<<GridLogMessage<<"Promoted back from subspace: err "<<norm2(err)<<std::endl;
  std::cout<<GridLogMessage<<"c_src "<<norm2(c_src)<<std::endl;
  std::cout<<GridLogMessage<<"prom  "<<norm2(prom)<<std::endl;

  LinOpDw.Op(prom,tmp);
  blockProject(c_proj,tmp,subspace);
  std::cout<<GridLogMessage<<" Called Big Dirac Op "<<norm2(tmp)<<std::endl;

  LittleDiracOp.M(c_src,c_res);
  std::cout<<GridLogMessage<<" Called Little Dirac Op c_src "<< norm2(c_src) << "  c_res "<< norm2(c_res) <<std::endl;

  std::cout<<GridLogMessage<<"Little dop : "<<norm2(c_res)<<std::endl;
  //  std::cout<<GridLogMessage<<" Little "<< c_res<<std::endl;

  std::cout<<GridLogMessage<<"Big dop in subspace : "<<norm2(c_proj)<<std::endl;
  //  std::cout<<GridLogMessage<<" Big "<< c_proj<<std::endl;
  c_proj = c_proj - c_res;
  std::cout<<GridLogMessage<<" ldop error: "<<norm2(c_proj)<<std::endl;
  //  std::cout<<GridLogMessage<<" error "<< c_proj<<std::endl;


  /**********
   * Some solvers
   **********
   */

  ///////////////////////////////////////
  // Coarse grid solver test
  ///////////////////////////////////////

  std::cout<<GridLogMessage<<"******************* "<<std::endl;
  std::cout<<GridLogMessage<<" Coarse Grid Solve -- Level 3 "<<std::endl;
  std::cout<<GridLogMessage<<"******************* "<<std::endl;
  TrivialPrecon<CoarseVector> simple;
  NonHermitianLinearOperator<LittleDiracOperator,CoarseVector> LinOpCoarse(LittleDiracOp);
  //  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(1.0e-4, 100, LinOpCoarse,simple,10,10); 
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(1.0e-2, 100, LinOpCoarse,simple,30,30); 
  L2PGCR.Level(3);
  c_res=Zero();
  L2PGCR(c_src,c_res);

  ////////////////////////////////////////
  // Fine grid smoother
  ////////////////////////////////////////
  std::cout<<GridLogMessage<<"******************* "<<std::endl;
  std::cout<<GridLogMessage<<" Fine Grid Smoother -- Level 2 "<<std::endl;
  std::cout<<GridLogMessage<<"******************* "<<std::endl;
  TrivialPrecon<LatticeFermionD> simple_fine;

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> SmootherGCR(0.01,1,ShiftedLinOpDw,simple_fine,6,6);
  SmootherGCR.Level(2);
  
  LatticeFermionD f_src(FGrid);
  LatticeFermionD f_res(FGrid);

  f_src = one;  // 1 in every element for vector 1.
  f_res=Zero();
  SmootherGCR(f_src,f_res);

  typedef MGPreconditioner<vSpinColourVector,  vTComplex,nbasis> TwoLevelMG;

  TwoLevelMG TwoLevelPrecon(Aggregates,
			    LinOpDw,
			    simple_fine,
			    SmootherGCR,
			    LinOpCoarse,
			    L2PGCR);
  
  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermion> L1PGCR(1.0e-8,1000,LinOpDw,TwoLevelPrecon,32,32);
  L1PGCR.Level(1);

  f_res=Zero();
  L1PGCR(f_src,f_res);

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;

  Grid_finalize();
  return 0;
}
