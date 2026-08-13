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
RealD FineSmootherShift = 0.1;
int   FineSmootherOrder = 8;
int   FineSmootherTol   = 0;
//RealD CoarseSmootherShift = 0.1;
//int   CoarseSmootherOrder = 8;
//int   CoarseSmootherTol   = 0;
RealD CoarseSolverShift   = 0.002;
RealD CoarseSolverTol     = 0.03;
int   CoarseSolverOrder   = 200;
int   CoarseMmax          = 20;   // coarse GCR restart length (was hardcoded 20)
RealD mass=0.00078;
void ParseEnvironment(void)
{
    
  if(getenv("MASS") ) mass = atof(getenv("MASS"));
  if(getenv("FineSmootherShift"))   FineSmootherShift   = atof(getenv("FineSmootherShift"));
  if(getenv("FineSmootherOrder"))   FineSmootherOrder   = atoi(getenv("FineSmootherOrder"));
  if(getenv("CoarseSolverShift"))   CoarseSolverShift   = atof(getenv("CoarseSolverShift"));
  if(getenv("CoarseSolverTol"))     CoarseSolverTol     = atof(getenv("CoarseSolverTol"));
  if(getenv("CoarseSolverOrder"))   CoarseSolverOrder   = atoi(getenv("CoarseSolverOrder"));
  if(getenv("CoarseMmax"))          CoarseMmax          = atoi(getenv("CoarseMmax"));
  if(getenv("DiagInvPrec"))
  {
    std::cout << GridLogMessage << "WARNING: DiagInvPrec option REMOVED (diagonal-inverse preconditioning wrecks fine->coarse null-vector inheritance); IGNORED" << std::endl;
  }

  //  if(getenv("CoarseSmootherShift")) CoarseSmootherShift = atof(getenv("CoarseSmootherShift"));
  //  if(getenv("CoarseSmootherOrder")) CoarseSmootherOrder = atoi(getenv("CoarseSmootherOrder"));
  
  std::cout << GridLogMessage << "PARAM: FineSmootherShift   "<<FineSmootherShift<<std::endl;
  std::cout << GridLogMessage << "PARAM: FineSmootherOrder   "<<FineSmootherOrder<<std::endl;
  //  std::cout << GridLogMessage << "PARAM: CoarseSmootherShift "<<CoarseSmootherShift<<std::endl;
  //  std::cout << GridLogMessage << "PARAM: CoarseSmootherOrder "<<CoarseSmootherOrder<<std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverShift   "<<CoarseSolverShift<<std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverTol     "<<CoarseSolverTol<<std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverOrder   "<<CoarseSolverOrder<<std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseMmax          "<<CoarseMmax<<std::endl;

  std::cout << GridLogMessage << "PARAM: MASS "<<mass<<std::endl;
}

template <class T> void readFile(T& out, std::string const fname){  
  #ifdef HAVE_LIME
    // Ref: https://github.com/paboyle/Grid/blob/feature/scidac-wp1/tests/debug/Test_general_coarse_hdcg_phys48.cc#L111
    std::cout << Grid::GridLogMessage << "Reads at: " << fname << std::endl;
    Grid::emptyUserRecord record;
    // Grid::ScidacReader SR(out.Grid()->IsBoss());
    Grid::ScidacReader SR;
    SR.open(fname);
    SR.readScidacFieldRecord(out, record);
    SR.close();
  #endif
}
template <class Field>
void saveSubspace(std::vector<Field> &subspace, std::string const fname){
  #ifdef HAVE_LIME
    std::cout << Grid::GridLogMessage << "Saving subspace (" << subspace.size() << " vectors) to: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacWriter SW(subspace[0].Grid()->IsBoss());
    SW.open(fname);
    for (int k = 0; k < (int)subspace.size(); k++)
      SW.writeScidacFieldRecord(subspace[k], record);
    SW.close();
  #endif
}

template <class Field>
void loadSubspace(std::vector<Field> &subspace, std::string const fname){
  #ifdef HAVE_LIME
    std::cout << Grid::GridLogMessage << "Loading subspace (" << subspace.size() << " vectors) from: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacReader SR;
    SR.open(fname);
    for (int k = 0; k < (int)subspace.size(); k++)
      SR.readScidacFieldRecord(subspace[k], record);
    SR.close();
  #endif
}

template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
  int nApp;
  int nAppDag;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV), nApp(0), nAppDag(0) {};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    //    std::cout << GridLogMessage<< "Op: PVdag M "<<std::endl;
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    nApp++;
  }
  void AdjOp     (const Field &in, Field &out){
    //    std::cout << GridLogMessage<<"AdjOp: Mdag PV "<<std::endl;
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
    nAppDag++;
  }

  void clear() {
    nApp = 0;
    nAppDag = 0;
  }

  void getApplications() {
    std::cout << GridLogMessage << "# applications of PVdagM: " << nApp << std::endl;
    std::cout << GridLogMessage << "# applications of PVdagM^dag: " << nAppDag << std::endl;
    std::cout << GridLogMessage << "# applications total: " << nApp + nAppDag << std::endl;
  }

  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    //    std::cout <<GridLogMessage<< "HermOp: Mdag PV PVdag M"<<std::endl;
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
    //    std::cout << "HermOp done "<<norm2(out)<<std::endl;
  }
};
template<class Matrix,class Field>
class MdagPVLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  MdagPVLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV){};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    //    std::cout <<GridLogMessage<< "Op: PVdag M "<<std::endl;
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
  }
  void AdjOp     (const Field &in, Field &out){
    //    std::cout <<GridLogMessage<< "AdjOp: Mdag PV "<<std::endl;
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    //    std::cout << GridLogMessage<<"HermOp: PVdag M Mdag PV "<<std::endl;
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
    //    std::cout << "HermOp done "<<norm2(out)<<std::endl;
  }
};
template<class Matrix,class Field>
class ShiftedPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
  int nApp;
  int nAppDag;
public:
  RealD shift;
  ShiftedPVdagMLinearOperator(RealD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV) , nApp(0), nAppDag(0){};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    //    std::cout << "Op: PVdag M "<<std::endl;
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    nApp++;
    out = out + shift * in;
  }
  void AdjOp     (const Field &in, Field &out){
    //    std::cout << "AdjOp: Mdag PV "<<std::endl;
    Field tmp(in.Grid());
    _PV.M(tmp,out);
    _Mat.Mdag(in,tmp);
    nAppDag++;
    out = out + shift * in;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
  void HermOp(const Field &in, Field &out){
    //    std::cout << "HermOp: Mdag PV PVdag M"<<std::endl;
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
  void clear() {
    nApp = 0;
    nAppDag = 0;
  }
  void getApplications() {
    std::cout << GridLogMessage << "# applications of ShiftedPVdagM: " << nApp << std::endl;
    std::cout << GridLogMessage << "# applications of ShiftedPVdagM^dag: " << nAppDag << std::endl;
    std::cout << GridLogMessage << "# applications total: " << nApp + nAppDag << std::endl;
  }
};
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
  std::string name;

  int    level;  void Level(int lv) {level = lv; };

  MGPreconditioner(Aggregates &Agg,
		   FineOperator &Fine,
		   FineSmoother &PreSmoother,
		   FineSmoother &PostSmoother,
		   CoarseOperator &CoarseOperator_,
		   CoarseSolver &CoarseSolve_,
		   std::string _name = std::string("unnamed"))
    : _Aggregates(Agg),
      _FineOperator(Fine),
      _PreSmoother(PreSmoother),
      _PostSmoother(PostSmoother),
      _CoarseOperator(CoarseOperator_),
      _CoarseSolve(CoarseSolve_),
      name(_name),
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

template<class PVdagM_t, class ShiftedPVdagM_t, class Subspace, class LittleDiracOperator, class CoarseVector, class TwoLevelMG>
void runMG(
  GridCartesian *FGrid,
  GridCartesian *Coarse5d, 
  NextToNearestStencilGeometry5D geom,
  PVdagM_t PVdagM,
  ShiftedPVdagM_t ShiftedPVdagM,
  // std::vector<LatticeFermion> subspace
  Subspace AggregatesPD
) {

  // typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  // typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  // typedef LittleDiracOperator::CoarseVector CoarseVector;
  ParseEnvironment();

  std::vector<LatticeFermion> subspace = AggregatesPD.subspace;
  int nbasis = subspace.size();
  const int cb = 0 ;

  LatticeFermion err(FGrid);
  LatticeFermion prom(FGrid);
  LatticeFermion tmp(FGrid);

  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  CoarseVector c_proj(Coarse5d);
  Complex one(1.0);

  LatticeFermionD f_src(FGrid);
  LatticeFermionD f_res(FGrid);
  // typedef MGPreconditioner<vSpinColourVector,  vTComplex,nbasis> TwoLevelMG;
  TrivialPrecon<CoarseVector> simple;
  TrivialPrecon<LatticeFermionD> simple_fine;

  // Subspace AggregatesPD(Coarse5d,FGrid,cb);

  // Orthonormalize subspace and compute nulliness

  ShiftedPVdagM.shift = CoarseSolverShift;
  int nonherm = 0;
  LittleDiracOperator LittleDiracOpPV(geom,FGrid,Coarse5d,nonherm);
  LittleDiracOpPV.CoarsenOperator(ShiftedPVdagM, AggregatesPD);
  ShiftedPVdagM.shift = FineSmootherShift;

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"Testing coarsened operator "<<std::endl;

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

  // PVdagM.Op(prom,tmp);
  // blockProject(c_proj,tmp,subspace);
  // std::cout<<GridLogMessage<<" Called Big Dirac Op "<<norm2(tmp)<<std::endl;

  // LittleDiracOpPV.M(c_src,c_res);
  // std::cout<<GridLogMessage<<" Called Little Dirac Op c_src "<< norm2(c_src) << "  c_res "<< norm2(c_res) <<std::endl;

  // std::cout<<GridLogMessage<<"Little dop : "<<norm2(c_res)<<std::endl;
  // //  std::cout<<GridLogMessage<<" Little "<< c_res<<std::endl;
  // std::cout<<GridLogMessage<<"Big dop in subspace : "<<norm2(c_proj)<<std::endl;
  // //  std::cout<<GridLogMessage<<" Big "<< c_proj<<std::endl;
  // c_proj = c_proj - c_res;
  // std::cout<<GridLogMessage<<" ldop error: "<<norm2(c_proj)<<std::endl;
  // //  std::cout<<GridLogMessage<<" error "<< c_proj<<std::endl;

  ///////////////////////////////////////
  // Coarse grid solver test
  ///////////////////////////////////////

  std::cout<<GridLogMessage<<"******************* "<<std::endl;
  std::cout<<GridLogMessage<<" Coarse Grid Solve -- Level 2 "<<std::endl;
  std::cout<<GridLogMessage<<"******************* "<<std::endl;
  NonHermitianLinearOperator<LittleDiracOperator,CoarseVector> LinOpCoarse(LittleDiracOpPV);
  // DiagonalInverse preconditioning REMOVED (library support withdrawn: it
  // wrecks the collinearity that makes fine->coarse null-vector inheritance
  // free).  TrivialPrecon reproduces the former DiagInvPrec=0 path exactly.
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(CoarseSolverTol, (CoarseSolverOrder+CoarseMmax-1)/CoarseMmax, LinOpCoarse,simple,CoarseMmax,CoarseMmax);
  L2PGCR.SetZeroGuess(1);              // callers zero Csol / c_res
  L2PGCR.Level(2);
  L2PGCR.Name("Couter");
  c_res=Zero();
  L2PGCR(c_src,c_res);


  ////////////////////////////////////////
  // Fine grid smoother
  ////////////////////////////////////////
  //  NonHermitianLinearOperator<PVdagM_t,LatticeFermionD> LinOpSmooth(PVdagM);
  
  //  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> SmootherGCR(0.05,1,ShiftedPVdagM,simple_fine,8,8);
  // Force 10 iters exactly, no early termination
  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> SmootherGCR(FineSmootherTol,1,
									    ShiftedPVdagM,simple_fine,
									    FineSmootherOrder,FineSmootherOrder);
  SmootherGCR.Level(1);
  SmootherGCR.Name("Fsmoother");
  SmootherGCR.SetZeroGuess(1);         // pre/post slots + direct call all zero their guess

  f_src = one;  // 1 in every element for vector 1.
  f_res=Zero();
  SmootherGCR(f_src,f_res);

  TwoLevelMG TwoLevelPrecon(AggregatesPD,
			    PVdagM,
			    simple_fine,
			    SmootherGCR,
			    LinOpCoarse,
			    L2PGCR,
			    "PVdagM");
  
  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermion> L1PGCR(1.0e-8,1000,PVdagM,TwoLevelPrecon,32,32);
  L1PGCR.SetZeroGuess(1);              // f_res=Zero() before the solve
  L1PGCR.Level(1);
  L1PGCR.Name("Fouter");

  std::cout<<GridLogMessage<<"******************* "<<std::endl;
  std::cout<<GridLogMessage<<" Running Multi Grid Solver "<<std::endl;
  std::cout<<GridLogMessage<<"******************* "<<std::endl;
  f_res=Zero();
  L1PGCR(f_src,f_res);

  std::cout << GridLogMessage << "Fine Grid Smoother -- Level 2 operator uses: " << std::endl;
  PVdagM.getApplications();
  PVdagM.clear();
  ShiftedPVdagM.getApplications();
  ShiftedPVdagM.clear();

}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  // TODO read in more parameters: nbasis, GCR iters, smoother order, m
  // Might be impossible because nbasis needs to be a constant to be a template parameter
  // Usage : $ ./Example_pvdagm <nbasis> <smooth> <outerIters> <m>
  // std::string nbasisStr  = argv[1];
  // std::string smoothStr  = argv[2];
  // std::string outerStr   = argv[3];
  // std::string mStr       = argv[4];
  // int nbasis = std::stoi(nbasisStr);
  // int smooth = std::stoi(smoothStr);

  const int Ls=24;
  RealD M5=1.8;

  
  //  const int nbasis = 40;
  const int nbasis = 60;
  
  std::cout << GridLogMessage << "Mass: " << mass << ", Ls: " << Ls << ", running Mobius kernel with b=1.5, c=0.5" << std::endl;
  std::cout << GridLogMessage << "nbasis: " << nbasis << std::endl;

  std::vector<int> lat_size {48, 48, 48, 96};

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Construct a coarsened grid
  // Coordinate clatt = GridDefaultLatt();
  Coordinate clatt = lat_size;
  Coordinate Block({4,4,4,4});
  std::cout << GridLogMessage << "Lattice size: " << lat_size << std::endl;
  for(int d=0;d<clatt.size();d++){
    clatt[d] = lat_size[d]/Block[d];
  }

  std::cout << GridLogMessage << "constructing coarse grid" << std::endl;
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  LatticeFermion    src(FGrid); random(RNG5,src);
  LatticeFermion result(FGrid); result=Zero();
  LatticeFermion    ref(FGrid); ref=Zero();
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid);


  std::cout << GridLogMessage << "Reading in gauge field" << std::endl;
  FieldMetaData header;
  // std::string file("/sdcc/u/poare/PETSc-Grid/ckpoint_lat.4000");
  std::string file("/ccs/home/poare/ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  /*

  // DWF, m=0.01
  // std::string eigenPath = "/hpcgpfs01/work/lqcd/staging/RBC/ckpoint_lat.4000/ks_evecs/PVdagM_Nm80_Nk40_Niter5000_337342/";

  // DWF, m=0.001
  // std::string eigenPath = "/hpcgpfs01/work/lqcd/staging/RBC/ckpoint_lat.4000/ks_evecs/PVdagM_Nm80_Nk40_Niter5000_m0p001_339143/";

  // Mobius, m=0.001
  // std::string eigenPath = "/hpcgpfs01/work/lqcd/staging/RBC/ckpoint_lat.4000/ks_evecs/PVdagM_Nm80_Nk40_Niter5000_346851/";

  // Frontier path
  std::string eigenPath = "/ccs/home/poare/lqcd/multigrid/spectra/ckpoint_lat.1000/...";
  

  std::cout << GridLogMessage << "Loading eigenvalues" << std::endl;
  std::ifstream evalFile(eigenPath + "evals.txt");
  std::string str;
  std::vector<ComplexD> evals;
  while (std::getline(evalFile, str)) {
      std::cout << GridLogMessage << "Reading line: " << str << std::endl;
      int i1 = str.find("(") + 1;
      int i2 = str.find(",") + 1;
      int i3 = str.find(")");
      std::cout << "i1,i2,i3 = " << i1 << "," << i2 << "," << i3 << std::endl;
      std::string reStr = str.substr(i1, i2 - i1);
      std::string imStr = str.substr(i2, i3 - i2);
      std::cout << GridLogMessage << "Parsed re = " << reStr << " and im = " << imStr << std::endl;
      // ComplexD z (std::stof(reStr), std::stof(imStr));
      ComplexD z (std::stod(reStr), std::stod(imStr));
      evals.push_back(z);
  }
  std::cout << GridLogMessage << "Eigenvalues: " << evals << std::endl;

  int Nevecs = 20;
  std::vector<LatticeFermion> evecs;
  LatticeFermion evec (FGrid);
  for (int i = 0; i < Nevecs; i++) {
    std::string evecPath = eigenPath + "evec" + std::to_string(i);
    readFile(evec, evecPath);
    evecs.push_back(evec);
  }
  std::cout << GridLogMessage << "Evecs loaded" << std::endl;

  */
  // TODO uncomment when evecs are computed!

  // DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  // DomainWallFermionD Dpv(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,M5);

  // Mobius
  RealD b=1.5;// Scale factor b+c=2, b-c=1
  RealD c=0.5;
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  MobiusFermionD Dpv(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,M5,b,c);

  const int cb = 0 ;
  LatticeFermion prom(FGrid);

  // assert(nbasis <= Nevecs);     // need to have enough evecs

  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNearestStencilGeometry5D geom(Coarse5d);

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;

  // typedef PVdagMLinearOperator<DomainWallFermionD,LatticeFermionD> PVdagM_t;
  // typedef MdagPVLinearOperator<DomainWallFermionD,LatticeFermionD> MdagPV_t;
  // typedef ShiftedPVdagMLinearOperator<DomainWallFermionD,LatticeFermionD> ShiftedPVdagM_t;
  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD> PVdagM_t;
  typedef MdagPVLinearOperator<MobiusFermionD,LatticeFermionD> MdagPV_t;
  typedef ShiftedPVdagMLinearOperator<MobiusFermionD,LatticeFermionD> ShiftedPVdagM_t;

  PVdagM_t PVdagM(Ddwf,Dpv);
  MdagPV_t MdagPV(Ddwf,Dpv);
  //  ShiftedPVdagM_t ShiftedPVdagM(2.0,Ddwf,Dpv); // 355
  //  ShiftedPVdagM_t ShiftedPVdagM(1.0,Ddwf,Dpv); // 246
  //  ShiftedPVdagM_t ShiftedPVdagM(0.5,Ddwf,Dpv); // 183
  //  ShiftedPVdagM_t ShiftedPVdagM(0.25,Ddwf,Dpv); // 145
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 134
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 127 -- NULL space via inverse iteration
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 57 -- NULL space via inverse iteration; 3 iterations
  //  ShiftedPVdagM_t ShiftedPVdagM(0.25,Ddwf,Dpv); // 57 , tighter inversion
  //  ShiftedPVdagM_t ShiftedPVdagM(0.25,Ddwf,Dpv); // nbasis 20 -- 49 iters
  //  ShiftedPVdagM_t ShiftedPVdagM(0.25,Ddwf,Dpv); // nbasis 20 -- 70 iters; asymmetric 
  //  ShiftedPVdagM_t ShiftedPVdagM(0.25,Ddwf,Dpv); // 58; Loosen coarse, tighten fine
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 56 ... 
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 51 ...  with 24 vecs
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 31 ...  with 24 vecs and 2^4 blocking
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 43 ...  with 16 vecs and 2^4 blocking, sloppier
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 35  ...  with 20 vecs and 2^4 blocking
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 35  ...  with 20 vecs and 2^4 blocking, looser coarse
  //  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv); // 64  ...  with 20 vecs, Christoph setup, and 2^4 blocking, looser coarse

  ShiftedPVdagM_t ShiftedPVdagM(FineSmootherShift,Ddwf,Dpv); // 

  // Run power method on HOA??
  PowerMethod<LatticeFermion>       PM;

  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  CoarseVector c_proj(Coarse5d);
  Complex one(1.0);

  std::vector<LatticeFermion> subspace(nbasis,FGrid);

  LatticeFermionD f_src(FGrid);
  LatticeFermionD f_res(FGrid);
  typedef MGPreconditioner<vSpinColourVector,  vTComplex,nbasis> TwoLevelMG;
  TrivialPrecon<CoarseVector> simple;
  TrivialPrecon<LatticeFermionD> simple_fine;

  // Warning: This routine calls PVdagM.Op, not PVdagM.HermOp
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;

  // Breeds right singular vectors with call to HermOp (V)
  // int chebyOrd = 500;
  // V.CreateSubspaceChebyshev(RNG5,PVdagM,
	// 		    nbasis,
	// 		    4000.0,0.003,
	// 		    chebyOrd);
  // AggregatesPD.CreateSubspaceChebyshev(RNG5,
	// 			       PVdagM,
	// 			       nbasis,
	// 			       4000.0,
	// 			       0.003,
	// 			       chebyOrd);
  
  // Subspace testing (uncomment blocks when needed)

  // - nbasis = 20, m=0.01, 35 outer iterations
  // - nbasis = 40, m=0.01, 23 outer iterations
  std::cout << GridLogMessage << "*** GCR setup ***" << std::endl;

  // Subspace cache: save after generation, reload on subsequent runs to skip expensive setup.
  // Set SUBSPACE_FILE to override the default path.
  std::string subspace_file = "/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/subspace_nb"
                              + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));

  // Check if subspace file exists (boss rank checks, result broadcast via GlobalSum).
  uint64_t file_exists = 0;
  if ( UGrid->IsBoss() ) {
    std::ifstream f(subspace_file);
    file_exists = f.good() ? 1 : 0;
  }
  UGrid->GlobalSum(file_exists);

  Subspace AggregatesGCR(Coarse5d,FGrid,cb);

  if ( file_exists ) {
    std::cout << GridLogMessage << "*** Loading subspace from disk ***" << std::endl;
    loadSubspace(AggregatesGCR.subspace, subspace_file);
    // Insurance: GLOBAL (whole-lattice) orthonormalise, matching what
    // CreateSubspaceGCR applies to generated subspaces (Aggregates.h:196), so a
    // reloaded file ends in the same state.  This replaces the block
    // Orthogonalise() previously called here -- that is redundant (CoarsenOperator
    // block-GS's the subspace internally) and would leave a loaded file block-
    // orthonormal while a generated one is globally orthonormal.  Global GS is
    // span-preserving, so the coarse operator is unchanged.
    AggregatesGCR.GlobalOrthonormalise();
    std::cout << GridLogMessage << "Subspace loaded and globally orthonormalised." << std::endl;
  } else {
    std::cout << GridLogMessage << "*** GCR subspace generation ***" << std::endl;
    AggregatesGCR.CreateSubspaceGCR(RNG5,PVdagM,nbasis);
    std::cout << GridLogMessage << "Subspace generation: PVdagM operator uses:" << std::endl;
    PVdagM.getApplications();
    PVdagM.clear();
    saveSubspace(AggregatesGCR.subspace, subspace_file);
    std::cout << GridLogMessage << "Subspace saved to: " << subspace_file << std::endl;
  }
  
  std::cout << GridLogMessage << "Basis construction operator uses: " << std::endl;
  PVdagM.getApplications();
  PVdagM.clear();

  std::cout << GridLogMessage << "Calling runMG " << std::endl;
  runMG<PVdagM_t, ShiftedPVdagM_t, Subspace, LittleDiracOperator, CoarseVector, TwoLevelMG>(
    FGrid, 
    Coarse5d, 
    geom, 
    PVdagM, 
    ShiftedPVdagM, 
    AggregatesGCR
  );

  //////////////////////////////////
  // Standard CG
  //////////////////////////////////
#if 0
  {
  std::cout << "**************************************"<<std::endl;
  std::cout << "Calling red black CG"<<std::endl;
  std::cout << "**************************************"<<std::endl;
    ConjugateGradient<LatticeFermionD>  CGfine(1.0e-8,30000,false);
    SchurDiagMooeeOperator<MobiusFermionD, LatticeFermion> HermOpEO(Ddwf);
      
    LatticeFermion result(FrbGrid); result=Zero();
    LatticeFermion    src(FrbGrid); random(RNG5,src);
    result=Zero();

    CGfine(HermOpEO, src, result);
  }
  {
    std::cout << "**************************************"<<std::endl;
    std::cout << "Calling MdagM CG"<<std::endl;
    std::cout << "**************************************"<<std::endl;
      
    LatticeFermion result(FGrid); result=Zero();
    LatticeFermion    src(FGrid); random(RNG5,src);
    result=Zero();

    MdagMLinearOperator<MobiusFermionD, LatticeFermionD> HermOp(Ddwf);
    ConjugateGradient<LatticeFermionD>  CGfine(1.0e-8,100000,false);
    CGfine(HermOp, src, result);
  }
  {
    std::cout << "**************************************"<<std::endl;
    std::cout << "Calling PVdagM GCR"<<std::endl;
    std::cout << "**************************************"<<std::endl;

    LatticeFermion result(FGrid); result=Zero();
    LatticeFermion    src(FGrid); random(RNG5,src);
    result=Zero();
    PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> GCR(1.0e-8,3000,PVdagM,simple_fine,50,50);
    GCR.Name("Fbaseline");
    GCR.SetZeroGuess(1);               // result=Zero() above
    GCR(src,result);
  }
#endif  

  
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;

  Grid_finalize();
  return 0;
}
