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

// copied here from Test_general_coarse_pvdagm.cc

// copied here from Test_general_coarse_pvdagm.cc

#include <cstdlib>

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>

#include <Grid/parallelIO/IldgIOtypes.h>
#include <Grid/parallelIO/IldgIO.h>

using namespace std;
using namespace Grid;

template <class T> void writeFile(T& in, std::string const fname){  
  #ifdef HAVE_LIME
    // Ref: https://github.com/paboyle/Grid/blob/feature/scidac-wp1/tests/debug/Test_general_coarse_hdcg_phys48.cc#L111
    std::cout << Grid::GridLogMessage << "Writes to: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacWriter WR(in.Grid()->IsBoss());
    WR.open(fname);
    WR.writeScidacFieldRecord(in,record,0); // Lexico
    WR.close();
  #endif
}

/**
 * Writes the eigensystem of a Krylov Schur object to a directory. 
 * 
 * Parameters
 * ----------
 * std::string path
 *    Directory to write to. 
 */
template <class Field>
void writeEigensystem(KrylovSchur<Field> KS, std::string outDir) {
  int Nk = KS.getNk();
  std::cout << GridLogMessage << "Writing output to directory: " << outDir << std::endl;
  
  // Write evals
  std::string evalPath = outDir + "/evals.txt";
  std::ofstream fEval;
  fEval.open(evalPath);
  Eigen::VectorXcd evals = KS.getEvals();
  for (int i = 0; i < Nk; i++) {
    // write Eigenvalues
    fEval << i << " " << evals(i);
    if (i < Nk - 1) { fEval << "\n"; }
  }
  fEval.close();
  
  // Write evecs
  std::vector<Field> evecs = KS.getEvecs();
  for (int i = 0; i < Nk; i++) {
    std::string fName = outDir + "/evec" + std::to_string(i);
    writeFile(evecs[i], fName);     // using method from Grid/HMC/ComputeWilsonFlow.cc
  }
}

// Hermitize a DWF operator by squaring it
template<class Matrix,class Field>
class SquaredLinearOperator : public LinearOperatorBase<Field> {

  public:
  Matrix &_Mat;

  public:
    SquaredLinearOperator(Matrix &Mat): _Mat(Mat) {};

    void OpDiag (const Field &in, Field &out) {    assert(0);  }
    void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
    void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
    void Op     (const Field &in, Field &out){
      // std::cout << "Op is overloaded as HermOp" << std::endl;
      HermOp(in, out);
    }
    void AdjOp     (const Field &in, Field &out){
      HermOp(in, out);
    }
    void _Op     (const Field &in, Field &out){
      // std::cout << "Op: M "<<std::endl;
      _Mat.M(in, out);
    }
    void _AdjOp     (const Field &in, Field &out){
      // std::cout << "AdjOp: Mdag "<<std::endl;
      _Mat.Mdag(in, out);
    }
    void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
    void HermOp(const Field &in, Field &out){
      // std::cout << "HermOp: Mdag M Mdag M"<<std::endl;
      Field tmp(in.Grid());
      _Op(in,tmp);
      _AdjOp(tmp,out);
    }
};

template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV){};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    std::cout << "Op: PVdag M "<<std::endl;
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
  }
  void AdjOp     (const Field &in, Field &out){
    std::cout << "AdjOp: Mdag PV "<<std::endl;
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
  void HermOp(const Field &in, Field &out){
    std::cout << "HermOp: Mdag PV PVdag M"<<std::endl;
    Field tmp(in.Grid());
    //    _Mat.M(in,tmp);
    //    _PV.Mdag(tmp,out);
    //    _PV.M(out,tmp);
    //    _Mat.Mdag(tmp,out);
    Op(in,tmp);
    AdjOp(tmp,out);
    //    std::cout << "HermOp done "<<norm2(out)<<std::endl;
  }
};

template<class Matrix,class Field>
class ShiftedPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
  RealD shift;
public:
  ShiftedPVdagMLinearOperator(RealD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV){};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    std::cout << "Op: PVdag M "<<std::endl;
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    out = out + shift * in;
  }
  void AdjOp     (const Field &in, Field &out){
    std::cout << "AdjOp: Mdag PV "<<std::endl;
    Field tmp(in.Grid());
    _PV.M(tmp,out);
    _Mat.Mdag(in,tmp);
    out = out + shift * in;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
  void HermOp(const Field &in, Field &out){
    std::cout << "HermOp: Mdag PV PVdag M"<<std::endl;
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

template<class Matrix, class Field>
class ShiftedComplexPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
  ComplexD shift;
public:
ShiftedComplexPVdagMLinearOperator(ComplexD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV){};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    std::cout << "Op: PVdag M "<<std::endl;
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    out = out + shift * in;
  }
  void AdjOp     (const Field &in, Field &out){
    std::cout << "AdjOp: Mdag PV "<<std::endl;
    Field tmp(in.Grid());
    _PV.M(tmp,out);
    _Mat.Mdag(in,tmp);
    out = out + shift * in;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
  void HermOp(const Field &in, Field &out){
    std::cout << "HermOp: Mdag PV PVdag M"<<std::endl;
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
  
  void resetShift(ComplexD newShift) {
    shift = newShift;
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

  // Usage : $ ./Example_spec_kryschur <Nm> <Nk> <maxiter> <Nstop> <inFile> <outDir>
  // assert (argc == 5);
  std::string NmStr      = argv[1];
  std::string NkStr      = argv[2];
  std::string maxIterStr = argv[3];
  std::string NstopStr   = argv[4];
  std::string file       = argv[5];
  std::string outDir     = argv[6];

  RitzFilter RF;
  if (argc == 8) {
    std::string rf       = argv[7];
    RF = selectRitzFilter(rf);
  } else {
    RF = EvalReSmall;
  }
  std::cout << "Sorting eigenvalues using " << rfToString(RF) << std::endl;

  //const int Ls=16;
  const int Ls = 8;

//   GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  //std::vector<int> lat_size {16, 16, 16, 32};
  std::vector<int> lat_size {8, 8, 8, 8};
  std::cout << "Lattice size: " << lat_size << std::endl;
  GridCartesian * UGrid = SpaceTimeGrid::makeFourDimGrid(lat_size, 
								          GridDefaultSimd(Nd,vComplex::Nsimd()),
								          GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Construct a coarsened grid
  // poare TODO: replace this with the following line?
  Coordinate clatt = lat_size;
//   Coordinate clatt = GridDefaultLatt();              // [PO] initial line before I edited it
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/2;
    //    clatt[d] = clatt[d]/4;
  }
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
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

  FieldMetaData header;
  // std::string file ("/sdcc/u/poare/PETSc-Grid/ckpoint_EODWF_lat.125");
  // std::string file("/Users/patrickoare/libraries/PETSc-Grid/ckpoint_EODWF_lat.125");
  NerscIO::readConfiguration(Umu,header,file);

  // RealD mass=0.01;
  RealD mass=0.001;
  RealD M5=1.8;

  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionD Dpv(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,M5);

  // const int nbasis = 20;            // size of approximate basis for low-mode space
  const int nbasis = 3;            // size of approximate basis for low-mode space
  const int cb = 0 ;
  LatticeFermion prom(FGrid);

  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNearestStencilGeometry5D geom(Coarse5d);

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;

  typedef PVdagMLinearOperator<DomainWallFermionD,LatticeFermionD> PVdagM_t;
  typedef ShiftedPVdagMLinearOperator<DomainWallFermionD,LatticeFermionD> ShiftedPVdagM_t;
  typedef ShiftedComplexPVdagMLinearOperator<DomainWallFermionD,LatticeFermionD> ShiftedComplexPVdagM_t;
  PVdagM_t PVdagM(Ddwf, Dpv);
  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv);
  SquaredLinearOperator<DomainWallFermionD, LatticeFermionD> Dsq (Ddwf);
  NonHermitianLinearOperator<DomainWallFermionD, LatticeFermionD> DLinOp (Ddwf);

  // int Nm = 200;
  // int Nk = 110;
  // int maxIter = 2000;
  // int Nstop = 100;

  int Nm = std::stoi(NmStr);
  int Nk = std::stoi(NkStr);
  int maxIter = std::stoi(maxIterStr);
  int Nstop = std::stoi(NstopStr);

  std::cout << GridLogMessage << "Runnning Krylov Schur. Nm = " << Nm << ", Nk = " << Nk << ", maxIter = " << maxIter 
                  << ", Nstop = " << Nstop << std::endl;
  // Arnoldi Arn(PVdagM, FGrid, 1e-8);
  // Arn(src, maxIter, Nm, Nk, Nstop);
  KrylovSchur KrySchur (PVdagM, FGrid, 1e-8, RF);
  KrySchur(src, maxIter, Nm, Nk, Nstop);

  std::cout<<GridLogMessage << "*******************************************" << std::endl;
  std::cout<<GridLogMessage << "***************** RESULTS *****************" << std::endl;
  std::cout<<GridLogMessage << "*******************************************" << std::endl;

  std::cout << GridLogMessage << "Krylov Schur eigenvalues: " << std::endl << KrySchur.getEvals() << std::endl;
  //std::cout << GridLogMessage << "Lanczos eigenvalues: " << std::endl << levals << std::endl;

  // KrySchur.writeEigensystem(outDir);

  writeEigensystem(KrySchur, outDir);

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;

  Grid_finalize();
  return 0;
}
