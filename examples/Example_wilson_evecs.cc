/*************************************************************************************

    Script for studying the Wilson eigenvectors resulting from the Krylov-Schur process. 

    Usage : 
      $ ./Example_spec_kryschur <Nm> <Nk> <maxiter> <Nstop> <inFile> <outDir> <?rf>

      Nm = Maximum size of approximation subspace.
      Nk = Size of truncation subspace
      maxiter = Maximum number of iterations.
      Nstop   = Stop when Nstop eigenvalues have converged. 
      inFile  = Gauge configuration to read in.
      outDir  = Directory to write output to.
      rf      = (Optional) RitzFilter to sort with. Takes in any string in 
                  {EvalNormSmall, EvalNormLarge, EvalReSmall, EvalReLarge, EvalImSmall, EvalImLarge}
    
    Output:
      ${outDir}/evals.txt  = Contains all eigenvalues. Each line is formatted as `$idx $eval $ritz`, where:
                              - $idx is the index of the eigenvalue.
                              - $eval is the eigenvalue, formated as "(re,im)".
                              - $ritz is the Ritz estimate of the eigenvalue (deviation from being a true eigenvalue)
      ${outDir}/evec${idx} = Eigenvector $idx written out in SCIDAC format (if LIME is enabled).

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_padded_cell.cc

    Copyright (C) 2023

    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: Patrick Oare <poare@bnl.edu>

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
  std::vector<RealD> ritz  = KS.getRitzEstimates();
  for (int i = 0; i < Nk; i++) {
    // write eigenvalues and Ritz estimates
    fEval << i << " " << evals(i) << " " << ritz[i];
    if (i < Nk - 1) { fEval << "\n"; }
  }
  fEval.close();
  
  // Write evecs
  int Nevecs = Nk;          // don't write all of them
  std::vector<Field> evecs = KS.getEvecs();
  for (int i = 0; i < Nevecs; i++) {
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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  // Usage : $ ./Example_wilson_evecs ${inFile}
  std::string file       = argv[1];

  const int Ls=16;

//   GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  //std::vector<int> lat_size {16, 16, 16, 32};
  std::vector<int> lat_size {32, 32, 32, 32};
  std::cout << "Lattice size: " << lat_size << std::endl;
  GridCartesian * UGrid = SpaceTimeGrid::makeFourDimGrid(lat_size, 
								          GridDefaultSimd(Nd,vComplex::Nsimd()),
								          GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  // GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  // GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
  GridCartesian * FGrid = UGrid;
  GridRedBlackCartesian * FrbGrid = UrbGrid;

  std::vector<int> seeds4({1,2,3,4});
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermion    src(FGrid); random(RNG4, src);
  LatticeGaugeField Umu(UGrid);

  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, file);

  std::cout << GridLogMessage << "Loaded configuration" << std::endl;

  // RealD mass = 0.01;
  RealD M5 = 1.8;

  // Wilson mass
  RealD mass = -1.6;

  std::cout << GridLogMessage << "masses specified" << std::endl;

  std::vector<Complex> boundary = {1,1,1,-1};
  WilsonFermionD::ImplParams Params(boundary);

  // DomainWallFermionD Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5);
  // NonHermitianLinearOperator<DomainWallFermionD, LatticeFermionD> DLinOp (Ddwf);

  // WilsonFermionD Dwilson(Umu, *FGrid, *FrbGrid, mass);
  WilsonFermionD Dwilson(Umu, *UGrid, *UrbGrid, mass, Params);
  NonHermitianLinearOperator<WilsonFermionD, LatticeFermionD> DLinOp (Dwilson);

  std::cout << GridLogMessage << "Dirac operator defined" << std::endl;

  std::string eigenPath = "/home/poare/lqcd/multigrid/spectra/32cube-rho0.124-tau4/U_smr_3.000000/Nm72_Nk24_8111835.aurora-pbs-0001.hostmgmt.cm.aurora.alcf.anl.gov/";

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

  int Nevecs = 24;
  std::vector<LatticeFermion> evecs;
  LatticeFermion evec (FGrid);
  for (int i = 0; i < Nevecs; i++) {
    std::string evecPath = eigenPath + "evec" + std::to_string(i);
    readFile(evec, evecPath);
    evecs.push_back(evec);
  }
  std::cout << GridLogMessage << "Evecs loaded" << std::endl;

  // Compute < evec | D - \lambda | evec >
  std::cout << GridLogMessage << "Testing eigenvectors" << std::endl;
  LatticeFermion Devec (FGrid);
  ComplexD ritz;
  for (int i = 0; i < Nevecs; i++) {
    Devec = Zero();
    DLinOp.Op(evecs[i], Devec);
    ritz = std::sqrt(norm2(Devec - evals[i] * evecs[i]));
    std::cout << GridLogMessage << "i = " << i << ", || (D - lambda) |vi> || = " << ritz << std::endl;
  }
  // Eigen::MatrixXcd Dw_evecs;
  // Dw_evecs = Eigen::MatrixXcd::Zero(Nevecs, Nevecs);
  // for (int i = 0; i < Nevecs; i++) {
  //   Linop.Op(evecs[i], Devec);
  //   for (int j = 0; j < Nevecs; j++) {

  //   }
  // }

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;

  Grid_finalize();
  return 0;
}
