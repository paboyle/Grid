/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./examples/Example_krylov_schur.cc

    Copyright (C) 2023

    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: Patrick Oare <patrickoare@gmail.com>
    Author: Chulwoo Jung <chulwoo@bnl.gov>

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

#include <cstdlib>

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>

using namespace std;
using namespace Grid;

namespace Grid {

struct LanczosParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
		  		RealD, mass , 
		  		RealD, mstep , 
				Integer, Nstop,
                                Integer, Nk,
                                Integer, Np,
                                Integer, ReadEvec,
                                Integer, maxIter,
                                Integer, Nblock,
                                Integer, verify,
		  		RealD, shift , 
	  			RealD, resid,
	  			RealD, ChebyLow,
	  			RealD, ChebyHigh,
	  			Integer, ChebyOrder)

  LanczosParameters() {
    ////////////////////////////// Default values
      mass = 0;
    /////////////////////////////////
  }

  template <class ReaderClass >
  LanczosParameters(Reader<ReaderClass> & TheReader){
    initialize(TheReader);
  }

  template < class ReaderClass > 
  void initialize(Reader<ReaderClass> &TheReader){
//    std::cout << GridLogMessage << "Reading HMC\n";
    read(TheReader, "HMC", *this);
  }


  void print_parameters() const {
//    std::cout << GridLogMessage << "[HMC parameters] Trajectories            : " << Trajectories << "\n";
//    std::cout << GridLogMessage << "[HMC parameters] Start trajectory        : " << StartTrajectory << "\n";
//    std::cout << GridLogMessage << "[HMC parameters] Metropolis test (on/off): " << std::boolalpha << MetropolisTest << "\n";
//    std::cout << GridLogMessage << "[HMC parameters] Thermalization trajs    : " << NoMetropolisUntil << "\n";
//    std::cout << GridLogMessage << "[HMC parameters] Starting type           : " << StartingType << "\n";
//    MD.print_parameters();
  }
  
};

}

template <class T> void writeFile(T& in, std::string const fname){
#if 1
  // Ref: https://github.com/paboyle/Grid/blob/feature/scidac-wp1/tests/debug/Test_general_coarse_hdcg_phys48.cc#L111
  std::cout << Grid::GridLogMessage << "Writes to: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacWriter WR(in.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(in,record,0);
  WR.close();
#endif
  // What is the appropriate way to throw error?
}


typedef WilsonFermionD WilsonOp;
typedef typename WilsonFermionD::FermionField FermionField;

template<class Matrix,class Field>
class InvertNonHermitianLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  RealD _stp;
public:
  InvertNonHermitianLinearOperator(Matrix &Mat,RealD stp=1e-8): _Mat(Mat),_stp(stp){};
  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {
    assert(0);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    assert(0);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
//    _Mat.MdirAll(in,out);
    assert(0);
  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.Mdag(in,tmp);
    MdagMLinearOperator<Matrix,Field> HermOp(_Mat);
    ConjugateGradient<Field> CG(_stp,10000);
    CG(HermOp,tmp,out);
//    out = out + shift * in;
  }
  void AdjOp     (const Field &in, Field &out){
    _Mat.Mdag(in,out);
//    out = out + shift * in;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    assert(0);
  }
  void HermOp(const Field &in, Field &out){
    assert(0);
  }
};

template<class Field>
void testSchurFromHess(Arnoldi<Field>& Arn, Field& src, int Nlarge, int Nm, int Nk) {

  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout << GridLogMessage << "Testing Schur reordering, Nm = " << Nm << ", Nk = " << Nk << std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  std::cout << GridLogMessage << "Running Arnoldi for 1 iteration to get a Hessenberg." << std::endl;
  Arn(src, 1, Nlarge, Nm, Nlarge);
  Eigen::MatrixXcd Hess = Arn.getHessenbergMat();
  std::cout << GridLogMessage << "Hessenberg for use: " << std::endl << Hess << std::endl;

  ComplexSchurDecomposition schur (Hess, true);
  bool isDecomposed = schur.checkDecomposition();
  std::cout << "Schur decomp holds? " << isDecomposed << std::endl;

  std::cout << GridLogMessage << "S = " << std::endl << schur.getMatrixS() << std::endl;
  std::cout << GridLogMessage << "Swapping S(3, 3) with S(4, 4)" << std::endl;
  schur.swapEvals(3);
  std::cout << GridLogMessage << "S after swap = " << std::endl << schur.getMatrixS() << std::endl;
  std::cout << "Schur decomp still holds? " << schur.checkDecomposition() << std::endl;

  // Now move last diagonal element all the way to the front.
  std::cout << GridLogMessage << "Moving last eval to front. S at start = " << std::endl << schur.getMatrixS() << std::endl;
  for (int i = 0; i < Nk - 1; i++) {
    int swapIdx = Nk - 2 - i;
    schur.swapEvals(swapIdx);
    std::cout << GridLogMessage << "S after swap of index " << swapIdx << " = " << std::endl << schur.getMatrixS() << std::endl;
    std::cout << "Schur decomp still holds? " << schur.checkDecomposition() << std::endl;
  }

  std::cout << GridLogMessage << "Testing Schur reorder" << std::endl;
  schur.schurReorder(Nk);
  std::cout << GridLogMessage << "S after reorder = " << std::endl << schur.getMatrixS() << std::endl;
  std::cout << "Schur decomp still holds? " << schur.checkDecomposition() << std::endl;

}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=16;

//   GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
//  std::vector<int> lat_size {32, 32, 32, 32};
//  std::cout << "Lattice size: " << lat_size << std::endl;
  GridCartesian * UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								          GridDefaultSimd(Nd,vComplex::Nsimd()),
								          GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

//  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
//  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
  GridCartesian         * FGrid   = UGrid;
  GridRedBlackCartesian * FrbGrid = UrbGrid;

  // Construct a coarsened grid
  // poare TODO: replace this with the following line?
  Coordinate clatt = GridDefaultLatt();
//   Coordinate clatt = GridDefaultLatt();              // [PO] initial line before I edited it
  for(int d=0;d<clatt.size();d++){
  std::cout << GridLogMessage<< clatt[d] <<std::endl;
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

  LatticeFermion result(FGrid); result=Zero();
  LatticeFermion    ref(FGrid); ref=Zero();
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid);

  FieldMetaData header;
  std::string file("config");
//  std::string file("Users/patrickoare/libraries/PETSc-Grid/ckpoint_lat.4000");
  NerscIO::readConfiguration(Umu,header,file);

  LanczosParameters LanParams;
  {
    XmlReader  HMCrd("LanParams.xml");
    read(HMCrd,"LanczosParameters",LanParams);
  }

  std::cout << GridLogMessage<< LanParams <<std::endl;
  {
    XmlWriter HMCwr("LanParams.xml.out");
    write(HMCwr,"LanczosParameters",LanParams);
  }


  RealD mass=0.01;
  RealD M5=1.8;

  // PowerMethod<LatticeFermion> PM; PM(PVdagM, src);
  int Nm = 50;
  int Nk = 12; 
  int Np = 38; 
  // int Nk = Nm+1;     // if just running once
  int maxIter = 10000;
  int Nstop = 10;
  RealD resid = 1.0e-5;

  std::vector<Complex> boundary = {1,1,1,-1};
  WilsonOp::ImplParams Params(boundary);

//  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
//  DomainWallFermionD Dpv(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,M5);

  mass=LanParams.mass;
  std::cout << GridLogIRL<< "mass "<<mass<<std::endl;
  WilsonOp WilsonOperator(Umu,*UGrid,*UrbGrid,mass,Params);

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

//  typedef PVdagMLinearOperator<DomainWallFermionD,LatticeFermionD> PVdagM_t;
//  typedef ShiftedPVdagMLinearOperator<DomainWallFermionD,LatticeFermionD> ShiftedPVdagM_t;
//  typedef ShiftedComplexPVdagMLinearOperator<DomainWallFermionD,LatticeFermionD> ShiftedComplexPVdagM_t;
//  PVdagM_t PVdagM(Ddwf, Dpv);
//  ShiftedPVdagM_t ShiftedPVdagM(0.1,Ddwf,Dpv);
//  SquaredLinearOperator<DomainWallFermionD, LatticeFermionD> Dsq (Ddwf);
//  NonHermitianLinearOperator<DomainWallFermionD, LatticeFermionD> DLinOp (Ddwf);


  NonHermitianLinearOperator<WilsonOp,FermionField> Dwilson(WilsonOperator); /// <-----
  InvertNonHermitianLinearOperator<WilsonOp,FermionField> Iwilson(WilsonOperator); /// <-----
  MdagMLinearOperator<WilsonOp,FermionField> HermOp(WilsonOperator); /// <-----
  Gamma5HermitianLinearOperator <WilsonOp,LatticeFermion> HermOp2(WilsonOperator); /// <----

  // PowerMethod<LatticeFermion> PM; PM(PVdagM, src);
  resid=LanParams.resid;
  Nstop=LanParams.Nstop;
  Nk=LanParams.Nk;
  Np=LanParams.Np;
  maxIter=LanParams.maxIter;
  Nm = Nk + Np;
  int Nu=16;
  std::vector<LatticeFermion> src(Nu,FGrid); 
  for(int i=0;i<Nu;i++){
     random(RNG5,src[i]);
#if 0
     LatticeFermion src_e(FrbGrid);
     pickCheckerboard(Even, src_e, src[i]);
     src_e=Zero();	
     setCheckerboard(src[i],src_e);
     std::cout << GridLogMessage << "src["<<i<<"] even sites zeroed " << std::endl;
#endif
  }

  if(LanParams.ReadEvec) {
    std::string evecs_file="evec_in";
    std::cout << GridLogIRL<< "Reading evecs from "<<evecs_file<<std::endl;
    emptyUserRecord record;
    Grid::ScidacReader RD;
    RD.open(evecs_file);
    RD.readScidacFieldRecord(src[0],record);
    RD.close();
  }

  Coordinate origin ({0,0,0,0});
  auto tmpSrc = peekSite(src[0], origin);
  std::cout << "[DEBUG] Source at origin = " <<  tmpSrc << std::endl;
  LatticeFermion src2 = src[0];

  // Run KrylovSchur and Arnoldi on a Hermitian matrix
  std::cout << GridLogMessage << "Running Krylov Schur" << std::endl;
    RealD shift=LanParams.shift;
#if 0
    KrylovSchur KrySchur (Dwilson, UGrid, resid,EvalImNormSmall);
//    KrySchur(src[0], maxIter, Nm, Nk, Nstop);
    KrySchur(src[0], maxIter, Nm, Nk, Nstop,&shift);
    std::cout << GridLogMessage << "KrylovSchur evec.size= " << KrySchur.evecs.size()<< std::endl;
#else
    int Nblock=4;
    Nblock=LanParams.Nblock;
    bool if_verify=false;
    if(LanParams.verify) if_verify=true;
//    BlockKrylovSchur KrySchur (Dwilson, UGrid, resid,EvalImNormSmall);
    HarmonicBlockKrylovSchur KrySchur (Dwilson, UGrid, resid,shift,EvalNormSmall);
    KrySchur(src, maxIter, Nm, Nk, Nstop,Nblock,true,if_verify);
    std::cout << GridLogMessage << "BlockKrylovSchur evec.size= " << KrySchur.evecs.size()<< std::endl;
#endif

  src[0]=KrySchur.evecs[0];
  std::cout << GridLogMessage << "KrySchur.evecs= "<< KrySchur.evecs.size() <<std::endl;
  for (int i=1;i<Nk;i++) src[0]+=KrySchur.evecs[i];
  for (int i=0;i<Nstop;i++) 
  {
	std::string evfile ("./evec_"+std::to_string(mass)+"_"+std::to_string(i));
        auto evdensity = localInnerProduct(KrySchur.evecs[i],KrySchur.evecs[i] );
        writeFile(evdensity,evfile);

  }

  {
        std::string evfile ("./evec_"+std::to_string(mass)+"_sum");
//        auto evdensity = localInnerProduct(evec[i],evec[i] );
        writeFile(src[0],evfile);
  }


  /*
  std::cout << GridLogMessage << "Running Arnoldi" << std::endl;
  // Arnoldi Arn (Dsq, FGrid, 1e-8);
  Arnoldi Arn (DLinOp, FGrid, 1e-8);
  testSchurFromHess<LatticeFermion>(Arn, src, 10, 6, 4);

  Arnoldi Arn2 (DLinOp, FGrid, 1e-8);
  testSchurFromHess<LatticeFermion>(Arn2, src, 16, 12, 8);
  */

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;

  Grid_finalize();
  return 0;
}
