/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_dwf_lanczos.cc

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
 ;

typedef WilsonFermionD FermionOp;
typedef typename WilsonFermionD::FermionField FermionField;


RealD AllZero(RealD x) { return 0.; }

namespace Grid {

#if 0
template<typename Field>
class RationalHermOp : public LinearFunction<Field> {
public:
  using LinearFunction<Field>::operator();  
//  OperatorFunction<Field>   & _poly;
  LinearOperatorBase<Field> &_Linop;
  RealD _massDen, _massNum;
   
  FunctionHermOp(LinearOperatorBase<Field>& linop, RealD massDen,RealD massNum)
    :  _Linop(linop) ,_massDen(massDen),_massNum(massNum) {};
      
  void operator()(const Field& in, Field& out) {
//    _poly(_Linop,in,out);
  } 
};
#endif

template<class Matrix,class Field>
class InvG5LinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  RealD _num;
  RealD _Tol;
  Integer _MaxIt;
  Gamma g5;

public:
  InvG5LinearOperator(Matrix &Mat,RealD num): _Mat(Mat),_num(num), _Tol(1e-12),_MaxIt(10000), g5(Gamma::Algebra::Gamma5) {};

  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {
    assert(0);
    _Mat.Mdiag(in,out);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    assert(0);
    _Mat.Mdir(in,out,dir,disp);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
    assert(0);
    _Mat.MdirAll(in,out);
  };
  void Op     (const Field &in, Field &out){
    assert(0);
    _Mat.M(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    assert(0);
    _Mat.Mdag(in,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
     Field tmp(in.Grid());
     MdagMLinearOperator<Matrix,Field> denom(_Mat);
     ConjugateGradient<Field> CG(_Tol,_MaxIt); 
     _Mat.M(in,tmp);
     tmp += _num*in;
     _Mat.Mdag(tmp,out);
     CG(denom,out,tmp);
     out = g5*tmp;
  }
};


struct LanczosParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
		  		RealD, mass , 
				RealD, resid,
	  			RealD, ChebyLow,
	  			RealD, ChebyHigh,
	  			Integer, ChebyOrder)
//                                  Integer, StartTrajectory,
//                                  Integer, Trajectories, /* @brief Number of sweeps in this run */
//                                  bool, MetropolisTest,
//                                  Integer, NoMetropolisUntil,
//                                  std::string, StartingType,
//                                  Integer, SW,
//				  RealD, Kappa,
//                                  IntegratorParameters, MD)

  LanczosParameters() {
    ////////////////////////////// Default values
      mass = 0;
//    MetropolisTest    = true;
//    NoMetropolisUntil = 10;
//    StartTrajectory   = 0;
//    SW                = 2;
//    Trajectories      = 10;
//    StartingType      = "HotStart";
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

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid =
      SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian* FGrid = UGrid;
  GridRedBlackCartesian* FrbGrid = UrbGrid;
//  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n", UGrid, UrbGrid, FGrid, FrbGrid);

  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid);
  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5rb(FrbGrid);
  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid);
//  SU<Nc>::HotConfiguration(RNG4, Umu);

  FieldMetaData header;
  std::string file("./config");

  int precision32 = 0;
  int tworow      = 0;
//  NerscIO::writeConfiguration(Umu,file,tworow,precision32);
  NerscIO::readConfiguration(Umu,header,file);

/*
  std::vector<LatticeColourMatrix> U(4, UGrid);
  for (int mu = 0; mu < Nd; mu++) {
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
  }
*/

  int Nstop = 5;
  int Nk = 10;
  int Np = 90;
  int Nm = Nk + Np;
  int MaxIt = 10000;
  RealD resid = 1.0e-5;

  RealD mass = -1.0;

  LanczosParameters LanParams;
#if 1
  {
    XmlReader  HMCrd("LanParams.xml");
    read(HMCrd,"LanczosParameters",LanParams);
  }
#else
  {
    LanParams.mass = mass;
  }
#endif
  std::cout << GridLogMessage<< LanParams <<std::endl;
  { 
    XmlWriter HMCwr("LanParams.xml.out");
    write(HMCwr,"LanczosParameters",LanParams);
  }

  mass=LanParams.mass;
  resid=LanParams.resid;


while ( mass > - 5.0){
  FermionOp WilsonOperator(Umu,*FGrid,*FrbGrid,2.+mass);
  InvG5LinearOperator<FermionOp,LatticeFermion> HermOp(WilsonOperator,-2.); /// <-----
  //SchurDiagTwoOperator<FermionOp,FermionField> HermOp(WilsonOperator);
//  Gamma5HermitianLinearOperator <FermionOp,LatticeFermion> HermOp2(WilsonOperator); /// <-----

  std::vector<double> Coeffs{0, 0, 1.};
  Polynomial<FermionField> PolyX(Coeffs);
  Chebyshev<FermionField> Cheby(LanParams.ChebyLow,LanParams.ChebyHigh,LanParams.ChebyOrder);

       FunctionHermOp<FermionField> OpCheby(Cheby,HermOp);
//     InvHermOp<FermionField> Op(WilsonOperator,HermOp);
     PlainHermOp<FermionField> Op     (HermOp);
//     PlainHermOp<FermionField> Op2     (HermOp2);

  ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby, Op, Nstop, Nk, Nm, resid, MaxIt);

  std::vector<RealD> eval(Nm);
  FermionField src(FGrid);
  gaussian(RNG5, src);
  std::vector<FermionField> evec(Nm, FGrid);
  for (int i = 0; i < 1; i++) {
    std::cout << i << " / " << Nm << " grid pointer " << evec[i].Grid()
              << std::endl;
  };

  int Nconv;
  IRL.calc(eval, evec, src, Nconv);

  std::cout << mass <<" : " << eval << std::endl;

  Gamma g5(Gamma::Algebra::Gamma5) ;
  ComplexD dot;
  FermionField tmp(FGrid);
  for (int i = 0; i < Nstop ; i++) {
    tmp = g5*evec[i];
    dot = innerProduct(tmp,evec[i]);
    std::cout << mass << " : " << eval[i]  << " " << real(dot) << " " << imag(dot)  << std::endl ;
  }
  src  = evec[0]+evec[1]+evec[2];
  mass += -0.1;
}

  Grid_finalize();
}
