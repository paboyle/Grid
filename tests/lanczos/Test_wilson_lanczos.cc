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
#include <Grid/parallelIO/IldgIOtypes.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczos.h>


using namespace std;
using namespace Grid;
 ;

typedef WilsonFermionD FermionOp;
typedef typename WilsonFermionD::FermionField FermionField;


RealD AllZero(RealD x) { return 0.; }

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


namespace Grid {

struct LanczosParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
		  		RealD, mass , 
		  		RealD, mstep , 
				Integer, Nstop,
                                Integer, Nk,
                                Integer, Np,
                                Integer, ReadEvec,
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

  int Ndir=4;
  auto mpi_layout  = GridDefaultMpi();
  std::vector<int> nblock(4,1);
  std::vector<int> mpi_split(4,1);
//Interested in avoiding degeneracy only for now
  nblock[3]=2;

  int mrhs=1;
  for(int i =0;i<Ndir;i++){
      mpi_split[i] = mpi_layout[i] / nblock[i];
      mrhs *= nblock[i];
  }


  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
      GridDefaultMpi());

  GridCartesian * SGrid = new GridCartesian(GridDefaultLatt(),
                                                    GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                    mpi_split,
                                                    *UGrid);

  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

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
//  SU<Nc>::ColdConfiguration(Umu);

  FieldMetaData header;
  std::string file("./config");

//  int precision32 = 0;
//  int tworow      = 0;
//  NerscIO::writeConfiguration(Umu,file,tworow,precision32);
  NerscIO::readConfiguration(Umu,header,file);

/*
  std::vector<LatticeColourMatrix> U(4, UGrid);
  for (int mu = 0; mu < Nd; mu++) {
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
  }
*/

  int Nstop = 10;
  int Nu = 1;
  int Nk = 20;
  int Np = 80;
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
  Nstop=LanParams.Nstop;
  Nu = mrhs;
  Nk=LanParams.Nk;
  Np=LanParams.Np;
  Nm = Nk + Np;

//  FermionField src(FGrid);
  std::vector<FermionField> src(Nu,FGrid);
  for(int i =0;i<Nu;i++) gaussian(RNG5, src[i]);

  if(LanParams.ReadEvec) {
    std::string evecs_file="evec_in";
    std::cout << GridLogIRL<< "Reading evecs from "<<evecs_file<<std::endl;
    emptyUserRecord record;
    Grid::ScidacReader RD;
    RD.open(evecs_file);
    RD.readScidacFieldRecord(src[0],record);
    RD.close();
  }

  std::vector<Complex> boundary = {1,1,1,-1};
//  std::vector<Complex> boundary = {1,1,1,1};
  FermionOp::ImplParams Params(boundary);

  GridCartesian         * SFGrid   = SGrid;
  GridRedBlackCartesian * SFrbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(SFGrid);
//  GridRedBlackCartesian * SFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(JP.Ls,SGrid);

  LatticeGaugeField s_Umu(SGrid);
  Grid_split  (Umu,s_Umu);



while ( mass > - 2.0){
  FermionOp WilsonOperator(Umu,*FGrid,*FrbGrid,mass,Params);
  MdagMLinearOperator<FermionOp,FermionField> HermOp(WilsonOperator); /// <-----
  FermionOp WilsonSplit(s_Umu,*SFGrid,*SFrbGrid,mass,Params);
  MdagMLinearOperator<FermionOp,FermionField> SHermOp(WilsonSplit); /// <-----
  //SchurDiagTwoOperator<FermionOp,FermionField> HermOp(WilsonOperator);
  Gamma5HermitianLinearOperator <FermionOp,LatticeFermion> HermOp2(WilsonOperator); /// <-----

  std::vector<double> Coeffs{0, 1.};
  Polynomial<FermionField> PolyX(Coeffs);
//  Chebyshev<FermionField> Cheby(0.5, 60., 31);
//                                  RealD, ChebyLow,
//                                RealD, ChebyHigh,
//                                Integer, ChebyOrder)

  Chebyshev<FermionField> Cheby(LanParams.ChebyLow,LanParams.ChebyHigh,LanParams.ChebyOrder);

  FunctionHermOp<FermionField> OpCheby(Cheby,HermOp);
     PlainHermOp<FermionField> Op     (HermOp);
     PlainHermOp<FermionField> Op2     (HermOp2);

//  ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby, Op2, Nstop, Nk, Nm, resid, MaxIt);
//  SimpleLanczos<FermionField> IRL(Op,Nstop, Nk, Nm, resid, MaxIt);
    ImplicitlyRestartedBlockLanczos<FermionField> IRBL(HermOp, SHermOp,
                                                     FrbGrid,SFrbGrid,mrhs,
                                                     Cheby,
                                                     Nstop, Nstop*2,
                                                     Nu, Nk, Nm,
                                                     resid, MaxIt,
                                                     IRBLdiagonaliseWithEigen);
  IRBL.split_test=1;
  std::vector<RealD> eval(Nm);
  std::vector<FermionField> evec(Nm, FGrid);
  for (int i = 0; i < 1; i++) {
    std::cout << i << " / " << Nm << " grid pointer " << evec[i].Grid()
              << std::endl;
  };

  int Nconv;
//  IRL.calc(eval, evec, src, Nconv);
  IRBL.calc(eval, evec, src, Nconv,LanczosType::irbl);

  std::cout << mass <<" : " << eval << std::endl;

  Gamma g5(Gamma::Algebra::Gamma5) ;
  ComplexD dot;
  FermionField tmp(FGrid);
  FermionField sav(FGrid);
  sav=evec[0];
  for (int i = 0; i < Nstop ; i++) {
    tmp = g5*evec[i];
    dot = innerProduct(tmp,evec[i]);
    std::cout << mass << " : " << eval[i]  << " " << real(dot) << " " << imag(dot)  << std::endl ;
//    if ( i<1)
    {
	std::string evfile ("./evec_"+std::to_string(mass)+"_"+std::to_string(i));
        auto evdensity = localInnerProduct(evec[i],evec[i] );
	writeFile(evdensity,evfile);
    }
    if (i>0) sav += evec[i];
  }
  {
	std::string evfile ("./evec_"+std::to_string(mass)+"_sum");
//        auto evdensity = localInnerProduct(evec[i],evec[i] );
	writeFile(sav,evfile);
  }
  for(int i =0;i<Nu;i++) src[i]=evec[i];
  for(int i=Nu;i<Nstop;i++) src[i%Nu] +=evec[i];
//  src  = evec[0]+evec[1]+evec[2];
//  src  += evec[3]+evec[4]+evec[5];
//  src  += evec[6]+evec[7]+evec[8];
  mass += LanParams.mstep;
}

  Grid_finalize();
}
