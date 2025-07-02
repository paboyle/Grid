/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_dwf_G5R5.cc

Copyright (C) 2015

Author: Chulwoo Jung <chulwoo@bnl.gov>
From Duo and Bob's Chirality study

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

//typedef WilsonFermionD FermionOp;
typedef DomainWallFermionD FermionOp;
typedef typename DomainWallFermionD::FermionField FermionField;

template <class T> void writeFile(T& in, std::string const fname){  
#ifdef HAVE_LIME
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

RealD AllZero(RealD x) { return 0.; }

namespace Grid {

struct LanczosParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
		  		RealD, mass , 
		  		RealD, M5 , 
	  			Integer, Ls,
	  			Integer, Nstop,
	  			Integer, Nk,
	  			Integer, Np,
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

  int Ls=16;
  RealD M5=1.8;
  RealD mass = 0.01;

  mass=LanParams.mass;
  Ls=LanParams.Ls;
  M5=LanParams.M5;

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid =
      SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
//  GridCartesian* FGrid = UGrid;
//  GridRedBlackCartesian* FrbGrid = UrbGrid;
  GridCartesian * FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);
//  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n", UGrid, UrbGrid, FGrid, FrbGrid);

  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5rb(FrbGrid); RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid);

  FieldMetaData header;
  std::string file("./config");

  int precision32 = 0;
  int tworow      = 0;
  NerscIO::readConfiguration(Umu,header,file);

/*
  std::vector<LatticeColourMatrix> U(4, UGrid);
  for (int mu = 0; mu < Nd; mu++) {
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
  }
*/
  int Nk = 20;
  int Nstop = Nk;
  int Np = 80;

  Nstop=LanParams.Nstop;
  Nk=LanParams.Nk;
  Np=LanParams.Np;

  int Nm = Nk + Np;
  int MaxIt = 100;
  RealD resid = 1.0e-4;


//while ( mass > - 5.0){
  FermionOp Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  MdagMLinearOperator<FermionOp,FermionField> HermOp(Ddwf); /// <-----
//  Gamma5HermitianLinearOperator <FermionOp,LatticeFermion> HermOp2(WilsonOperator); /// <-----
  Gamma5R5HermitianLinearOperator<FermionOp, LatticeFermion> G5R5Herm(Ddwf);
//  Gamma5R5HermitianLinearOperator
  std::vector<double> Coeffs{0, 1.};
  Polynomial<FermionField> PolyX(Coeffs);

  Chebyshev<FermionField> Cheby(LanParams.ChebyLow,LanParams.ChebyHigh,LanParams.ChebyOrder);

  FunctionHermOp<FermionField> OpCheby(Cheby,HermOp);
  PlainHermOp<FermionField> Op     (HermOp);
  PlainHermOp<FermionField> Op2     (G5R5Herm);

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

  std::cout << mass <<" : " << eval        << std::endl;
  std::cout << " #evecs "   << evec.size() << std::endl;
  std::cout << " Nconv  "   << Nconv       << std::endl;
  std::cout << " Nm     "   << Nm          << std::endl;
  if ( Nconv > evec.size() ) Nconv = evec.size();
  
#if 0
  Gamma g5(Gamma::Algebra::Gamma5) ;
  ComplexD dot;
  FermionField tmp(FGrid);
//  RealD eMe,eMMe;
  for (int i = 0; i < Nstop ; i++) {
//    tmp = g5*evec[i];
    dot = innerProduct(evec[i],evec[i]);
//    G5R5(tmp,evec[i]);
    G5R5Herm.HermOpAndNorm(evec[i],tmp,eMe,eMMe);
    std::cout <<"Norm "<<M5<<" "<< mass << " : " << i << " " << real(dot) << " " << imag(dot)  << " "<< eMe << " " <<eMMe<< std::endl ;
    for (int j = 0; j < Nstop ; j++) {
      dot = innerProduct(tmp,evec[j]);
      std::cout <<"G5R5 "<<M5<<" "<< mass << " : " << i << " " <<j<<" " << real(dot) << " " << imag(dot)  << std::endl ;
    }
  }
//  src  = evec[0]+evec[1]+evec[2];
//  mass += -0.1;
#endif

  //**********************************************************************
  //orthogonalization
  //calculat the matrix
  cout << "Start orthogonalization " << endl;
  cout << "calculate the matrix element" << endl;
  vector<LatticeFermion> G5R5Mevec(Nconv, FGrid);
  vector<LatticeFermion> finalevec(Nconv, FGrid);
  vector<RealD> eMe(Nconv), eMMe(Nconv);
  for(int i = 0; i < Nconv; i++){
    cout << "calculate the matrix element["<<i<<"]" << endl;
    G5R5Herm.HermOpAndNorm(evec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
  }
  cout << "Re<evec, G5R5M(evec)>: " << endl;
  cout << eMe << endl;
  cout << "<G5R5M(evec), G5R5M(evec)>" << endl;
  cout << eMMe << endl;
  vector<vector<ComplexD>> VevecG5R5Mevec(Nconv);
  Eigen::MatrixXcd evecG5R5Mevec = Eigen::MatrixXcd::Zero(Nconv, Nconv);
  for(int i = 0; i < Nconv; i++){
    VevecG5R5Mevec[i].resize(Nconv);
    for(int j = 0; j < Nconv; j++){
      VevecG5R5Mevec[i][j] = innerProduct(evec[i], G5R5Mevec[j]);
      evecG5R5Mevec(i, j) = VevecG5R5Mevec[i][j];
    }
  }
  //calculate eigenvector
  cout << "Eigen solver" << endl;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(evecG5R5Mevec);
  vector<RealD> eigeneval(Nconv);
  vector<vector<ComplexD>> eigenevec(Nconv);
  for(int i = 0; i < Nconv; i++){
    eigeneval[i] = eigensolver.eigenvalues()[i];
    eigenevec[i].resize(Nconv);
    for(int j = 0; j < Nconv; j++){
      eigenevec[i][j] = eigensolver.eigenvectors()(i, j);
    }
  }
  //rotation
  cout << "Do rotation" << endl;
  for(int i = 0; i < Nconv; i++){
    finalevec[i] = finalevec[i] - finalevec[i];
    for(int j = 0; j < Nconv; j++){
      finalevec[i] = eigenevec[j][i]*evec[j] + finalevec[i];
    }
  }
  //normalize again;
  for(int i = 0; i < Nconv; i++){
    RealD tmp_RealD = norm2(finalevec[i]);
    tmp_RealD = 1./pow(tmp_RealD, 0.5);
    finalevec[i] = finalevec[i]*tmp_RealD;
  }

  //check
  for(int i = 0; i < Nconv; i++){
    G5R5Herm.HermOpAndNorm(finalevec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
  }

  //**********************************************************************
  //sort the eigenvectors
  vector<LatticeFermion> finalevec_copy(Nconv, FGrid);
  for(int i = 0; i < Nconv; i++){
    finalevec_copy[i] = finalevec[i];
  }
  vector<RealD> eMe_copy(eMe);
  for(int i = 0; i < Nconv; i++){
    eMe[i] = fabs(eMe[i]);
    eMe_copy[i] = eMe[i];
  }
  sort(eMe_copy.begin(), eMe_copy.end());
  for(int i = 0; i < Nconv; i++){
    for(int j = 0; j < Nconv; j++){
      if(eMe[j] == eMe_copy[i]){
        finalevec[i] = finalevec_copy[j];
      }
    }
  }
  for(int i = 0; i < Nconv; i++){
    G5R5Herm.HermOpAndNorm(finalevec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
  }
  cout << "Re<evec, G5R5M(evec)>: " << endl;
  cout << eMe << endl;
  cout << "<G5R5M(evec), G5R5M(evec)>" << endl;
  cout << eMMe << endl;

  

//  vector<LatticeFermion> finalevec(Nconv, FGrid);
// temporary, until doing rotation
//  for(int i = 0; i < Nconv; i++)
//	  finalevec[i]=evec[i];
  //**********************************************************************
  //calculate chirality matrix
  vector<LatticeFermion> G5evec(Nconv, FGrid);
  vector<vector<ComplexD>> chiral_matrix(Nconv);
  vector<vector<RealD>> chiral_matrix_real(Nconv);
  for(int i = 0; i < Nconv; i++){
//    G5evec[i] = G5evec[i] - G5evec[i];
    G5evec[i] = Zero();
    for(int j = 0; j < Ls/2; j++){
      axpby_ssp(G5evec[i], 1., finalevec[i], 0., G5evec[i], j, j);
    }
    for(int j = Ls/2; j < Ls; j++){
      axpby_ssp(G5evec[i], -1., finalevec[i], 0., G5evec[i], j, j);
    }
  }
  
  for(int i = 0; i < Nconv; i++){
    Ddwf.M(finalevec[i], G5R5Mevec[i]);
    for(int j = 0; j < Nconv; j++){
      std::cout << "<"<<j<<"|Ddwf|"<<i<<"> = "<<innerProduct(finalevec[j],G5R5Mevec[i])<<std::endl;
    }
  }
  for(int i = 0; i < Nconv; i++){
    RealD t1,t2;
    G5R5Herm.HermOpAndNorm(finalevec[i], G5R5Mevec[i], t1, t2);
    for(int j = 0; j < Nconv; j++){
      std::cout << "<"<<j<<"|G5R5 M|"<<i<<"> = "<<innerProduct(finalevec[j],G5R5Mevec[i])<<std::endl;
    }
  }
  
  for(int i = 0; i < Nconv; i++){
    chiral_matrix_real[i].resize(Nconv);
    chiral_matrix[i].resize(Nconv);

    std::string evfile("./evec_density");
    evfile = evfile+"_"+std::to_string(i);
    auto evdensity = localInnerProduct(finalevec[i],finalevec[i] );
    writeFile(evdensity,evfile);

    for(int j = 0; j < Nconv; j++){
      chiral_matrix[i][j] = innerProduct(finalevec[i], G5evec[j]);
      std::cout <<" chiral_matrix_real signed "<<i<<" "<<j<<" "<< chiral_matrix_real[i][j] << std::endl;
      chiral_matrix_real[i][j] = abs(chiral_matrix[i][j]);
      std::cout <<" chiral_matrix_real "<<i<<" "<<j<<" "<< chiral_matrix_real[i][j] << std::endl;
      if ( chiral_matrix_real[i][j] > 0.8 ) {
	auto g5density = localInnerProduct(finalevec[i], G5evec[j]);
	std::string chfile("./chiral_density_");
	chfile = chfile +std::to_string(i)+"_"+std::to_string(j);
	writeFile(g5density,chfile);
      }
    }
  }
  for(int i = 0; i < Nconv; i++){
    if(chiral_matrix[i][i].real() < 0.){
      chiral_matrix_real[i][i] = -1. * chiral_matrix_real[i][i];
    }
  }

  FILE *fp = fopen("lego-plot.py","w"); assert(fp!=NULL);
#define PYTHON_LINE(A)  fprintf(fp,A"\n");
  PYTHON_LINE("import matplotlib.pyplot as plt");
  PYTHON_LINE("import numpy as np");
  PYTHON_LINE("");
  PYTHON_LINE("fig = plt.figure()");
  PYTHON_LINE("ax = fig.add_subplot(projection='3d')");
  PYTHON_LINE("");
  PYTHON_LINE("x, y = np.random.rand(2, 100) * 4");
  fprintf(fp,"hist, xedges, yedges = np.histogram2d(x, y, bins=%d, range=[[0, %d], [0, %d]])\n",Nconv,Nconv-1,Nconv-1);
  PYTHON_LINE("");
  PYTHON_LINE("# Construct arrays for the anchor positions of the 16 bars");
  PYTHON_LINE("xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25, indexing=\"ij\")");
  PYTHON_LINE("xpos = xpos.ravel()");
  PYTHON_LINE("ypos = ypos.ravel()");
  PYTHON_LINE("zpos = 0");
  PYTHON_LINE("");
  PYTHON_LINE("# Construct arrays with the dimensions for the 16 bars.");
  PYTHON_LINE("dx = dy = 0.5 * np.ones_like(zpos)");
  PYTHON_LINE("dz = np.array([");
  for(int i = 0; i < Nconv; i++){
    fprintf(fp,"\t[ ");
    for(int j = 0; j < Nconv; j++){
      fprintf(fp,"%lf ",chiral_matrix_real[i][j]);
      if(j<Nconv-1) fprintf(fp,",");
      else          fprintf(fp," ");
    }
    fprintf(fp,"]");
    if(i<Nconv-1) fprintf(fp,",\n");
    else          fprintf(fp,"\n");
  }
	      
  PYTHON_LINE("\t])");
  PYTHON_LINE("dz = dz.ravel()");
  PYTHON_LINE("ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')");
  PYTHON_LINE("plt.show()");
  fclose(fp);
  
  Grid_finalize();
}
