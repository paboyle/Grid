/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/solver/Test_zMADWF_prec.cc

Copyright (C) 2015

Author: Christopher Kelly <ckelly@phys.columbia.edu>

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

//This test computes the zMobius approximation to the Mobius action and uses it within the MADWF context to accelerate an inversion
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

struct TestParams{
  bool load_config;
  std::string config_file;

  double mass;

  std::string outer_precon;
  std::string inner_precon;
  
  int Ls_outer;
  double b_plus_c_outer;
  double resid_outer;
  
  int Ls_inner;
  double b_plus_c_inner; //irrelevant for ZMobius
  double resid_inner;
  bool zmobius_inner;
  double lambda_max; //upper bound of H_T eigenvalue range required to generate zMobius approximation
  
  TestParams(): load_config(true), config_file("ckpoint_lat.1000"), mass(0.01),
		Ls_outer(24), b_plus_c_outer(2.0), resid_outer(1e-8),
		Ls_inner(12), b_plus_c_inner(1.0), resid_inner(1e-8), zmobius_inner(true), lambda_max(1.42), outer_precon("Standard"), inner_precon("Standard")
  {}
  
  void write(const std::string &file) const{
    XmlWriter wr(file);
#define DOIT(A) wr.writeDefault(#A, A)
    DOIT(load_config);
    DOIT(config_file);
    DOIT(mass);
    DOIT(outer_precon);
    DOIT(inner_precon);
    DOIT(Ls_outer);
    DOIT(b_plus_c_outer);
    DOIT(resid_outer);
    DOIT(Ls_inner);
    DOIT(b_plus_c_inner);
    DOIT(resid_inner);
    DOIT(zmobius_inner);
    DOIT(lambda_max);
#undef DOIT
  }
  void read(const std::string &file){
    XmlReader rd(file);
#define DOIT(A) rd.readDefault(#A, A)
    DOIT(load_config);
    DOIT(config_file);
    DOIT(mass);
    DOIT(outer_precon);
    DOIT(inner_precon);
    DOIT(Ls_outer);
    DOIT(b_plus_c_outer);
    DOIT(resid_outer);
    DOIT(Ls_inner);
    DOIT(b_plus_c_inner);
    DOIT(resid_inner);
    DOIT(zmobius_inner);
    DOIT(lambda_max);
#undef DOIT
  }
};

struct RunParamsPrecStd{
  typedef SchurRedBlackDiagMooeeSolve<LatticeFermionD> SchurSolverType;

  template<typename Action>
  using HermOpType = SchurDiagMooeeOperator<Action, LatticeFermionD>;
};

struct RunParamsPrecDiagTwo{
  typedef SchurRedBlackDiagTwoSolve<LatticeFermionD> SchurSolverType;

  template<typename Action>
  using HermOpType = SchurDiagTwoOperator<Action, LatticeFermionD>;
};


struct CGincreaseTol : public MADWFinnerIterCallbackBase{
  ConjugateGradient<LatticeFermionD> &cg_inner;  
  RealD outer_resid;

  CGincreaseTol(ConjugateGradient<LatticeFermionD> &cg_inner,
		RealD outer_resid): cg_inner(cg_inner), outer_resid(outer_resid){}
  
  void operator()(const RealD current_resid){
    std::cout << "CGincreaseTol with current residual " << current_resid << " changing inner tolerance " << cg_inner.Tolerance << " -> ";
    while(cg_inner.Tolerance < current_resid) cg_inner.Tolerance *= 2;    
    //cg_inner.Tolerance = outer_resid/current_resid;
    std::cout << cg_inner.Tolerance << std::endl;
  }
};

template<typename RunParamsOuter, typename RunParamsInner>
void run(const TestParams &params){
  RealD bmc = 1.0; //use Shamir kernel
  std::vector<ComplexD> gamma_inner;

  std::cout << "Compute parameters" << std::endl;
  if(params.zmobius_inner){
    Approx::computeZmobiusGamma(gamma_inner, params.b_plus_c_inner, params.Ls_inner, params.b_plus_c_outer, params.Ls_outer, params.lambda_max);
  }else{
    Approx::zolotarev_data *zdata = Approx::higham(1.0,params.Ls_inner);
    gamma_inner.resize(params.Ls_inner);
    for(int s=0;s<params.Ls_inner;s++) gamma_inner[s] = zdata->gamma[s];
    Approx::zolotarev_free(zdata);
  }
  std::cout << "gamma:\n";
  for(int s=0;s<params.Ls_inner;s++) std::cout << s << " " << gamma_inner[s] << std::endl;


  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);


  GridCartesian* FGrid_outer = SpaceTimeGrid::makeFiveDimGrid(params.Ls_outer, UGrid);
  GridCartesian* FGrid_inner = SpaceTimeGrid::makeFiveDimGrid(params.Ls_inner, UGrid);

  GridRedBlackCartesian* FrbGrid_outer = SpaceTimeGrid::makeFiveDimRedBlackGrid(params.Ls_outer, UGrid);
  GridRedBlackCartesian* FrbGrid_inner = SpaceTimeGrid::makeFiveDimRedBlackGrid(params.Ls_inner, UGrid);


  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});

  GridParallelRNG RNG5_outer(FGrid_outer);
  RNG5_outer.SeedFixedIntegers(seeds5);

  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD src4(UGrid); random(RNG4,src4);

  LatticeFermionD result_outer(FGrid_outer);
  result_outer = Zero();
  LatticeGaugeFieldD Umu(UGrid);

  if(params.load_config){
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, params.config_file);

    for(int i=0;i<Nd;i++){
      assert(header.dimension[i] == GridDefaultLatt()[i]);
    }
  }else{    
    SU3::HotConfiguration(RNG4, Umu);
  }
    
  std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt()
            << "   Ls: " << params.Ls_outer << std::endl;

  RealD M5 = 1.8;

  RealD b_outer = (params.b_plus_c_outer + bmc)/2.;
  RealD c_outer = (params.b_plus_c_outer - bmc)/2.;

  RealD b_inner = (params.b_plus_c_inner + bmc)/2.;
  RealD c_inner = (params.b_plus_c_inner - bmc)/2.;

  MobiusFermionD D_outer(Umu, *FGrid_outer, *FrbGrid_outer, *UGrid, *UrbGrid, params.mass, M5, b_outer, c_outer);
  ZMobiusFermionD D_inner(Umu, *FGrid_inner, *FrbGrid_inner, *UGrid, *UrbGrid, params.mass, M5, gamma_inner, b_inner, c_inner);

  LatticeFermionD src_outer(FGrid_outer);
  D_outer.ImportPhysicalFermionSource(src4,src_outer); //applies D_- 

  //Solve using a regular even-odd preconditioned CG for the Hermitian operator
  //M y = x
  //Mprec y'_o = x'_o     where Mprec = Doo - Doe Dee^-1 Deo    and  x'_o = -Doe Dee^-1 x_e + x_o
  //y_o = y'_o

  //(Mprec^dag Mprec) y'_o = Mprec^dag x'_o 
  //y'_o = (Mprec^dag Mprec)^-1 Mprec^dag x'_o 

  //We can get Mprec^dag x'_o from x_o  from SchurRedBlackDiagMooeeSolve::RedBlackSource
  ConjugateGradient<LatticeFermionD> CG_outer(params.resid_outer, 10000);
  typename RunParamsOuter::SchurSolverType SchurSolver_outer(CG_outer);
  
  LatticeFermionD tmp_e_outer(FrbGrid_outer);
  LatticeFermionD src_o_outer(FrbGrid_outer);
  SchurSolver_outer.RedBlackSource(D_outer, src_outer, tmp_e_outer, src_o_outer);
  
  LatticeFermionD result_o_outer(FrbGrid_outer);
  result_o_outer = Zero();

  GridStopWatch CGTimer;
  
  typename RunParamsOuter::template HermOpType<MobiusFermionD> HermOpEO_outer(D_outer);

  CGTimer.Start();
  CG_outer(HermOpEO_outer, src_o_outer, result_o_outer);
  CGTimer.Stop();

  std::cout << GridLogMessage << "Total outer CG time : " << CGTimer.Elapsed()
            << std::endl;

  CGTimer.Reset();

  //Solve for y using MADWF with internal preconditioning

  //typedef PauliVillarsSolverRBprec<LatticeFermionD, typename RunParamsOuter::SchurSolverType> PVtype;
  //PVtype PV_outer(SchurSolver_outer);

  typedef PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD> PVtype;
  PVtype PV_outer(Umu, CG_outer);

  ConjugateGradient<LatticeFermionD> CG_inner(params.resid_inner, 10000, 0);

  CGincreaseTol update(CG_inner, params.resid_outer);

  typename RunParamsInner::SchurSolverType SchurSolver_inner(CG_inner);

  ZeroGuesser<LatticeFermion> Guess;
  MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, typename RunParamsInner::SchurSolverType, ZeroGuesser<LatticeFermion> > madwf(D_outer, D_inner, PV_outer, SchurSolver_inner, Guess, params.resid_outer, 100, &update);
  
  LatticeFermionD result_MADWF(FGrid_outer);
  result_MADWF = Zero();

  CGTimer.Start();
  madwf(src4, result_MADWF);
  CGTimer.Stop();

  LatticeFermionD result_o_MADWF(FrbGrid_outer);
  pickCheckerboard(Odd, result_o_MADWF, result_MADWF);

  std::cout << GridLogMessage << "Total MADWF time : " << CGTimer.Elapsed()
            << std::endl;

  LatticeFermionD diff = result_o_MADWF - result_o_outer;
  std::cout <<GridLogMessage<< "Odd-parity MADWF result norm " << norm2(result_o_MADWF) 
	    << " Regular result norm " << norm2(result_o_outer) 
	    << " Norm of diff " << norm2(diff)<<std::endl;


  //std::cout << GridLogMessage << "######## Dhop calls summary" << std::endl;
  //D_outer.Report();
}



int main(int argc, char** argv) {
  std::cout << "Init" << std::endl;
  Grid_init(&argc, &argv);

  TestParams params;
  
  if( GridCmdOptionExists(argv,argv+argc,"--params") ){
    std::string pfile = GridCmdOptionPayload(argv,argv+argc,"--params");
    if(pfile == "TEMPLATE"){
      params.write("params.templ");
      return 0;
    }else{
      params.read(pfile);
    }
  }

  if(params.outer_precon == "Standard" && params.inner_precon == "Standard" ){
    run<RunParamsPrecStd, RunParamsPrecStd>(params);
  }else if(params.outer_precon == "DiagTwo" && params.inner_precon == "Standard"){
    run<RunParamsPrecDiagTwo, RunParamsPrecStd>(params);
  }else if(params.outer_precon == "Standard" && params.inner_precon == "DiagTwo"){
    run<RunParamsPrecStd, RunParamsPrecDiagTwo>(params);
  }else if(params.outer_precon == "DiagTwo" && params.inner_precon == "DiagTwo"){
    run<RunParamsPrecDiagTwo, RunParamsPrecDiagTwo>(params);
  }else assert(0);

  Grid_finalize();
}
