    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./forces/Test_gpdwf_force_1f_2f.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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

using namespace std;
using namespace Grid;

//Here we test the G-parity action and force between the 1f (doubled-lattice) and 2f approaches 


void copyConjGauge(LatticeGaugeFieldD &Umu_1f, const LatticeGaugeFieldD &Umu_2f, const int nu){
  GridBase* UGrid_2f = Umu_2f.Grid();
  GridBase* UGrid_1f = Umu_1f.Grid();

  Replicate(Umu_2f,Umu_1f);

  int L_2f = UGrid_2f->FullDimensions()[nu];
  int L_1f = UGrid_1f->FullDimensions()[nu]; 
  assert(L_1f == 2 * L_2f);

  //Coordinate grid for reference
  LatticeInteger xcoor_1f(UGrid_1f);
  LatticeCoordinate(xcoor_1f,nu);

  //Copy-conjugate the gauge field
  //First C-shift the lattice by Lx/2
  {
    LatticeGaugeField Umu_shift = conjugate( Cshift(Umu_1f,nu,L_2f) );
    Umu_1f = where( xcoor_1f >= Integer(L_2f), Umu_shift, Umu_1f );

    //We use the in built APBC 
    //Make the gauge field antiperiodic in nu-direction
    //decltype(PeekIndex<LorentzIndex>(Umu_1f,nu)) Unu(UGrid_1f);
    //Unu = PeekIndex<LorentzIndex>(Umu_1f,nu);
    //Unu = where(xcoor_1f == Integer(2*L_2f-1), -Unu, Unu);
    //PokeIndex<LorentzIndex>(Umu_1f,Unu,nu);
  }
}

template<typename FermionField2f, typename FermionField1f>
void convertFermion1f_from_2f(FermionField1f &out_1f, const FermionField2f &in_2f, const int nu, bool is_4d){
  GridBase* FGrid_1f = out_1f.Grid();
  GridBase* FGrid_2f = in_2f.Grid();

  int nuoff = is_4d ? 0 : 1;   //s in 0 direction

  Integer L_2f = FGrid_2f->FullDimensions()[nu+nuoff];
  Integer L_1f = FGrid_1f->FullDimensions()[nu+nuoff];
  assert(L_1f == 2 * L_2f);
  
  auto in_f0_2fgrid = PeekIndex<GparityFlavourIndex>(in_2f,0); //flavor 0 on 2f Grid
  FermionField1f in_f0_1fgrid(FGrid_1f);
  Replicate(in_f0_2fgrid, in_f0_1fgrid); //has flavor 0 on both halves

  auto in_f1_2fgrid = PeekIndex<GparityFlavourIndex>(in_2f,1); //flavor 1 on 2f Grid
  FermionField1f in_f1_1fgrid(FGrid_1f);
  Replicate(in_f1_2fgrid, in_f1_1fgrid); //has flavor 1 on both halves

  LatticeInteger xcoor_1f(FGrid_1f);
  LatticeCoordinate(xcoor_1f,nu+nuoff);
  
  out_1f = where(xcoor_1f < L_2f, in_f0_1fgrid, in_f1_1fgrid);
}

template<typename GparityAction, typename StandardAction>
class RatioActionSetupBase{
protected:
  TwoFlavourEvenOddRatioPseudoFermionAction<WilsonImplD> *pf_1f;
  TwoFlavourEvenOddRatioPseudoFermionAction<GparityWilsonImplD> *pf_2f;

  GparityAction* action_2f;
  GparityAction* action_PV_2f;
  StandardAction* action_1f;
  StandardAction* action_PV_1f;

  ConjugateGradient<typename StandardAction::FermionField> CG_1f;
  ConjugateGradient<typename GparityAction::FermionField> CG_2f;

  RatioActionSetupBase(): CG_1f(1.0e-8,10000), CG_2f(1.0e-8,10000){}

  void setupPseudofermion(){
    pf_1f = new TwoFlavourEvenOddRatioPseudoFermionAction<WilsonImplD>(*action_PV_1f, *action_1f, CG_1f, CG_1f);
    pf_2f = new TwoFlavourEvenOddRatioPseudoFermionAction<GparityWilsonImplD>(*action_PV_2f, *action_2f, CG_2f, CG_2f);
  }

public:
  GparityAction & action2f(){ return *action_2f; }
  StandardAction & action1f(){ return *action_1f; }

  void refreshAction(LatticeGaugeField &Umu_2f, typename GparityAction::FermionField &eta_2f,
		     LatticeGaugeField &Umu_1f, typename StandardAction::FermionField &eta_1f){  
    pf_1f->refresh(Umu_1f, eta_1f);
    pf_2f->refresh(Umu_2f, eta_2f);

    //Compare PhiOdd
    RealD norm_1f = norm2(pf_1f->getPhiOdd());
    RealD norm_2f = norm2(pf_2f->getPhiOdd());
    
    std::cout << "Test PhiOdd 2f: " << norm_2f << " 1f: " << norm_1f << std::endl;
  }

  void computeAction(RealD &S_2f, RealD &S_1f, LatticeGaugeField &Umu_2f, LatticeGaugeField &Umu_1f){
    S_1f = pf_1f->S(Umu_1f);
    S_2f = pf_2f->S(Umu_2f);
  }

  void computeDeriv(LatticeGaugeField &deriv_2f, LatticeGaugeField &deriv_1f, LatticeGaugeField &Umu_2f, LatticeGaugeField &Umu_1f){    
    pf_1f->deriv(Umu_1f, deriv_1f);
    pf_2f->deriv(Umu_2f, deriv_2f);
  }

};




template<typename GparityAction, typename StandardAction>
struct setupAction{};

template<>
struct setupAction<GparityWilsonTMFermionD, WilsonTMFermionD>: public RatioActionSetupBase<GparityWilsonTMFermionD, WilsonTMFermionD>{
  typedef GparityWilsonTMFermionD GparityAction;
  typedef WilsonTMFermionD StandardAction;
  
  setupAction(GridCartesian* UGrid_2f, GridRedBlackCartesian* UrbGrid_2f,  GridCartesian* FGrid_2f, GridRedBlackCartesian* FrbGrid_2f,
	      GridCartesian* UGrid_1f, GridRedBlackCartesian* UrbGrid_1f,  GridCartesian* FGrid_1f, GridRedBlackCartesian* FrbGrid_1f,
	      LatticeGaugeField &Umu_2f, LatticeGaugeField &Umu_1f, int nu): RatioActionSetupBase(){
    RealD mass=-1.8;   
    //Use same DSDR twists as https://arxiv.org/pdf/1208.4412.pdf
    RealD epsilon_f = 0.02; //numerator (in determinant)
    RealD epsilon_b = 0.5; 

    std::vector<int> twists(Nd,0);
    twists[nu] = 1; //GPBC in y
    twists[3] = 1; //APBC
    GparityAction::ImplParams params_2f;  params_2f.twists = twists;
    action_2f = new GparityWilsonTMFermionD(Umu_2f,*UGrid_2f,*UrbGrid_2f, mass, epsilon_f, params_2f);
    action_PV_2f = new GparityWilsonTMFermionD(Umu_2f,*UGrid_2f,*UrbGrid_2f, mass, epsilon_b, params_2f);

    DomainWallFermionD::ImplParams params_1f;  
    params_1f.boundary_phases[nu] = -1; 
    params_1f.boundary_phases[3] = -1; 

    action_1f = new WilsonTMFermionD(Umu_1f,*UGrid_1f,*UrbGrid_1f, mass, epsilon_f, params_1f);
    action_PV_1f = new WilsonTMFermionD(Umu_1f,*UGrid_1f,*UrbGrid_1f, mass, epsilon_b, params_1f);

    setupPseudofermion();
  }

  static bool is4d(){ return true; }
};


template<>
struct setupAction<GparityDomainWallFermionD, DomainWallFermionD>: public RatioActionSetupBase<GparityDomainWallFermionD, DomainWallFermionD>{
  typedef GparityDomainWallFermionD GparityAction;
  typedef DomainWallFermionD StandardAction;
  
  setupAction(GridCartesian* UGrid_2f, GridRedBlackCartesian* UrbGrid_2f,  GridCartesian* FGrid_2f, GridRedBlackCartesian* FrbGrid_2f,
	      GridCartesian* UGrid_1f, GridRedBlackCartesian* UrbGrid_1f,  GridCartesian* FGrid_1f, GridRedBlackCartesian* FrbGrid_1f,
	      LatticeGaugeField &Umu_2f, LatticeGaugeField &Umu_1f, int nu): RatioActionSetupBase(){
    RealD mass=0.01;   
    RealD M5=1.8; 

    std::vector<int> twists(Nd,0);
    twists[nu] = 1; //GPBC in y
    twists[3] = 1; //APBC
    GparityDomainWallFermionD::ImplParams params_2f;  params_2f.twists = twists;
    action_2f = new GparityDomainWallFermionD(Umu_2f,*FGrid_2f,*FrbGrid_2f,*UGrid_2f,*UrbGrid_2f,mass,M5,params_2f);
    action_PV_2f = new GparityDomainWallFermionD(Umu_2f,*FGrid_2f,*FrbGrid_2f,*UGrid_2f,*UrbGrid_2f,1.0,M5,params_2f);

    DomainWallFermionD::ImplParams params_1f;  
    params_1f.boundary_phases[nu] = -1; 
    params_1f.boundary_phases[3] = -1; 

    action_1f = new DomainWallFermionD(Umu_1f,*FGrid_1f,*FrbGrid_1f,*UGrid_1f,*UrbGrid_1f,mass,M5,params_1f);
    action_PV_1f = new DomainWallFermionD(Umu_1f,*FGrid_1f,*FrbGrid_1f,*UGrid_1f,*UrbGrid_1f,1.0,M5,params_1f);

    setupPseudofermion();
  }

  static bool is4d(){ return false; }
};





//For EOFA we need a different pseudofermion type
template<>
struct setupAction<GparityDomainWallEOFAFermionD, DomainWallEOFAFermionD>{
  typedef GparityDomainWallEOFAFermionD GparityAction;
  typedef DomainWallEOFAFermionD StandardAction;

  ExactOneFlavourRatioPseudoFermionAction<WilsonImplD> *pf_1f;
  ExactOneFlavourRatioPseudoFermionAction<GparityWilsonImplD> *pf_2f;

  GparityAction* action_2f;
  GparityAction* action_PV_2f;
  StandardAction* action_1f;
  StandardAction* action_PV_1f;

  ConjugateGradient<typename StandardAction::FermionField> CG_1f;
  ConjugateGradient<typename GparityAction::FermionField> CG_2f;

public:
  GparityAction & action2f(){ return *action_2f; }
  StandardAction & action1f(){ return *action_1f; }

  void refreshAction(LatticeGaugeField &Umu_2f, typename GparityAction::FermionField &eta_2f,
		     LatticeGaugeField &Umu_1f, typename StandardAction::FermionField &eta_1f){  
    pf_1f->refresh(Umu_1f, eta_1f);
    pf_2f->refresh(Umu_2f, eta_2f);

    //Compare PhiOdd
    RealD norm_1f = norm2(pf_1f->getPhi());
    RealD norm_2f = norm2(pf_2f->getPhi());
    
    std::cout << "Test Phi 2f: " << norm_2f << " 1f: " << norm_1f << std::endl;
  }

  void computeAction(RealD &S_2f, RealD &S_1f, LatticeGaugeField &Umu_2f, LatticeGaugeField &Umu_1f){
    S_1f = pf_1f->S(Umu_1f);
    S_2f = pf_2f->S(Umu_2f);
  }

  void computeDeriv(LatticeGaugeField &deriv_2f, LatticeGaugeField &deriv_1f, LatticeGaugeField &Umu_2f, LatticeGaugeField &Umu_1f){    
    pf_1f->deriv(Umu_1f, deriv_1f);
    pf_2f->deriv(Umu_2f, deriv_2f);
  }


  setupAction(GridCartesian* UGrid_2f, GridRedBlackCartesian* UrbGrid_2f,  GridCartesian* FGrid_2f, GridRedBlackCartesian* FrbGrid_2f,
	      GridCartesian* UGrid_1f, GridRedBlackCartesian* UrbGrid_1f,  GridCartesian* FGrid_1f, GridRedBlackCartesian* FrbGrid_1f,
	      LatticeGaugeField &Umu_2f, LatticeGaugeField &Umu_1f, int nu): CG_1f(1.0e-8,10000), CG_2f(1.0e-8,10000){
    RealD mass=0.01;   
    RealD M5=1.8; 

    std::vector<int> twists(Nd,0);
    twists[nu] = 1; //GPBC in y
    twists[3] = 1; //APBC
    GparityAction::ImplParams params_2f;  params_2f.twists = twists;
    action_2f = new GparityAction(Umu_2f,*FGrid_2f,*FrbGrid_2f,*UGrid_2f,*UrbGrid_2f, mass, mass, 1.0, 0.0, -1, M5, params_2f);
    action_PV_2f = new GparityAction(Umu_2f,*FGrid_2f,*FrbGrid_2f,*UGrid_2f,*UrbGrid_2f, 1.0, mass, 1.0, -1.0, 1, M5, params_2f); //cf Test_dwf_gpforce_eofa.cc

    StandardAction::ImplParams params_1f;  
    params_1f.boundary_phases[nu] = -1; 
    params_1f.boundary_phases[3] = -1; 

    action_1f = new StandardAction(Umu_1f,*FGrid_1f,*FrbGrid_1f,*UGrid_1f,*UrbGrid_1f, mass, mass, 1.0, 0.0, -1, M5, params_1f);
    action_PV_1f = new StandardAction(Umu_1f,*FGrid_1f,*FrbGrid_1f,*UGrid_1f,*UrbGrid_1f, 1.0, mass, 1.0, -1.0, 1, M5, params_1f);

    OneFlavourRationalParams RationalParams(0.95, 100.0, 5000, 1.0e-12, 12);

    pf_1f = new ExactOneFlavourRatioPseudoFermionAction<WilsonImplD>(*action_1f, *action_PV_1f, CG_1f, CG_1f, CG_1f, CG_1f, CG_1f, RationalParams, true);
    pf_2f = new ExactOneFlavourRatioPseudoFermionAction<GparityWilsonImplD>(*action_2f, *action_PV_2f, CG_2f, CG_2f, CG_2f, CG_2f, CG_2f, RationalParams, true);
  }

  static bool is4d(){ return false; }
};


template<typename GparityAction, typename StandardAction>
void runTest(int argc, char** argv){
  Grid_init(&argc,&argv);

  const int nu = 1;
  Coordinate latt_2f   = GridDefaultLatt();
  Coordinate latt_1f   = latt_2f;
  latt_1f[nu] *= 2;

  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  const int Ls=8;

  GridCartesian         * UGrid_1f   = SpaceTimeGrid::makeFourDimGrid(latt_1f, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_1f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_1f);
  GridCartesian         * FGrid_1f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_1f);
  GridRedBlackCartesian * FrbGrid_1f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_1f);


  GridCartesian         * UGrid_2f   = SpaceTimeGrid::makeFourDimGrid(latt_2f, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_2f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_2f);
  GridCartesian         * FGrid_2f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_2f);
  GridRedBlackCartesian * FrbGrid_2f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_2f);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG5_2f(FGrid_2f);  RNG5_2f.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4_2f(UGrid_2f);  RNG4_2f.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu_2f(UGrid_2f);
  SU<Nc>::HotConfiguration(RNG4_2f,Umu_2f);

  LatticeGaugeField Umu_1f(UGrid_1f);
  copyConjGauge(Umu_1f, Umu_2f, nu);

  typedef typename GparityAction::FermionField GparityFermionField;
  typedef typename StandardAction::FermionField StandardFermionField;

  setupAction<GparityAction, StandardAction> setup(UGrid_2f, UrbGrid_2f, FGrid_2f, FrbGrid_2f,
						   UGrid_1f, UrbGrid_1f, FGrid_1f, FrbGrid_1f,
						   Umu_2f, Umu_1f, nu);
  GridBase* FGrid_2f_a = setup.action2f().FermionGrid();
  GridBase* FGrid_1f_a = setup.action1f().FermionGrid();
  GridBase* FrbGrid_2f_a = setup.action2f().FermionRedBlackGrid();
  GridBase* FrbGrid_1f_a = setup.action1f().FermionRedBlackGrid();
  bool is_4d = setup.is4d();

  //Check components by doing an inversion
  {
    setup.action2f().ImportGauge(Umu_2f);
    setup.action1f().ImportGauge(Umu_1f);

    GparityFermionField src_2f(FGrid_2f_a);
    gaussian(is_4d ? RNG4_2f : RNG5_2f, src_2f);
    
    StandardFermionField src_1f(FGrid_1f_a);
    convertFermion1f_from_2f(src_1f, src_2f, nu, is_4d);

    StandardFermionField src_o_1f(FrbGrid_1f_a);
    StandardFermionField result_o_1f(FrbGrid_1f_a);
    pickCheckerboard(Odd,src_o_1f,src_1f);
    result_o_1f=Zero();

    SchurDiagMooeeOperator<StandardAction,StandardFermionField> HermOpEO_1f(setup.action1f());
    ConjugateGradient<StandardFermionField> CG_1f(1.0e-8,10000);
    CG_1f(HermOpEO_1f,src_o_1f,result_o_1f);


    GparityFermionField src_o_2f(FrbGrid_2f_a);
    GparityFermionField result_o_2f(FrbGrid_2f_a);
    pickCheckerboard(Odd,src_o_2f,src_2f);
    result_o_2f=Zero();

    SchurDiagMooeeOperator<GparityAction,GparityFermionField> HermOpEO_2f(setup.action2f());
    ConjugateGradient<GparityFermionField> CG_2f(1.0e-8,10000);
    CG_2f(HermOpEO_2f,src_o_2f,result_o_2f);

    RealD norm_1f = norm2(result_o_1f);
    RealD norm_2f = norm2(result_o_2f);

    std::cout << "Test fermion inversion 2f: " << norm_2f << " 1f: " << norm_1f << std::endl;
  }

  //Generate eta
  RealD scale = std::sqrt(0.5);

  GparityFermionField eta_2f(FGrid_2f_a);    
  gaussian(is_4d ? RNG4_2f : RNG5_2f,eta_2f); eta_2f = eta_2f * scale;

  StandardFermionField eta_1f(FGrid_1f_a);    
  convertFermion1f_from_2f(eta_1f, eta_2f, nu, is_4d);
  
  setup.refreshAction(Umu_2f, eta_2f, Umu_1f, eta_1f);
 
  //Initial action is just |eta^2|
  RealD S_1f, S_2f;

  setup.computeAction(S_2f, S_1f, Umu_2f, Umu_1f);

  std::cout << "Test Initial action 2f: " << S_2f << " 1f: " << S_1f << " diff: " << S_2f - S_1f << std::endl;

  //Do a random gauge field refresh
  SU<Nc>::HotConfiguration(RNG4_2f,Umu_2f);
  copyConjGauge(Umu_1f, Umu_2f, nu);

  //Compute the action again
  setup.computeAction(S_2f, S_1f, Umu_2f, Umu_1f);
  
  std::cout << "Test Action after gauge field randomize 2f: " << S_2f << " 1f: " << S_1f << " diff: " << S_2f - S_1f << std::endl;

  //Compute the derivative and test the conjugate relation
  LatticeGaugeField deriv_2f(UGrid_2f);
  LatticeGaugeField deriv_1f(UGrid_1f);
  setup.computeDeriv(deriv_2f, deriv_1f, Umu_2f, Umu_1f);

  //Have to combine the two forces on the 1f by symmetrizing under the complex conjugate
  {
    RealD norm2_pre = norm2(deriv_1f);
    LatticeGaugeField deriv_1f_shift = conjugate( Cshift(deriv_1f, nu, latt_2f[nu]) );
    deriv_1f = deriv_1f + deriv_1f_shift;
    std::cout << "Test combine/symmetrize forces on 1f lattice, dS/dU : " << norm2_pre << " -> " << norm2(deriv_1f) << std::endl;
  }
  
  LatticeGaugeField deriv_1f_from_2f(UGrid_1f);  
  copyConjGauge(deriv_1f_from_2f, deriv_2f, nu);
  std::cout << "Test copy-conj 2f dS/dU to obtain equivalent 1f force : " << norm2(deriv_2f) << " -> " << norm2(deriv_1f_from_2f) << std::endl;
  
  LatticeGaugeField diff_deriv_1f = deriv_1f - deriv_1f_from_2f;

  std::cout << "Test dS/dU 1f constructed from 2f derivative: " << norm2(deriv_1f_from_2f) << "  dS/dU 1f actual: " << norm2(deriv_1f) << "  Norm of difference: " << norm2(diff_deriv_1f) << std::endl;

  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}


  

int main (int argc, char ** argv)
{
  std::string action = "DWF";
  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--action"){
      action = argv[i+1];
    }
  }

  if(action == "DWF"){
    runTest<GparityDomainWallFermionD, DomainWallFermionD>(argc, argv);
  }else if(action == "EOFA"){
    runTest<GparityDomainWallEOFAFermionD, DomainWallEOFAFermionD>(argc, argv);
  }else if(action == "DSDR"){
    runTest<GparityWilsonTMFermionD, WilsonTMFermionD>(argc,argv);
  }else{
    assert(0);
  }
}
